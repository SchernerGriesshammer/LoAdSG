#include "lobpcg.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace lobpcg_detail {
namespace {

class BlockStore {
public:
    BlockStore(AdaptiveSparseGrid* grid, int capacity) : active_(0) {
        columns_.reserve(capacity);
        for (int i = 0; i < capacity; ++i) {
            columns_.push_back(std::unique_ptr<VectorSparseG>(
                new VectorSparseG(grid)));
        }
    }

    int capacity() const { return static_cast<int>(columns_.size()); }
    int active() const { return active_; }
    void set_active(int active) { active_ = active; }

    VectorSparseG& operator[](int i) { return *columns_[i]; }
    const VectorSparseG& operator[](int i) const { return *columns_[i]; }

    void swap(BlockStore& other) noexcept {
        columns_.swap(other.columns_);
        std::swap(active_, other.active_);
    }

private:
    std::vector<std::unique_ptr<VectorSparseG>> columns_;
    int active_;
};

bool b_orthonormalize(VectorSparseG& column,
                      const std::vector<VectorSparseG*>& basis,
                      const LinearOperator& B, VectorSparseG& B_column,
                      double tolerance) {
    const double threshold = tolerance * tolerance;

    for (int pass = 0; pass < 2; ++pass) {
        B.apply(column, B_column);
        for (std::size_t j = 0; j < basis.size(); ++j) {
            const double coefficient = product(*basis[j], B_column);
            if (coefficient != 0.0) {
                column = column - coefficient * (*basis[j]);
            }
        }
    }

    B.apply(column, B_column);
    const double norm2 = product(column, B_column);
    if (!(norm2 > threshold)) {
        return false;
    }

    const double scale = 1.0 / std::sqrt(norm2);
    column = column * scale;
    B_column = B_column * scale;
    return true;
}

void linear_combination(VectorSparseG& destination,
                        const std::vector<VectorSparseG*>& columns,
                        const DenseMatrix& coefficients, int column,
                        int first, int count) {
    destination = 0.0;
    for (int j = 0; j < count; ++j) {
        const double coefficient = coefficients(first + j, column);
        if (coefficient != 0.0) {
            destination = destination + coefficient * (*columns[first + j]);
        }
    }
}

void assemble_gram(const std::vector<VectorSparseG*>& S,
                   const std::vector<VectorSparseG*>& AS,
                   const std::vector<VectorSparseG*>& BS,
                   DenseMatrix& S_A, DenseMatrix& S_B) {
    const int n = static_cast<int>(S.size());
    S_A.resize(n, n);
    S_B.resize(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            const double a = product(*S[i], *AS[j]);
            const double b = product(*S[i], *BS[j]);
            S_A(i, j) = a;
            S_A(j, i) = a;
            S_B(i, j) = b;
            S_B(j, i) = b;
        }
    }
}

} // namespace

bool run_lobpcg(VectorSparseG& x, double eps, int max_iterations,
                int block_size, double orthogonality_tolerance,
                const LinearOperator& A, const LinearOperator& B,
                const Preconditioner& precondition,
                const SmallGeneralizedEigenSolver& solver,
                const std::function<void(int, double, double)>& observer,
                double& eigenvalue, int& iterations) {
    AdaptiveSparseGrid* grid = x.getSparseGrid();

    const int capacity = std::max(1, block_size);
    const int subspace_capacity = 3 * capacity;

    BlockStore X(grid, capacity);
    BlockStore P(grid, capacity);
    BlockStore W(grid, capacity);
    BlockStore AX(grid, capacity);
    BlockStore BX(grid, capacity);
    BlockStore AW(grid, capacity);
    BlockStore AP(grid, capacity);
    BlockStore BW(grid, capacity);
    BlockStore BP(grid, capacity);
    BlockStore Xnew(grid, capacity);
    BlockStore Pnew(grid, capacity);
    BlockStore AXnew(grid, capacity);
    BlockStore BXnew(grid, capacity);
    BlockStore APnew(grid, capacity);
    BlockStore BPnew(grid, capacity);

    std::vector<VectorSparseG*> S;
    std::vector<VectorSparseG*> AS;
    std::vector<VectorSparseG*> BS;
    std::vector<VectorSparseG*> basis;
    S.reserve(subspace_capacity);
    AS.reserve(subspace_capacity);
    BS.reserve(subspace_capacity);
    basis.reserve(3 * capacity);

    std::vector<double> lambdas;
    std::vector<double> values;
    std::vector<char> w_ok;
    lambdas.reserve(capacity);
    values.reserve(subspace_capacity);
    w_ok.reserve(capacity);

    DenseMatrix S_A(subspace_capacity, subspace_capacity);
    DenseMatrix S_B(subspace_capacity, subspace_capacity);
    DenseMatrix RitzVectors(subspace_capacity, subspace_capacity);

    X[0] = x;
    X.set_active(1);
    P.set_active(0);

    basis.clear();
    if (!b_orthonormalize(X[0], basis, B, BX[0],
                          orthogonality_tolerance)) {
        iterations = 0;
        eigenvalue = 0.0;
        return false;
    }

    A.apply(X[0], AX[0]);
    AX.set_active(1);
    BX.set_active(1);
    AP.set_active(0);
    BP.set_active(0);

    double lambda = 0.0;
    double lambda_previous = 0.0;
    bool converged = false;
    int iteration = 0;

    for (iteration = 0; iteration < max_iterations; ++iteration) {
        lambdas.resize(X.active());
        W.set_active(X.active());

        double residual = 0.0;

        for (int i = 0; i < X.active(); ++i) {
            const double denominator = product(X[i], BX[i]);
            if (std::abs(denominator) < 1e-300) {
                iterations = iteration;
                eigenvalue = lambda;
                return false;
            }

            lambdas[i] = product(X[i], AX[i]) / denominator;
            lambda = lambdas[i];
            W[i] = AX[i] - lambdas[i] * BX[i];

            const double w_norm = std::sqrt(product(W[i], W[i]));
            const double a_norm = std::sqrt(product(AX[i], AX[i]));
            const double relative = w_norm / std::max(a_norm, 1e-300);
            residual = std::max(residual, relative);
        }

        const bool eigenvalue_converged =
            std::abs(lambda - lambda_previous) < eps;
        lambda_previous = lambda;

        if (observer) {
            observer(iteration, lambda, residual);
        }

        if (eigenvalue_converged) {
            converged = true;
            break;
        }

        for (int i = 0; i < W.active(); ++i) {
            precondition.apply(W[i], lambdas[i]);
        }

        basis.clear();
        for (int i = 0; i < X.active(); ++i) {
            basis.push_back(&X[i]);
        }
        for (int i = 0; i < P.active(); ++i) {
            basis.push_back(&P[i]);
        }

        w_ok.resize(W.active());
        bool any_w = false;

        for (int i = 0; i < W.active(); ++i) {
            w_ok[i] = b_orthonormalize(W[i], basis, B, BW[i],
                                       orthogonality_tolerance) ? 1 : 0;
            if (w_ok[i]) {
                basis.push_back(&W[i]);
                A.apply(W[i], AW[i]);
                any_w = true;
            }
        }

        if (!any_w) {
            break;
        }

        S.clear();
        AS.clear();
        BS.clear();

        for (int i = 0; i < X.active(); ++i) {
            S.push_back(&X[i]);
            AS.push_back(&AX[i]);
            BS.push_back(&BX[i]);
        }
        for (int i = 0; i < W.active(); ++i) {
            if (w_ok[i]) {
                S.push_back(&W[i]);
                AS.push_back(&AW[i]);
                BS.push_back(&BW[i]);
            }
        }
        for (int i = 0; i < P.active(); ++i) {
            S.push_back(&P[i]);
            AS.push_back(&AP[i]);
            BS.push_back(&BP[i]);
        }

        const int subspace_size = static_cast<int>(S.size());

        assemble_gram(S, AS, BS, S_A, S_B);
        solver.solve(S_A, S_B, values, RitzVectors);

        const int x_count = X.active();
        const int wp_count = subspace_size - x_count;

        Xnew.set_active(1);
        AXnew.set_active(1);
        BXnew.set_active(1);
        linear_combination(Xnew[0], S, RitzVectors, 0, 0, subspace_size);
        linear_combination(AXnew[0], AS, RitzVectors, 0, 0, subspace_size);
        linear_combination(BXnew[0], BS, RitzVectors, 0, 0, subspace_size);

        const double x_norm2 = product(Xnew[0], BXnew[0]);
        if (!(x_norm2 > orthogonality_tolerance * orthogonality_tolerance)) {
            break;
        }

        const double x_scale = 1.0 / std::sqrt(x_norm2);
        Xnew[0] = Xnew[0] * x_scale;
        AXnew[0] = AXnew[0] * x_scale;
        BXnew[0] = BXnew[0] * x_scale;

        if (wp_count > 0) {
            Pnew.set_active(1);
            APnew.set_active(1);
            BPnew.set_active(1);
            linear_combination(Pnew[0], S, RitzVectors, 0,
                               x_count, wp_count);
            linear_combination(APnew[0], AS, RitzVectors, 0,
                               x_count, wp_count);
            linear_combination(BPnew[0], BS, RitzVectors, 0,
                               x_count, wp_count);

            const double p_norm2 = product(Pnew[0], BPnew[0]);
            if (p_norm2 > orthogonality_tolerance * orthogonality_tolerance) {
                const double p_scale = 1.0 / std::sqrt(p_norm2);
                Pnew[0] = Pnew[0] * p_scale;
                APnew[0] = APnew[0] * p_scale;
                BPnew[0] = BPnew[0] * p_scale;
            } else {
                Pnew.set_active(0);
                APnew.set_active(0);
                BPnew.set_active(0);
            }
        } else {
            Pnew.set_active(0);
            APnew.set_active(0);
            BPnew.set_active(0);
        }

        X.swap(Xnew);
        AX.swap(AXnew);
        BX.swap(BXnew);
        P.swap(Pnew);
        AP.swap(APnew);
        BP.swap(BPnew);

        lambda = values[0];
    }

    x = X[0];
    eigenvalue = lambda;
    iterations = iteration;

    return converged;
}

} // namespace lobpcg_detail
