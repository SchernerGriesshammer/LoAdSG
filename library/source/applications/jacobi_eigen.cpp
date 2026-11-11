#include "dense_eigen.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace lobpcg_detail {
namespace {

bool cholesky_factor(const DenseMatrix& T, int n, double jitter,
                     std::vector<double>& L) {
    L.assign(static_cast<std::size_t>(n) * n, 0.0);

    for (int j = 0; j < n; ++j) {
        double diagonal = T(j, j) + jitter;
        for (int k = 0; k < j; ++k) {
            diagonal -= L[j * n + k] * L[j * n + k];
        }

        if (!(diagonal > 0.0)) {
            return false;
        }

        const double root = std::sqrt(diagonal);
        L[j * n + j] = root;

        for (int i = j + 1; i < n; ++i) {
            double value = T(i, j);
            for (int k = 0; k < j; ++k) {
                value -= L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = value / root;
        }
    }

    return true;
}

void forward_solve(const std::vector<double>& L, int n,
                   const double* rhs, double* x) {
    for (int i = 0; i < n; ++i) {
        double value = rhs[i];
        for (int k = 0; k < i; ++k) {
            value -= L[i * n + k] * x[k];
        }
        x[i] = value / L[i * n + i];
    }
}

void backward_solve_transposed(const std::vector<double>& L, int n,
                               const double* rhs, double* x) {
    for (int i = n - 1; i >= 0; --i) {
        double value = rhs[i];
        for (int k = i + 1; k < n; ++k) {
            value -= L[k * n + i] * x[k];
        }
        x[i] = value / L[i * n + i];
    }
}

void transform_to_standard(const DenseMatrix& S, const std::vector<double>& L,
                           int n, DenseMatrix& C) {
    std::vector<double> X(static_cast<std::size_t>(n) * n, 0.0);
    std::vector<double> rhs(n, 0.0);
    std::vector<double> solution(n, 0.0);

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            rhs[i] = S(i, j);
        }
        forward_solve(L, n, rhs.data(), solution.data());
        for (int i = 0; i < n; ++i) {
            X[i * n + j] = solution[i];
        }
    }

    C.resize(n, n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            rhs[i] = X[j * n + i];
        }
        forward_solve(L, n, rhs.data(), solution.data());
        for (int i = 0; i < n; ++i) {
            C(i, j) = solution[i];
        }
    }
}

void jacobi_rotations(DenseMatrix& C, DenseMatrix& Z, int n) {
    Z.resize(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Z(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    constexpr int max_sweeps = 100;
    constexpr double tolerance = 1e-15;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double off_diagonal = 0.0;
        double frobenius = 0.0;

        for (int i = 0; i < n; ++i) {
            frobenius += C(i, i) * C(i, i);
            for (int j = i + 1; j < n; ++j) {
                const double value = C(i, j);
                off_diagonal += value * value;
                frobenius += 2.0 * value * value;
            }
        }

        if (off_diagonal <= tolerance * tolerance * frobenius) {
            break;
        }

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const double apq = C(p, q);
                if (std::abs(apq) < 1e-300) {
                    continue;
                }

                const double app = C(p, p);
                const double aqq = C(q, q);
                const double tau = (aqq - app) / (2.0 * apq);
                const double sign = (tau >= 0.0) ? 1.0 : -1.0;
                const double t = sign / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = t * c;

                for (int k = 0; k < n; ++k) {
                    if (k == p || k == q) {
                        continue;
                    }
                    const double ckp = C(k, p);
                    const double ckq = C(k, q);
                    C(k, p) = c * ckp - s * ckq;
                    C(k, q) = s * ckp + c * ckq;
                    C(p, k) = C(k, p);
                    C(q, k) = C(k, q);
                }

                C(p, p) = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                C(q, q) = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                C(p, q) = 0.0;
                C(q, p) = 0.0;

                for (int k = 0; k < n; ++k) {
                    const double zkp = Z(k, p);
                    const double zkq = Z(k, q);
                    Z(k, p) = c * zkp - s * zkq;
                    Z(k, q) = s * zkp + c * zkq;
                }
            }
        }
    }
}

void back_transform(const std::vector<double>& L, const DenseMatrix& Z,
                    int n, DenseMatrix& V) {
    V.resize(n, n);
    std::vector<double> rhs(n, 0.0);
    std::vector<double> solution(n, 0.0);

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            rhs[i] = Z(i, j);
        }
        backward_solve_transposed(L, n, rhs.data(), solution.data());
        for (int i = 0; i < n; ++i) {
            V(i, j) = solution[i];
        }
    }
}

} // namespace

void JacobiGeneralizedEigenSolver::solve(const DenseMatrix& S,
                                         const DenseMatrix& T,
                                         std::vector<double>& eigenvalues,
                                         DenseMatrix& eigenvectors) const {

    if(S.rows() != S.cols() || S.rows() != T.rows() || T.rows() != T.cols()) {
        throw std::invalid_argument(
            "JacobiGeneralizedEigenSolver: Matrices have wrong dimensions");
    }
                                       
    const int n = S.rows();

    if (n <= 0) {
        eigenvalues.clear();
        eigenvectors.resize(0, 0);
        return;
    }

    double trace = 0.0;
    for (int i = 0; i < n; ++i) {
        trace += T(i, i);
    }

    std::vector<double> L;
    double jitter = 0.0;
    bool factorized = false;

    for (int attempt = 0; attempt < 4 && !factorized; ++attempt) {
        factorized = cholesky_factor(T, n, jitter, L);
        if (!factorized) {
            const double scale = trace / n + 1.0;
            jitter = (jitter == 0.0) ? 1e-12 * scale : jitter * 10.0;
        }
    }

    if (!factorized) {
        throw std::runtime_error(
            "JacobiGeneralizedEigenSolver: T is not positive definite");
    }

    DenseMatrix C;
    transform_to_standard(S, L, n, C);

    DenseMatrix Z;
    jacobi_rotations(C, Z, n);

    DenseMatrix V;
    back_transform(L, Z, n, V);

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&C](int a, int b) { return C(a, a) < C(b, b); });

    eigenvalues.resize(n);
    eigenvectors.resize(n, n);

    for (int j = 0; j < n; ++j) {
        eigenvalues[j] = C(order[j], order[j]);
        for (int i = 0; i < n; ++i) {
            eigenvectors(i, j) = V(i, order[j]);
        }
    }
}

} // namespace lobpcg_detail
