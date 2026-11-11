#ifndef LOBPCG_H
#define LOBPCG_H

#include <functional>
#include <memory>
#include <sys/time.h>
#include <stdexcept>

#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "cg_method.h"
#include "dense_eigen.h"

namespace lobpcg_detail {

struct LinearOperator {
    std::function<void(VectorSparseG&, VectorSparseG&)> apply;
};

struct Preconditioner {
    std::function<void(VectorSparseG&, double)> apply;
};

bool run_lobpcg(VectorSparseG& x, double eps, int max_iterations,
                int block_size, double orthogonality_tolerance,
                const LinearOperator& A, const LinearOperator& B,
                const Preconditioner& precondition,
                const SmallGeneralizedEigenSolver& solver,
                const std::function<void(int, double, double)>& observer,
                double& eigenvalue, int& iterations);

inline double elapsed_seconds(const timeval& begin, const timeval& end) {
    return static_cast<double>(end.tv_sec - begin.tv_sec)
         + static_cast<double>(end.tv_usec - begin.tv_usec) * 1e-6;
}

} // namespace lobpcg_detail

struct LobpcgOptions {
    int block_size = 1;
    int max_iterations = 500;
    double orthogonality_tolerance = 1e-12;
    bool use_preconditioner = true;
    std::function<void(int, double, double)> observer;
};

class Lobpcg {
public:
    using Options = LobpcgOptions;

    template<class Stencil_left, class Stencil_right>
    static bool solve(double eps, double cg_eps, VectorSparseG& x,
                      double& eigenvalue, int& iterations,
                      double& avg_cg_iterations, MatrixVectorHomogen& matrix,
                      Stencil_left lhs, Stencil_right rhs,
                      bool rightoperator, double& precondition,
                      double& duration,
                      const Options& options = Options{});
};

template<class Stencil_left, class Stencil_right>
bool Lobpcg::solve(double eps, double cg_eps, VectorSparseG& x,
                   double& eigenvalue, int& iterations,
                   double& avg_cg_iterations, MatrixVectorHomogen& matrix,
                   Stencil_left lhs, Stencil_right rhs,
                   bool rightoperator, double& precondition,
                   double& duration, const Options& options) {
    (void)cg_eps;

    if(options.block_size != 1) {
        throw std::invalid_argument("LOBPCG: Blocksize != 1 is currently not supported");
    }

    AdaptiveSparseGrid* grid = x.getSparseGrid();

    VectorSparseG z_left(grid);
    VectorSparseG z_right(grid);

    std::unique_ptr<Preconditioning<Stencil_left>> left_preconditioner;
    std::unique_ptr<Preconditioning<Stencil_right>> right_preconditioner;

    timeval begin_precondition;
    gettimeofday(&begin_precondition, 0);

    if (options.use_preconditioner) {
        left_preconditioner.reset(
            new Preconditioning<Stencil_left>(z_left, lhs));
        if (rightoperator) {
            right_preconditioner.reset(
                new Preconditioning<Stencil_right>(z_right, rhs, true));
        }
    }

    timeval end_precondition;
    gettimeofday(&end_precondition, 0);
    precondition = lobpcg_detail::elapsed_seconds(begin_precondition,
                                                  end_precondition);

    lobpcg_detail::LinearOperator A{
        [&](VectorSparseG& in, VectorSparseG& out) {
            matrix.multiplication<Stencil_left>(in, out, lhs);
        }
    };

    lobpcg_detail::LinearOperator B{
        [&](VectorSparseG& in, VectorSparseG& out) {
            if (rightoperator) {
                matrix.multiplication<Stencil_right>(in, out, rhs);
            } else {
                out = in;
            }
        }
    };

    lobpcg_detail::Preconditioner pre{
        [&](VectorSparseG& residual, double lambda) {
            if (!left_preconditioner) {
                return;
            }
            if (right_preconditioner) {
                left_preconditioner->apply_inverse_shift(
                    &residual, *right_preconditioner, lambda);
            } else {
                left_preconditioner->apply_inverse(&residual);
            }
        }
    };

    lobpcg_detail::JacobiGeneralizedEigenSolver solver;

    timeval begin_solve;
    gettimeofday(&begin_solve, 0);

    const bool converged = lobpcg_detail::run_lobpcg(
        x, eps, options.max_iterations, options.block_size,
        options.orthogonality_tolerance, A, B, pre, solver,
        options.observer, eigenvalue, iterations);

    timeval end_solve;
    gettimeofday(&end_solve, 0);
    duration = lobpcg_detail::elapsed_seconds(begin_solve, end_solve);

    avg_cg_iterations = 0.0;

    return converged;
}

#endif // LOBPCG_H
