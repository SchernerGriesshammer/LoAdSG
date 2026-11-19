//
// ============================================================================
//  3D Laplace eigenvalue problem on a regular sparse grid -- a worked example
// ============================================================================
//
// This file is the 3D counterpart of mainfiles/laplace_eigenvalue.cpp (the
// plain 2D example). It is modelled after
// mainfiles/schroedinger_main/solve_schroedinger_hydro_simple.cpp, but with
// the variable coefficient (the "-1/r" potential of the hydrogen atom) taken
// out completely. What remains is the simplest possible eigenvalue problem
// you can pose on this library:
//
//      -laplacian u = lambda * u      in  (0,1)^3
//                 u = 0               on  the boundary of the cube
//
// i.e. we are looking for the eigenpairs of the (negative) Laplace operator
// with homogeneous Dirichlet boundary conditions.
//
// -----------------------------------------------------------------------------
// From the PDE to a matrix eigenvalue problem
// -----------------------------------------------------------------------------
// The library discretises functions on (0,1)^d in a prewavelet basis on a
// (locally refined / here: regular) sparse grid. Testing the weak form of
// -laplacian u = lambda u with every prewavelet basis function phi_p gives
//
//      sum_q  a(phi_q, phi_p) * x_q   =   lambda * sum_q  m(phi_q, phi_p) * x_q
//
// for every basis function p, where x_q are the (unknown) coefficients of u
// in the prewavelet basis, a(.,.) is the bilinear form of the Laplacian
// (stiffness), and m(.,.) is the L2 scalar product (mass). Collecting all p
// into a vector gives the generalised matrix eigenvalue problem
//
//      A x = lambda * M x                                                (*)
//
// In this library:
//   * the stiffness matrix A is assembled on the fly by the "Poisson" stencil
//     class (source/stencils/Stencil.h) -- despite the generic name, this is
//     exactly the discrete Laplacian / stiffness matrix of -laplacian.
//   * the mass matrix M is assembled on the fly by the "HelmHoltz" stencil
//     class. Despite the name (historically this class was written for a
//     Helmholtz-type equation), for our purposes it is simply the L2 mass
//     matrix, i.e. m(phi_q, phi_p) = integral phi_q * phi_p.
//
// Both classes implement the StencilTemplate interface, so a generic
// "MatrixVectorHomogen" object can apply either of them to a vector of
// prewavelet coefficients (A*x or M*x) without ever assembling A or M as an
// explicit (dense or sparse) matrix -- this is the "matrix-free" /
// "on-the-fly" philosophy of the whole code base.
//
// Because there is no variable coefficient here (unlike in the hydrogen
// example, where a "-1/r" singular potential had to be added to A), A is
// simply the (positive definite!) discrete Laplacian. This means we do NOT
// need any "shift trick" to make A positive definite -- it already is -- and
// we can go straight for the smallest eigenvalue with a plain inverse power
// iteration.
//
// -----------------------------------------------------------------------------
// Inverse power iteration
// -----------------------------------------------------------------------------
// To find the SMALLEST eigenvalue of (*) we use inverse iteration:
//
//      x_{k+1} = A^{-1} M x_k ,           x_{k+1} <- x_{k+1} / ||x_{k+1}||
//
// Starting from (almost) any x_0, x_k converges to the eigenvector belonging
// to the smallest eigenvalue lambda_min, and lambda_min = 1 / ||x_{k+1}||
// (with x_k already normalised to ||x_k|| = 1 in the previous step).
// Intuitively: applying A^{-1} amplifies the components of x_k that belong to
// small eigenvalues the most, so repeated application "filters out" every
// eigenvector except the one for the smallest eigenvalue.
//
// Every application of A^{-1} (i.e. every "solve A x_{k+1} = M x_k for
// x_{k+1}") is itself done by the preconditioned conjugate gradient method,
// since A is symmetric positive definite.
//
// -----------------------------------------------------------------------------
// Shifted inverse iteration
// -----------------------------------------------------------------------------
// The convergence speed of inverse iteration is governed by the ratio
// lambda_1 / lambda_2 of the two smallest eigenvalues. With a shift sigma we
// instead iterate with (A - sigma*M)^{-1} M, i.e.
//
//      solve (A - sigma*M) x_{k+1} = M x_k ,     x_{k+1} <- x_{k+1} / ||x_{k+1}||
//
// which converges with the ratio (lambda_1 - sigma) / (lambda_2 - sigma):
// the closer sigma is to lambda_1, the faster. Since the inner solves are
// done by CG, A - sigma*M has to stay positive definite, so sigma must lie
// below lambda_1. sigma = 0 is the plain inverse iteration from above.
//
// The whole procedure -- outer shifted inverse iteration, inner CG solves,
// convergence monitoring -- is implemented in
//
//      Power::power_inverse_shifted<Stencil_left, Stencil_right>(...)
//
// in source/applications/power_method.h. We call it with
// Stencil_left = Poisson (our A), Stencil_right = HelmHoltz (our M) and the
// shift chosen below in this file. The returned eigenvalue is the Rayleigh
// quotient x^T A x / x^T M x, i.e. directly lambda of (*), no back-shift
// needed.
//
// -----------------------------------------------------------------------------
// Known analytic solution (used to check convergence)
// -----------------------------------------------------------------------------
// On the unit cube with homogeneous Dirichlet boundary conditions, the
// eigenpairs of -laplacian are known in closed form:
//
//      lambda_{n1,n2,n3} = (n1^2 + n2^2 + n3^2) * pi^2 ,     n1,n2,n3 = 1,2,3,...
//      u_{n1,n2,n3}(x,y,z) = sin(n1*pi*x) * sin(n2*pi*y) * sin(n3*pi*z)
//
// The smallest one is the "ground state" n1=n2=n3=1:
//
//      lambda_1 = 3 * pi^2  ~=  29.6088
//      u_1(x,y,z) = sin(pi*x) * sin(pi*y) * sin(pi*z)
//
// We use this to measure both the error in the eigenvalue and the error in
// the eigenfunction as the sparse grid is refined.
// ============================================================================

#include "../source/sparseEXPDE.h"
#include "../source/iterator/depthIterator.h"
#include "../source/applications/power_method.h"

#include <string>
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>


// This file is written for the 3D case. The spatial dimension is a single,
// global compile-time constant (#define DimensionSparseGrid ... in
// source/abbrevi.h) that applies to the *whole* library, not just to one main
// file -- exactly like it is already done for the hydrogen (3D) and helium
// (6D) examples. Before building this file, make sure
//      source/abbrevi.h   contains   #define DimensionSparseGrid 3
// The static_assert below turns a forgotten edit into a clear compile error
// instead of a silently wrong (e.g. 2D) run.
static_assert(DimensionSparseGrid == 3,
              "This example is written for the 3D Laplace eigenvalue problem. "
              "Set '#define DimensionSparseGrid 3' in source/abbrevi.h and rebuild.");


/// @brief stopping criterion for the outer (power) iteration: the iteration
/// stops once two successive estimates of the eigenvalue/Rayleigh quotient
/// differ by less than eps (see Power::power_inverse_shifted for the exact test).
const double eps = 1e-10;

/// @brief stopping criterion for the inner CG solves (A x = M x_old). This
/// can be fairly tight here because A is well conditioned (plain Laplacian,
/// no singular coefficient), unlike in the hydrogen example.
double cg_eps = 1e-12;

/// @brief shift sigma of the inverse iteration, see the file header. It must
/// be smaller than the smallest discrete eigenvalue (~3*pi^2 ~= 29.6, the
/// discrete one is slightly larger); the next eigenvalue is 6*pi^2 ~= 59.2,
/// so the convergence ratio drops from ~0.5 (shift 0) to ~0.13 (shift 25).
/// Set it to 0.0 for the unshifted inverse iteration.
const double shift = 25.0;

/// @brief first and last refinement level of the regular sparse grid. Every
/// increase in level roughly doubles the resolution in every coordinate
/// direction while only increasing the number of degrees of freedom (DOFs)
/// log-linearly -- this is the whole point of using a *sparse* grid instead
/// of a full tensor grid, which would need 2^(level*DimensionSparseGrid) DOFs.
const int level_start = 2;
const int level_end = 9;

/// @brief exact smallest eigenvalue 3*pi^2 of -laplacian on the unit cube
/// with homogeneous Dirichlet boundary conditions (ground state n1=n2=n3=1).
const double eig_exp = double(DimensionSparseGrid) * M_PI * M_PI;


/// @brief mpi appropriate output: only rank 0 prints, every other rank stays
/// silent. In a purely serial (non-MPI) build, rank is always 0, so this
/// simply prints unconditionally.
/// @param s output string
/// @param endline whether to end with a linebreak
void mpi_cout(string s, bool endline = true) {
    int rank = 0;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == 0) {
        cout << s;
        if (endline) cout << endl;
        cout << flush;
    }
}


/// @brief rank of this process as assigned by the MPI launcher (mpirun/srun),
/// read from the environment. Used when the library is built WITHOUT MPI
/// (MY_MPI_ON undefined): then MPI_Comm_rank is unavailable and every process
/// started by mpirun would otherwise believe it is rank 0 and print.
/// @return launcher rank, or 0 when not started by an MPI launcher
int launcher_rank() {
    for (const char *var : {"OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK",
                            "MV2_COMM_WORLD_RANK", "SLURM_PROCID"}) {
        if (const char *val = getenv(var)) return atoi(val);
    }
    return 0;
}


/// @brief the exact ground-state eigenfunction sin(pi*x)*sin(pi*y)*sin(pi*z),
/// evaluated at a grid point. Note that I.coordinate(d) already returns the
/// coordinate in (0,1), so no coordinate transform is needed here (unlike in
/// the hydrogen example, where unitSQ2atom() rescales (0,1) to (-border,border)).
/// @param I grid point
/// @return value of the exact eigenfunction at I
double known_eigenfunction(IndexDimension I) {
    double val = 1.0;
    for (int d = 0; d < DimensionSparseGrid; d++) {
        val *= sin(M_PI * I.coordinate(d));
    }
    return val;
}


int main(int argc, char **argv) {

    int rank = 0;
    int num_tasks = 1;

#ifdef MY_MPI_ON
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#else
    rank = launcher_rank();
#endif

    // Only rank 0 writes to the screen. Detaching cout's buffer on all other
    // ranks also silences output printed from inside the library (e.g. the
    // per-iteration line of Power::power_inverse_shifted).
    if (rank != 0) cout.rdbuf(nullptr);

    // Standard OpenMP boilerplate: just report how many threads we got, the
    // actual OpenMP parallelism is hidden inside the stencil / matrix-vector
    // routines of the library.
#pragma omp parallel
    {
#pragma omp single
        {
            int total_threads = omp_get_num_threads();
            mpi_cout("number of mpi_tasks " + to_string(num_tasks));
            mpi_cout("number of omp_threads " + to_string(total_threads));
        }
    }

    if (rank == 0) {
        cout << endl;
        cout << "Computing the smallest eigenvalue of -laplacian u = lambda u" << endl;
        cout << "on the unit cube (0,1)^" << DimensionSparseGrid
             << " with homogeneous Dirichlet boundary conditions." << endl;
        cout << "No variable coefficient, no coordinate transform. "
             << "Shifted inverse iteration with (A - shift*M), shift < smallest eigenvalue." << endl;
        cout << "The expected (smallest) eigenvalue is " << DimensionSparseGrid
             << "*pi^2 = " << setprecision(12) << eig_exp << setprecision(6) << endl;
        cout << "With epsilon (outer, power iteration) = " << eps << endl;
        cout << "With cg_epsilon (inner, CG solves)     = " << cg_eps << endl;
        cout << "With shift                             = " << shift << endl;
        cout << "level start " << level_start << endl;
        cout << "level end   " << level_end << endl;
        cout << endl;
        cout << "level\tDOFS\t\teig\t\t\terr\t\trate\terr_rel\t\tL2(u)\t\trate\tpow_it\tcg_it" << endl;
        cout << endl;
    }

    double err_old = 0.0;
    double l2_old = 0.0;

    for (int level = level_start; level < level_end + 1; level++) {

        // -------------------------------------------------------------
        // 1) Build a *regular* sparse grid of the given level.
        //    "Regular" here means: every grid point allowed by the sparse
        //    grid construction up to this level is present (no adaptive
        //    refinement based on an error indicator) -- the simplest choice,
        //    and exactly what solve_schroedinger_hydro_simple.cpp does too.
        // -------------------------------------------------------------
        AdaptiveSparseGrid grid;
        IndexDimension centerPoint; // the default-constructed IndexDimension
                                    // is the single point in the exact center
                                    // of the domain, i.e. (0.5,0.5,0.5).
        grid.AddRecursiveSonsOfPoint(centerPoint, level);

        // MultiLevelAdaptiveSparseGrid groups the grid points of "grid" into
        // depth-ordered full subgrids; this bookkeeping structure is what
        // lets MatrixVectorHomogen apply the stencils efficiently.
        MultiLevelAdaptiveSparseGrid mgrid(&grid);

        // "symmetry = true" tells MatrixVectorHomogen that our bilinear forms
        // are symmetric (true for both Poisson and HelmHoltz), which enables
        // an optimisation that only computes "half" of the stencil entries.
        //
        // The last argument switches on MPI-based distribution of the
        // stencil work across ranks. We only want that if we really run with
        // more than one MPI rank -- turning it on unconditionally (as e.g.
        // solve_schroedinger_hydro_simple.cpp does with "MatrixVectorHomogen
        // m(grid, mgrid, 1, 1);") only works in a build that actually defines
        // MY_MPI_ON; a purely serial build would hit an MPI-only code path
        // and abort. Guarding it with "num_tasks > 1" keeps this example
        // runnable both as a serial test and as a real MPI job.
        MatrixVectorHomogen m(grid, mgrid, true, num_tasks > 1);

        // -------------------------------------------------------------
        // 2) The two stencils making up the generalised eigenvalue problem
        //    A x = lambda * M x. Both are assembled matrix-free, i.e. there
        //    is never an explicit matrix in memory -- "lhs"/"rhs" only know
        //    how to compute, for any two grid points, the local contribution
        //    of their bilinear form.
        // -------------------------------------------------------------
        Poisson lhs(grid);    // A: stiffness matrix of -laplacian
        HelmHoltz rhs(grid);  // M: mass matrix (L2 scalar product)

        // -------------------------------------------------------------
        // 3) A good initial guess for the power iteration.
        //    Instead of starting from an arbitrary vector (e.g. all ones),
        //    we cheat a little and start from the *exact* solution,
        //    evaluated nodally and then transformed into the prewavelet
        //    basis via calcPrewByNodal. This is exactly what
        //    solve_schroedinger_hydro_simple.cpp does with its
        //    known_eigenfunction() as well -- it makes the power iteration
        //    converge in very few steps and demonstrates the nodal <->
        //    prewavelet basis transform used throughout the library.
        //    (Power iteration would also converge from a generic starting
        //    vector, just more slowly.)
        // -------------------------------------------------------------
        VectorSparseG u_nodal_start(grid);
        VectorSparseG prew_u(grid);

        for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++) {
            IndexDimension I = grid.getIndexOfTable(i);
            u_nodal_start.setValue(i, known_eigenfunction(I));
        }
        calcPrewByNodal(prew_u, u_nodal_start);

        // -------------------------------------------------------------
        // 4) Solve the generalised eigenvalue problem for the smallest
        //    eigenvalue by shifted inverse power iteration (see the file
        //    header for an explanation of the algorithm).
        // -------------------------------------------------------------
        double eigenvalue_power = 0.0; // out: the computed smallest eigenvalue
        int pow_iterations = 0;        // out: number of outer (power) iterations
        double cg_iterations = 0.0;    // out: average number of inner CG iterations
        double time_power = 0.0;       // out: wall-clock time of the iteration

        Power::power_inverse_shifted<Poisson, HelmHoltz>(
                eps, cg_eps, shift, prew_u,
                eigenvalue_power, pow_iterations, cg_iterations,
                m, lhs, rhs,
                /*rightoperator=*/true,
                time_power);

        // -------------------------------------------------------------
        // 5) Compare against the known analytic solution.
        // -------------------------------------------------------------
        double err = abs(eigenvalue_power - eig_exp);
        double err_rel = err / eig_exp;

        // The eigenvector computed above is only determined up to an
        // arbitrary nonzero scalar factor (and sign): if u solves the
        // eigenvalue problem, so does c*u for any c != 0. power_inverse_shifted
        // normalises its result in its own (somewhat arbitrary) way, which
        // does not match the normalisation of known_eigenfunction(). To
        // compare the two pointwise, we rescale the computed eigenfunction
        // so that it matches the exact one *at the center point*
        // (0.5,0.5,0.5), where known_eigenfunction equals
        // sin(pi/2)^3 = 1. This fixes both the unknown scale and the sign.
        VectorSparseG u_nodal(grid), u_exact(grid), d(grid);
        u_nodal = 0.0;
        calcNodalByPrew(prew_u, u_nodal);

        for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++) {
            IndexDimension I = grid.getIndexOfTable(i);
            u_exact.setValue(i, known_eigenfunction(I));
        }

        double u_center = u_nodal.getValue(centerPoint);
        double scale = (abs(u_center) > 1e-30) ? known_eigenfunction(centerPoint) / u_center : 1.0;

        d = scale * u_nodal - u_exact;

        L2 L2_(grid, mgrid);
        double l2 = L2_.getValue(d);

        // convergence rate err(level-1)/err(level); undefined on the first level
        string err_rate = (err_old > 0.0 && err > 0.0) ? to_string(err_old / err).substr(0, 4) : "-";
        string l2_rate = (l2_old > 0.0 && l2 > 0.0) ? to_string(l2_old / l2).substr(0, 4) : "-";

        if (rank == 0) {
            cout << level
                 << "\t" << grid.getDOFS()
                 << "\t\t" << scientific << setprecision(12) << eigenvalue_power
                 << "\t" << setprecision(3) << err
                 << "\t" << err_rate
                 << "\t" << scientific << setprecision(3) << err_rel
                 << "\t" << l2
                 << "\t" << l2_rate
                 << "\t" << pow_iterations
                 << "\t" << fixed << setprecision(1) << cg_iterations
                 << defaultfloat << setprecision(6) << endl;
        }

        err_old = err;
        l2_old = l2;
    }

#ifdef MY_MPI_ON
    MPI_Finalize();
#endif
    return 0;
}
