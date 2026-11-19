//
// Simple Laplace eigenvalue problem on a regular sparse grid.
//
// Solve for the smallest eigenvalue of
//
//      -laplacian u = lambda * u   in  (0,1)^d
//                 u = 0            on  boundary
//
// In the (prewavelet) Galerkin discretisation this is the generalised
// eigenvalue problem
//
//      A x = lambda * M x
//
// with A the Poisson (stiffness) matrix and M the HelmHoltz (mass) matrix.
// The smallest eigenvalue is computed by inverse power iteration, where every
// inner system A x_new = M x_old is solved with the preconditioned CG method.
//
// Exact solution on the unit square, no coordinate transformation needed:
//      lambda_1 = d * pi^2       (= 2*pi^2 in 2D)
//      u_1(x)   = prod_d sin(pi * x_d)
//
// Structure follows mainfiles/schroedinger_main/solve_schroedinger_hydro_simple.cpp
// and solve_he_params.cpp, but without variable coefficient / Monte Carlo
// stencil and without the shift trick: the Laplacian is already positive
// definite and we are after its *smallest* eigenvalue.
//

#include "../source/sparseEXPDE.h"
#include "../source/iterator/depthIterator.h"
#include "../source/applications/power_method.h"

#include <string>
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <iomanip>


/// @brief stopping criterion for the power iteration (on the Rayleigh quotient)
const double eps = 1e-10;

/// @brief stopping criterion for the inner CG solves
double cg_eps = 1e-14;

/// @brief maximal number of power iterations
const int max_power_iterations = 200;

/// @brief first and last refinement level of the regular sparse grid
const int level_start = 2;
const int level_end = 10;

/// @brief exact smallest eigenvalue d * pi^2 on (0,1)^d
const double eig_exp = double(DimensionSparseGrid) * M_PI * M_PI;


/// @brief mpi appropriate output
/// @param s output string
/// @param endline whether to end with linebreak
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


/// @brief exact eigenfunction prod_d sin(pi*x_d), normalised to 1 in the center
/// @param I grid point
/// @return value at I
double known_eigenfunction(IndexDimension I) {
    double val = 1.0;
    for (int d = 0; d < DimensionSparseGrid; d++) {
        val *= sin(M_PI * I.coordinate(d));
    }
    return val;
}


/// @brief inverse power iteration for the smallest eigenvalue of A x = lambda M x
/// @param x in: start vector, out: eigenvector (prewavelet coefficients)
/// @param eigenvalue out: smallest eigenvalue (Rayleigh quotient)
/// @param iterations out: number of power iterations
/// @param avg_cg_iterations out: average number of CG iterations per power step
/// @param duration out: wall clock time of the iteration in seconds
/// @param matrix matrix vector multiplication object
/// @param lhs stencil of A
/// @param rhs stencil of M
template<class Stencil_left, class Stencil_right>
void power_inverse_laplace(VectorSparseG &x, double &eigenvalue, int &iterations,
                           double &avg_cg_iterations, double &duration,
                           MatrixVectorHomogen &matrix, Stencil_left lhs, Stencil_right rhs) {

    AdaptiveSparseGrid *grid = x.getSparseGrid();

    int rank = 0;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    VectorSparseG z(grid);
    VectorSparseG x_old(grid), x_new(grid), y(grid);
    VectorSparseG Ax(grid), Mx(grid);

    // Preconditioner for the inner CG solves with A
    Preconditioning<Stencil_left> preconditioning(z, lhs);

    struct timeval begin, end;
    gettimeofday(&begin, 0);

    x_old = x / sqrt(product(x, x));
    x_new = x_old; // < initial guess for the first CG solve

    int cg_iterations = 0;
    double cg_time = 0.0;
    double cg_iterations_sum = 0.0;

    double lambda = 0.0;
    double lambda_old = 0.0;

    int k = 0;
    for (int i = 0; i < max_power_iterations; i++) {
        k = i + 1;

        // y := M * x_old
        matrix.multiplication<Stencil_right>(x_old, y, rhs);

        // solve A * x_new = y  (x_new enters as initial guess)
        CG::solveHomogenP<Stencil_left>(cg_eps, x_new, y, &cg_iterations, x_old,
                                        matrix, lhs, &cg_time, preconditioning);
        cg_iterations_sum += cg_iterations;

        // Rayleigh quotient lambda = (x,Ax)/(x,Mx)
        matrix.multiplication<Stencil_left>(x_new, Ax, lhs);
        matrix.multiplication<Stencil_right>(x_new, Mx, rhs);
        lambda = product(x_new, Ax) / product(x_new, Mx);

        x_old = x_new / sqrt(product(x_new, x_new));

        double val = abs(lambda - lambda_old);
        lambda_old = lambda;

        if (rank == 0) {
            cout << "\tpower it " << setw(4) << k
                 << "\tlambda " << scientific << setprecision(12) << lambda
                 << "\tdiff " << val
                 << "\tcg_it " << cg_iterations
                 << "\tcg_eps " << cg_eps << defaultfloat << endl;
        }

        if (val < eps) break;
    }

    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    duration = seconds + microseconds * 1e-6;

    eigenvalue = lambda;
    iterations = k;
    avg_cg_iterations = cg_iterations_sum / double(k);
    x = x_new / sqrt(product(x_new, x_new));
}


int main(int argc, char **argv) {

    int rank = 0;
    int num_tasks = 1;

#ifdef MY_MPI_ON
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif

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
        cout << "Computing the smallest eigenvalue of -laplacian u = lambda u on (0,1)^"
             << DimensionSparseGrid << " with homogeneous Dirichlet boundary." << endl;
        cout << "Grid dimension " << DimensionSparseGrid << endl;
        cout << "Regular sparse grid." << endl;
        cout << "The expected eigenvalue is " << DimensionSparseGrid << "*pi^2 = "
             << setprecision(12) << eig_exp << setprecision(6) << endl;
        cout << "With epsilon = " << eps << endl;
        cout << "With cg_epsilon = " << cg_eps << endl;
        cout << "level start " << level_start << endl;
        cout << "level end " << level_end << endl;
        cout << endl;
        cout << "level\tDOFS\t\teig\t\t\terr\t\trate\terr_rel\t\tL2(u)\t\trate\tpow_it\tcg_it\ttime" << endl;
        cout << endl;
    }

    double err_old = 0.0;
    double l2_old = 0.0;

    for (int level = level_start; level < level_end + 1; level++) {

        // Regular grid construction
        AdaptiveSparseGrid grid;
        IndexDimension centerPoint;
        grid.AddRecursiveSonsOfPoint(centerPoint, level);

        MultiLevelAdaptiveSparseGrid mgrid(&grid);
        // symmetry on; distribute the stencil work with MPI only if we really run
        // with more than one task (the MPI path is not available in a serial build)
        MatrixVectorHomogen m(grid, mgrid, true, num_tasks > 1);

        // Matrix construction: A = Poisson (stiffness), M = HelmHoltz (mass)
        Poisson lhs(grid);
        HelmHoltz rhs(grid);

        // fields for solving
        double eigenvalue_power = 0.0;
        int pow_iterations = 0;
        double cg_iterations = 0.0;
        double time_power = 0.0;

        VectorSparseG prew_u(grid);
        prew_u = 1.0;

        // SOLVING
        power_inverse_laplace<Poisson, HelmHoltz>(prew_u, eigenvalue_power, pow_iterations,
                                                  cg_iterations, time_power, m, lhs, rhs);

        double err = abs(eigenvalue_power - eig_exp);
        double err_rel = err / eig_exp;

        // Error of the eigenfunction: scale the computed eigenvector such that it
        // attains the value of the exact eigenfunction in the center point (0.5,..,0.5),
        // which removes the sign/scaling ambiguity of the eigenvector.
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

        // convergence rate err(level-1)/err(level), undefined on the first level
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
                 << "\t" << setprecision(2) << time_power
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
