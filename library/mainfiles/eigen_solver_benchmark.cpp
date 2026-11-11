#include "../source/sparseEXPDE.h"
#include "../source/iterator/depthIterator.h"
#include "../source/applications/power_method.h"
#include "../source/applications/lobpcg.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sys/time.h>

namespace {

constexpr double border = 7.5;
constexpr double cut_coeff = 1e4;
constexpr double target_energy = -2.0;
constexpr double search_point = 200.0;
constexpr double shift_margin = 2.5;

constexpr double varc_factor = -8.0 * border;
constexpr double rhs_factor = 8.0 * border * border;
constexpr double shiftvalue =
    search_point + shift_margin - target_energy * rhs_factor;

constexpr int grid_level = 4;
constexpr int matvec_process_count = 2;
constexpr int minimum_process_count = matvec_process_count + 1;
constexpr int potential_mc_sample_factor = 100;
constexpr int random_seed = 649454244;
constexpr double eigenvalue_tolerance = 1e-5;
constexpr double cg_tolerance = 1e-10;
constexpr bool use_symmetry = true;
constexpr bool use_mpi_distribution = true;
constexpr bool use_mass_operator_on_right = true;

using AssembledOperator = LocalStiffnessMatricesDynamicDistribution;
using PotentialFunction = double (*)(double*);
using PotentialStencil = StencilMC<PotentialFunction>;
using MassStencil = HelmHoltz;

double potential_at(double* coordinates) {
    double radius_squared = 0.0;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        const double c = 2.0 * coordinates[d] - 1.0;
        radius_squared += c * c;
    }

    const double radius = std::sqrt(radius_squared);
    const double inverse_radius =
        (radius > 1e-15) ? std::min(1.0 / radius, cut_coeff) : cut_coeff;

    return varc_factor * inverse_radius + shiftvalue;
}

double elapsed_seconds(const timeval& begin, const timeval& end) {
    return static_cast<double>(end.tv_sec - begin.tv_sec)
         + static_cast<double>(end.tv_usec - begin.tv_usec) * 1e-6;
}

struct EigenResult {
    const char* solver_name = "";
    double eigenvalue = 0.0;
    double physical_energy = 0.0;
    double residual = 0.0;
    int iterations = 0;
    double preconditioning_seconds = 0.0;
    double solve_seconds = 0.0;
    double wall_seconds = 0.0;
    bool converged = false;
};

struct ComparisonResult {
    int degrees_of_freedom = 0;
    EigenResult inverse_power;
    EigenResult lobpcg;
};

class ProblemOperators {
public:
    explicit ProblemOperators(AdaptiveSparseGrid& grid,
                              int assembly_process_count)
        : laplacian_stencil_(grid),
          mass_stencil_(grid),
          potential_stencil_(grid, &potential_at,
                             potential_mc_sample_factor, random_seed),
          hamiltonian_(grid, potential_stencil_, assembly_process_count),
          mass_operator_(grid, mass_stencil_, assembly_process_count) {
        hamiltonian_.addStencil(laplacian_stencil_);
    }

    ProblemOperators(const ProblemOperators&) = delete;
    ProblemOperators& operator=(const ProblemOperators&) = delete;
    ProblemOperators(ProblemOperators&&) = delete;
    ProblemOperators& operator=(ProblemOperators&&) = delete;

    AssembledOperator& hamiltonian() { return hamiltonian_; }
    const AssembledOperator& hamiltonian() const { return hamiltonian_; }

    AssembledOperator& mass_operator() { return mass_operator_; }
    const AssembledOperator& mass_operator() const { return mass_operator_; }

private:
    Poisson laplacian_stencil_;
    MassStencil mass_stencil_;
    PotentialStencil potential_stencil_;
    AssembledOperator hamiltonian_;
    AssembledOperator mass_operator_;
};

void initialize_vector(VectorSparseG& vector) {
    vector = 1.0;
}

void finalize_result(EigenResult& result, VectorSparseG& eigenvector,
                     AdaptiveSparseGrid& grid,
                     MatrixVectorHomogen& matrix_vector,
                     ProblemOperators& operators) {
    result.physical_energy = (result.eigenvalue - shiftvalue) / rhs_factor;

    VectorSparseG left_product(grid);
    VectorSparseG right_product(grid);
    VectorSparseG residual_vector(grid);

    matrix_vector.multiplication(eigenvector, left_product,
                                 operators.hamiltonian());
    matrix_vector.multiplication(eigenvector, right_product,
                                 operators.mass_operator());

    residual_vector = left_product - result.eigenvalue * right_product;
    result.residual = std::sqrt(product(residual_vector, residual_vector));
}

EigenResult run_inverse_power(AdaptiveSparseGrid& grid,
                              MatrixVectorHomogen& matrix_vector,
                              ProblemOperators& operators) {
    EigenResult result;
    result.solver_name = "inverse power";

    VectorSparseG eigenvector(grid);
    initialize_vector(eigenvector);

    timeval begin;
    gettimeofday(&begin, 0);

    double unused_cg_iterations = 0.0;

    result.converged = Power::power_inverse_riccarda<
        AssembledOperator,
        AssembledOperator>(
            eigenvalue_tolerance, cg_tolerance, eigenvector,
            result.eigenvalue, result.iterations,
            unused_cg_iterations, matrix_vector,
            operators.hamiltonian(), operators.mass_operator(),
            use_mass_operator_on_right, result.preconditioning_seconds,
            result.solve_seconds);

    timeval end;
    gettimeofday(&end, 0);
    result.wall_seconds = elapsed_seconds(begin, end);

    finalize_result(result, eigenvector, grid, matrix_vector, operators);

    return result;
}

EigenResult run_lobpcg(AdaptiveSparseGrid& grid,
                       MatrixVectorHomogen& matrix_vector,
                       ProblemOperators& operators) {
    EigenResult result;
    result.solver_name = "lobpcg";

    VectorSparseG eigenvector(grid);
    initialize_vector(eigenvector);

    timeval begin;
    gettimeofday(&begin, 0);

    double unused_cg_iterations = 0.0;

    result.converged = Lobpcg::solve<
        AssembledOperator,
        AssembledOperator>(
            eigenvalue_tolerance, cg_tolerance, eigenvector,
            result.eigenvalue, result.iterations,
            unused_cg_iterations, matrix_vector,
            operators.hamiltonian(), operators.mass_operator(),
            use_mass_operator_on_right, result.preconditioning_seconds,
            result.solve_seconds);

    timeval end;
    gettimeofday(&end, 0);
    result.wall_seconds = elapsed_seconds(begin, end);

    finalize_result(result, eigenvector, grid, matrix_vector, operators);

    return result;
}

[[nodiscard]] ComparisonResult solve_eigenvalue_problem(
    int assembly_process_count) {
    AdaptiveSparseGrid grid;
    IndexDimension center_point;
    grid.AddRecursiveSonsOfPoint(center_point, grid_level);

    MultiLevelAdaptiveSparseGrid multilevel_grid(&grid);
    MatrixVectorHomogen matrix_vector(
        grid, multilevel_grid, use_symmetry, use_mpi_distribution);

    ProblemOperators operators(grid, assembly_process_count);

    ComparisonResult comparison;
    comparison.degrees_of_freedom = grid.getDOFS();
    comparison.inverse_power =
        run_inverse_power(grid, matrix_vector, operators);
    comparison.lobpcg = run_lobpcg(grid, matrix_vector, operators);

    return comparison;
}

void print_result(const EigenResult& result) {
    std::cout << result.solver_name << ":\n";
    std::cout << "  converged: "
              << (result.converged ? "yes" : "no") << '\n';
    std::cout << "  physical energy: " << result.physical_energy << '\n';
    std::cout << "  absolute error: "
              << std::abs(result.physical_energy - target_energy) << '\n';
    std::cout << "  residual: " << result.residual << '\n';
    std::cout << "  iterations: " << result.iterations << '\n';
    std::cout << "  preconditioning_seconds: "
              << result.preconditioning_seconds << '\n';
    std::cout << "  solve_seconds: " << result.solve_seconds << '\n';
    std::cout << "  wall_seconds: " << result.wall_seconds << '\n';
}

void print_comparison(const ComparisonResult& comparison) {
    std::cout << "DOFs: " << comparison.degrees_of_freedom << '\n';
    std::cout << "target energy: " << target_energy << "\n\n";

    print_result(comparison.inverse_power);
    std::cout << '\n';
    print_result(comparison.lobpcg);

    std::cout << "\ncomparison:\n";
    std::cout << "  physical energy difference: "
              << std::abs(comparison.lobpcg.physical_energy
                          - comparison.inverse_power.physical_energy) << '\n';

    const double ratio =
        (comparison.lobpcg.solve_seconds > 0.0)
            ? comparison.inverse_power.solve_seconds
                  / comparison.lobpcg.solve_seconds
            : 0.0;
    std::cout << "  solve time ratio (power/lobpcg): " << ratio << '\n';
}

#ifdef MY_MPI_ON
class MpiSession {
public:
    MpiSession(int& argc, char**& argv) {
        MPI_Init(&argc, &argv);
    }

    ~MpiSession() {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) {
            MPI_Finalize();
        }
    }

    MpiSession(const MpiSession&) = delete;
    MpiSession& operator=(const MpiSession&) = delete;
};
#else
class MpiSession {
public:
    MpiSession(int&, char**&) {}
};
#endif

int process_rank() {
#ifdef MY_MPI_ON
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
#else
    return 0;
#endif
}

int process_count() {
#ifdef MY_MPI_ON
    int count = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    return count;
#else
    return 1;
#endif
}

[[noreturn]] void report_fatal(const char* message) {
    if (process_rank() == 0) {
        std::cerr << "Fatal error: " << message << '\n';
    }

#ifdef MY_MPI_ON
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif

    std::exit(EXIT_FAILURE);
}

int run_benchmark() {
    const int count = process_count();

    if (count < minimum_process_count) {
        if (process_rank() == 0) {
            std::cerr << "Run this with at least " << minimum_process_count
                      << " MPI processes.\n";
        }
        return EXIT_FAILURE;
    }

    const int assembly_process_count = count - matvec_process_count;
    const ComparisonResult comparison =
        solve_eigenvalue_problem(assembly_process_count);

    if (process_rank() == 0) {
        print_comparison(comparison);
    }

    return (comparison.inverse_power.converged &&
            comparison.lobpcg.converged) ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace

int main(int argc, char** argv) {
    MpiSession mpi(argc, argv);
    (void)mpi;

    try {
        return run_benchmark();
    } catch (const std::exception& error) {
        report_fatal(error.what());
    } catch (...) {
        report_fatal("unknown exception");
    }
}
