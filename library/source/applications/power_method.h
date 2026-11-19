//
// Created by rumi on 19.12.24
//

#ifndef GRUN_Power_METHOD_H
#define GRUN_Power_METHOD_H


#include <sys/time.h>
#include <unordered_set>

#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "../MatrixVectorMultiplication/MatrixVectorInhomogen.h"
#include "cg_method.h"


void saveLevelResults(int current_level, int iteration, double ev2, double val, double cg_time, int cg_iterations) {

    // 1. Job-ID für eindeutigen Dateinamen holen
    const char* raw_id = std::getenv("SLURM_JOB_ID");
    std::string job_id = (raw_id != nullptr) ? raw_id : "local";
    std::string filename = "results_job_" + job_id + ".txt";

    // 2. Zustand speichern, um Level-Wechsel zu erkennen
    static int last_level = -1; 
    bool level_changed = (last_level != -1 && last_level != current_level);
    last_level = current_level;

    // 3. Datei im Append-Modus öffnen
    std::ofstream file(filename, std::ios::app);

    if (file.is_open()) {
        // Falls das Level gewechselt hat: Leerzeile und Header einfügen
        if (level_changed) {
            file << "\n" << std::string(50, '-') << "\n";
            file << "Starting Level " << current_level << "\n\n";
        }

        // 4. Daten schreiben (Scientific Notation)
        file << std::scientific << std::setprecision(12)
             << "\tIt: " << std::setw(4) << iteration 
             << "\tEV: " << ev2 
             << "\tDiff: "<< val 
             << "\tCG: " << cg_iterations << std::endl;

        file.close();
    }
}

template<class Stencil_left, class Stencil_right>
double calculate_residual(VectorSparseG &x, MatrixVectorHomogen &matrix, Stencil_left lhs,Stencil_right rhs, double lambda) {
    
    // 1. Berechne Ax = Matrix * x
    VectorSparseG Ax(x.getSparseGrid());
    VectorSparseG Mx(x.getSparseGrid());
    matrix.multiplication<Stencil_left>(x, Ax, lhs);
    matrix.multiplication<Stencil_right>(x, Mx, rhs);
    
    double lambda2 = product(x, Ax) / product(x, Mx);
    
    // 3. Berechne r = Ax - lambda * x
    VectorSparseG r(x.getSparseGrid());
    // r = Ax - lambda * x
    // Da VectorSparseG vermutlich keine direkte Vektor-Subtraktion mit Skalar hat:
    r = Ax-lambda2*Mx; // Kopie

   
    // 4. Berechne ||r||_2 = sqrt(r^T * r)
    return sqrt(product(r, r));
}


class Power {
public:
   

    template<class Stencil_left, class Stencil_right>
    static bool power_max(double eps, VectorSparseG &x,
                            double &eigenvalue, int &iterations,
                            MatrixVectorHomogen& matrix,Stencil_left lhs, Stencil_right rhs,
                            bool rightoperator, double &precondition, double &duration);

    template<class Stencil_left, class Stencil_right>
    static bool power_inverse(double eps,double cg_eps, VectorSparseG &x,
                            double &eigenvalue, int &iterations,
                            double &avg_cg_iterations, MatrixVectorHomogen& matrix,
                            Stencil_left lhs, Stencil_right rhs,
                            bool rightoperator, double &precondition, double &duration);

    template<class Stencil_left, class Stencil_right>
    static  bool power_inverse_riccarda(double eps,double cg_eps, VectorSparseG &x,
                            double &eigenvalue, int &iterations,
                            double &avg_cg_iterations, MatrixVectorHomogen& matrix,
                            Stencil_left lhs, Stencil_right rhs,
                            bool rightoperator, double &precondition, double &duration);

    /// @brief shifted inverse power iteration for A x = lambda * M x.
    /// Every step solves (A - shift*M) y = M x_k by preconditioned CG and
    /// normalises x_{k+1} = y / ||y||. The iteration converges to the
    /// eigenpair whose eigenvalue lies closest to shift. Since CG needs a
    /// symmetric positive definite system, shift must lie BELOW the smallest
    /// eigenvalue of the problem.
    /// @param eps stopping criterion: change in eigenvalue and in x_k
    /// @param cg_eps stopping criterion of the inner CG solves
    /// @param shift the shift sigma, handed over by the caller
    /// @param x in: initial guess, out: normalised eigenvector
    /// @param eigenvalue out: Rayleigh quotient x^T A x / x^T M x
    /// @param iterations out: number of power iterations
    /// @param avg_cg_iterations out: average number of CG iterations per solve
    /// @param rightoperator if false, M is replaced by the identity
    /// @param duration out: wall-clock time in seconds
    /// @param max_iterations maximal number of power iterations
    /// @param verbose print one line per power iteration (rank 0 only)
    /// @param shift_preconditioner if true, CG is preconditioned with the
    /// diagonal of the shifted operator, diag(A) - shift*diag(M); if false,
    /// with diag(A) only. Falls back to diag(A) if the shifted diagonal is
    /// not positive (which can only happen for shift > smallest eigenvalue).
    template<class Stencil_left, class Stencil_right>
    static bool power_inverse_shifted(double eps, double cg_eps, double shift, VectorSparseG &x,
                            double &eigenvalue, int &iterations,
                            double &avg_cg_iterations, MatrixVectorHomogen& matrix,
                            Stencil_left lhs, Stencil_right rhs,
                            bool rightoperator, double &duration,
                            int max_iterations = 50, bool verbose = true,
                            bool shift_preconditioner = true);

    template<class Stencil_left, class Stencil_right>
    static double eigenvalue_exp(MatrixVectorHomogen &matrix, Stencil_left &lhs, Stencil_right &rhs, VectorSparseG &u, AdaptiveSparseGrid &grid);
   
    

    // template<class Stencil>
    // double precon(Stencil stencil, AdaptiveSparseGrid* grid);

    static const int maxIteration = 500;
    static double error[MaximumDepth][maxIteration];
    static int count;

template<class Stencil>
static double precon(VectorSparseG z, Stencil stencil);

    
};


template<class Stencil_left, class Stencil_right>
double 
Power::eigenvalue_exp(MatrixVectorHomogen &matrix, Stencil_left &lhs, Stencil_right &rhs, VectorSparseG &u, AdaptiveSparseGrid &grid) {

    VectorSparseG left(grid), right(grid);
    matrix.multiplication<Stencil_left>(u, left, lhs);
    matrix.multiplication<Stencil_right>(u, right, rhs);

    double scal_left = sqrt(product(left,left));
    double scal_right = sqrt(product(right,right));
    return scal_left / scal_right;
}



template<class Stencil>
double Power::precon(VectorSparseG z, Stencil stencil) {

    // Start measuring time of preconditioning
    struct timeval begin_pre, end_pre;
    gettimeofday(&begin_pre, 0);

    //dirichlet_grid.completeDirichletGrid();  ???
    Preconditioning<Stencil> P(z,stencil);

    gettimeofday(&end_pre, 0);
    long pre_seconds = end_pre.tv_sec - begin_pre.tv_sec;
    long pre_microseconds = end_pre.tv_usec - begin_pre.tv_usec;
    double pre_duration_def = pre_seconds + pre_microseconds*1e-6;
    
    return pre_duration_def;
}



//INLINE TEMPLATE FUNCTIONS
template<class Stencil_left, class Stencil_right>
bool
Power::power_max(double eps, VectorSparseG &x, double &eigenvalue, int &iterations,
                 MatrixVectorHomogen& matrix, Stencil_left lhs, Stencil_right rhs, bool rightoperator, double &precondition, double &duration) {


    AdaptiveSparseGrid* grid = (x.getSparseGrid());
    //ListOfDepthOrderedSubgrids list(*grid);

    VectorSparseG z(grid);
    VectorSparseG xalt(grid);
    VectorSparseG y(grid);
    VectorSparseG xneu(grid);
    xalt = x / sqrt(product(x,x));

    // solves  Ax = M * (lambda * x); for lambda and x

       // rho = Raiyleigh-Quotient
    double rho=0.0;
    double rho_old=1.0;



    // Preconditioning
    if (precondition != 0) {
        precondition = Power::precon<Stencil_right>(z, rhs);
    }


    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);


    // norm of vector y
    double norm_x = 0;
    double norm_x_old = 0;


    int k = 0;

    // parameters for cg iteration to solve rhs
    int cg_iterations = 0;
    double cg_time = 0.0;
    double cg_eps = eps * 1e-2;
    
    // cout << endl;

    // Iterate until maxiter
    for (int i = 0; i < maxIteration; i++) {
        // Remember iterations
        k = i;

        // y := A * xalt
        matrix.multiplication<Stencil_left>(xalt, y, lhs);
        if (rightoperator) {
            // solve (y = M * x_neu) for xneu
            CG::solveHomogen<Stencil_right>(cg_eps, xneu, y, &cg_iterations, xalt, matrix, rhs, &cg_time);
        } else {
            xneu = y;
        }

        norm_x = sqrt(product(xneu, xneu)); /// grid->getDOFS());

        if (abs(norm_x - norm_x_old) < eps)
            break;
        

        // Update norm_old and xalt
        norm_x_old = norm_x;

        // normalize vector x := xneu / |xneu|
        xalt = xneu / norm_x;


    }


    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    duration = seconds + microseconds*1e-6;


    eigenvalue = norm_x;
    iterations = k;
    x = xneu/norm_x;

    return true;
}



//INLINE TEMPLATE FUNCTIONS
template<class Stencil_left, class Stencil_right>
bool
Power::power_inverse(double eps,double cg_eps, VectorSparseG &x, double &eigenvalue, int &iterations, double &avg_cg_iterations,
                 MatrixVectorHomogen& matrix, Stencil_left lhs, Stencil_right rhs, bool rightoperator, double &precondition, double &duration) {

    double border = 7.5;

    double rhs_factor = 8.0*border*border;
    const double shiftfactor = 0.9;


    AdaptiveSparseGrid* grid = (x.getSparseGrid());


    VectorSparseG z(grid);
    VectorSparseG diff(grid);
    VectorSparseG xalt(grid);
    VectorSparseG y(grid);
    VectorSparseG xneu(grid);
    VectorSparseG xneu_norm(grid);
    
    xalt = x;


    // rho = Raiyleigh-Quotient
    double rho=0.0;
    double rho_old=1.0;

    // solves  Ax = M * (lambda * x); for smallest lambda and x

    // Preconditioning



    Preconditioning<Stencil_left>  preconditioning(z,lhs);
    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);


    // seminorm of vector xneu
    double norm_x = 0;
    double norm_x_old = 0;

    int k = 0;

    // (return) parameters for cg iteration to solve rhs
    int cg_iterations = 0;
    double cg_time = 0.0;
    //double cg_eps = 1e-10;
    double avg_cg_iter = 0;
                int rank;
   

    

    // Iterate until maxiter
    for (int i = 0; i < 50; i++) {
        // Remember iterations
        k = i;

        if (rightoperator) {
            // y := M * xalt
            matrix.multiplication<Stencil_right>(xalt, y, rhs);
        } else {
            y = xalt;
        }
        // solve (y = A * x_neu) for xneu
        CG::solveHomogenP<Stencil_left>(cg_eps, xneu, y, &cg_iterations, xalt, matrix, lhs, &cg_time,preconditioning);
        if(cg_iterations < 3 && cg_eps > 1e-10) cg_eps=0.1*cg_eps;


        norm_x = (product(xneu, xneu));
              rho = product(xalt,xneu)/norm_x;
        norm_x=sqrt(norm_x);
        
        xneu_norm = xneu/norm_x;

        rho = product(xalt,xneu);

        diff = xneu/norm_x - xalt;
        double norm_test = sqrt(product(diff,diff));
        
        if(abs(rho_old-rho)<eps && norm_test<eps) break;
               double val =abs(rho_old-rho);
        rho_old = rho;
 
        //if (abs(norm_x - norm_x_old) < eps)





        
        norm_x_old = norm_x;
        xalt = xneu / norm_x;

        avg_cg_iter += cg_iterations;

                int rank = 0;
        #ifdef MY_MPI_ON
    
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
#endif
         eigenvalue = 1.0 / norm_x;
         double res = calculate_residual<Stencil_left,Stencil_right>(xalt, matrix, lhs,rhs,eigenvalue);
        
    
        double ev=eigenvalue/(rhs_factor*shiftfactor);
        double ev2 = 1.0/rho;
        ev2=ev2/(rhs_factor*shiftfactor);
     
    if (rank==0)cout << "\tPower Iteration " << to_string(i)<< "\t abs(rho_old-rho) " << std::scientific <<val<< "\t res " << std::scientific <<res<< "\t diff " << to_string(norm_test) << "\t cg_time"  << to_string(cg_time) << "\t cg_iterations " <<to_string(cg_iterations) << "\t cg eps" << cg_eps <<"\t eigenvalue" << ev<<"\t rho" << rho << "\t eigenvalue rho " << to_string(ev2)<< endl;
    }

    avg_cg_iter /= k;


    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    duration = seconds + microseconds*1e-6;


    eigenvalue = 1.0 / norm_x;
    iterations = k;
    avg_cg_iterations = avg_cg_iter;
    x = xneu/norm_x;

    return true;
}



//INLINE TEMPLATE FUNCTIONS
template<class Stencil_left, class Stencil_right>
bool
Power::power_inverse_riccarda(double eps,double cg_eps, VectorSparseG &x, double &eigenvalue, int &iterations, double &avg_cg_iterations,
                 MatrixVectorHomogen& matrix, Stencil_left lhs, Stencil_right rhs, bool rightoperator, double &precondition, double &duration) {


        // 1. Dateinamen generieren (z.B. ergebnisse_level_0.txt)
    
    AdaptiveSparseGrid* grid = (x.getSparseGrid());
    int level = grid->getLevel();
    //std::string fileName = "ergebnisse_level_neu" + std::to_string(level) + ".txt";

    // 2. Datei im "Append"-Modus öffnen, damit neue Iterationen 
    // desselben Levels in die gleiche Datei geschrieben werden.
    //std::ofstream outFile(fileName, std::ios::app);

 



    //outFile << "DOFS: " << grid->getDOFS() << std::endl;

    VectorSparseG z_right(grid),z(grid);

    VectorSparseG diff(grid);
    VectorSparseG xalt(grid);
    VectorSparseG y(grid);
    VectorSparseG xneu(grid);
    VectorSparseG xneu_norm(grid);
    
    xalt = x;




    // solves  Ax = M * (lambda * x); for smallest lambda and x




    int rank;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    // Preconditioning
    struct timeval begin, end;
    gettimeofday(&begin, 0);
    Preconditioning<Stencil_left>  preconditioning_left(z,lhs);
    Preconditioning<Stencil_right>  preconditioning_right(z_right,rhs,true);
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double duration_def = seconds + microseconds*1e-6;
    precondition = duration_def;
  //if (rank == 0 ) cout << "time precond " << duration_def << endl;
    gettimeofday(&begin, 0);



    int k = 0;

    // (return) parameters for cg iteration to solve rhs
    int cg_iterations = 0;
    double cg_time = 0.0;





    double lambda_old=0.0;
    double lambda=0.0;


    double sigma =200.0 ;
    double finalsigma =sigma;

    // Iterate until maxiter
    for (int i = 0; i < maxIteration; i++)
    {
   
        // Remember iterations
        k = i;
    
        // y := M * xalt
        matrix.multiplication<Stencil_right>(xalt, y, rhs);


        // solve (A-sigma*B)xneu=y
        CG::solveHomogenShift2<Stencil_left,Stencil_right>(cg_eps, xneu, y, &cg_iterations, xalt, matrix, lhs,rhs, &cg_time,sigma, preconditioning_left, preconditioning_right);

        VectorSparseG Mx(grid);
        matrix.multiplication<Stencil_right>(xneu, Mx, rhs); // z_right = M * xneu


        VectorSparseG Ax(grid);
        matrix.multiplication<Stencil_left>(xneu,Ax, lhs); // z_right = M * xneu



        double m_norm = sqrt(product(xneu, xneu));
        xalt = xneu / m_norm;

        lambda = product(xneu,Ax)/product(xneu,Mx);






        double val = abs(lambda- lambda_old);

        lambda_old = lambda;
    
        if(rank==0) saveLevelResults(level,i,lambda,0.0,cg_time,cg_iterations);


        if(abs(val)<eps) break;

        iterations = k;
    }

    eigenvalue = lambda;
    x = xneu/sqrt(product(xneu,xneu));

    gettimeofday(&end, 0);
    seconds = end.tv_sec - begin.tv_sec;
    microseconds = end.tv_usec - begin.tv_usec;
    duration_def = seconds + microseconds*1e-6;
    duration = duration_def;
    //if (rank == 0 ) cout << "time power " << duration_def << endl;


    return true;

}


/// @brief preconditioned CG for the shifted system (A - shift*M) x = f.
/// A and M are applied separately, so no combined stencil is needed.
/// D_inv holds the inverse of the diagonal (Jacobi) preconditioner.
/// @return number of CG iterations
template<class Stencil_left, class Stencil_right>
int solve_shifted_system(double cg_eps, double shift, VectorSparseG &x, VectorSparseG &f,
                         MatrixVectorHomogen &matrix, Stencil_left &lhs, Stencil_right &rhs,
                         bool rightoperator, VectorSparseG &D_inv) {
    AdaptiveSparseGrid *grid = x.getSparseGrid();
    VectorSparseG r(grid), h(grid), d(grid), g(grid), Md(grid);

    // g := (A - shift*M) * v; for shift == 0 the M product is skipped
    auto apply = [&](VectorSparseG &v, VectorSparseG &out) {
        out = 0.0;
        matrix.multiplication<Stencil_left>(v, out, lhs);
        if (shift == 0.0) return;
        if (rightoperator) {
            Md = 0.0;
            matrix.multiplication<Stencil_right>(v, Md, rhs);
        } else {
            Md = v;
        }
        out = out - shift * Md;
    };

    // r = A_shift*x - f
    apply(x, r);
    r = r - f;
    h = r * D_inv;
    d = -1.0 * h;
    double delta = product(r, h);

    const int max_cg_iterations = 1000;
    int k = 0;
    for (int i = 1; i <= max_cg_iterations && delta > cg_eps; ++i) {
        grid->WorkOnHangingNodes = false;
        k = i;

        apply(d, g);
        double dg = product(d, g);
        if (dg <= 0.0) {
            int rank = 0;
#ifdef MY_MPI_ON
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
            if (rank == 0)
                cout << "WARNING: A - shift*M is not positive definite (shift too large?), CG stopped." << endl;
            break;
        }
        double tau = delta / dg;

        r = r + (tau * g);
        x = x + (tau * d);
        h = r * D_inv;

        double delta_prime = product(r, h);
        if (delta_prime < cg_eps) break;

        double beta = delta_prime / delta;
        delta = delta_prime;
        d = -1.0 * h + beta * d;
    }
    return k;
}


template<class Stencil_left, class Stencil_right>
bool
Power::power_inverse_shifted(double eps, double cg_eps, double shift, VectorSparseG &x,
                             double &eigenvalue, int &iterations, double &avg_cg_iterations,
                             MatrixVectorHomogen &matrix, Stencil_left lhs, Stencil_right rhs,
                             bool rightoperator, double &duration, int max_iterations, bool verbose,
                             bool shift_preconditioner) {

    int rank = 0;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    AdaptiveSparseGrid *grid = x.getSparseGrid();
    VectorSparseG z(grid), zM(grid), D_inv(grid);
    VectorSparseG xalt(grid), xneu(grid), y(grid), diff(grid), Ax(grid), Mx(grid);

    struct timeval begin, end;
    gettimeofday(&begin, 0);

    // Jacobi preconditioner: z := diag(A), and for the shift-aware variant
    // z := diag(A) - shift*diag(M). Every entry A_ii/M_ii is a Rayleigh
    // quotient and hence >= the smallest eigenvalue, so the shifted diagonal
    // is positive whenever shift lies below the smallest eigenvalue.
    Preconditioning<Stencil_left> preconditioning(z, lhs);
    if (shift_preconditioner && shift != 0.0) {
        if (rightoperator) {
            Preconditioning<Stencil_right> preconditioning_right(zM, rhs); // fills zM = diag(M)
        } else {
            zM = 1.0;
        }
        bool positive = true;
        for (unsigned long k = 0; k < grid->getMaximalOccupiedSecondTable(); k++) {
            if (grid->getActiveTable()[k] && z.getValue(k) - shift * zM.getValue(k) <= 0.0) {
                positive = false;
                break;
            }
        }
        if (positive) {
            z = z - shift * zM;
        } else if (rank == 0) {
            cout << "WARNING: diag(A) - shift*diag(M) is not positive (shift too large?), "
                 << "using diag(A) as preconditioner." << endl;
        }
    }
    D_inv = 0.0;
    for (unsigned long k = 0; k < grid->getMaximalOccupiedSecondTable(); k++) {
        if (grid->getActiveTable()[k]) D_inv.setValue(k, 1.0 / z.getValue(k));
    }

    xalt = x / sqrt(product(x, x));

    // Rayleigh quotient of the initial guess, used for the first warm start
    Ax = 0.0;
    matrix.multiplication<Stencil_left>(xalt, Ax, lhs);
    if (rightoperator) {
        Mx = 0.0;
        matrix.multiplication<Stencil_right>(xalt, Mx, rhs);
    } else {
        Mx = xalt;
    }
    double lambda = product(xalt, Ax) / product(xalt, Mx);
    double lambda_old = lambda;
    double sum_cg_iterations = 0.0;
    int k = 0;

    for (int i = 0; i < max_iterations; i++) {
        k = i + 1;

        // y := M * x_k
        if (rightoperator) {
            y = 0.0;
            matrix.multiplication<Stencil_right>(xalt, y, rhs);
        } else {
            y = xalt;
        }

        // solve (A - shift*M) x_{k+1} = y. If x_k is an eigenvector to
        // lambda, the solution is x_k / (lambda - shift), so this is the warm
        // start. Starting from x_k itself would be off by that factor
        // (e.g. ~30 for shift 0 here) and leave CG with a large error.
        double lambda_shifted = lambda_old - shift;
        if (lambda_shifted > 0.0) xneu = xalt / lambda_shifted;
        else xneu = xalt;
        int cg_iterations = solve_shifted_system<Stencil_left, Stencil_right>(
                cg_eps, shift, xneu, y, matrix, lhs, rhs, rightoperator, D_inv);
        sum_cg_iterations += cg_iterations;

        // normalise and keep the sign of the previous iterate
        xneu = xneu / sqrt(product(xneu, xneu));
        if (product(xneu, xalt) < 0.0) xneu = -1.0 * xneu;

        // Rayleigh quotient and residual ||A x - lambda M x||
        Ax = 0.0;
        matrix.multiplication<Stencil_left>(xneu, Ax, lhs);
        if (rightoperator) {
            Mx = 0.0;
            matrix.multiplication<Stencil_right>(xneu, Mx, rhs);
        } else {
            Mx = xneu;
        }
        lambda = product(xneu, Ax) / product(xneu, Mx);
        diff = Ax - lambda * Mx;
        double res = sqrt(product(diff, diff));

        diff = xneu - xalt;
        double norm_diff = sqrt(product(diff, diff));
        double change = abs(lambda - lambda_old);

        if (verbose && rank == 0)
            cout << "\tShifted Power Iteration " << i
                 << "\t shift " << shift
                 << "\t eigenvalue " << std::scientific << std::setprecision(12) << lambda
                 << "\t abs(lambda_old-lambda) " << std::setprecision(3) << change
                 << "\t diff " << norm_diff
                 << "\t res " << res
                 << "\t cg_iterations " << cg_iterations
                 << std::defaultfloat << std::setprecision(6) << endl;

        xalt = xneu;
        lambda_old = lambda;

        if (i > 0 && change < eps && norm_diff < eps) break;
    }

    gettimeofday(&end, 0);
    duration = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec) * 1e-6;

    eigenvalue = lambda;
    iterations = k;
    avg_cg_iterations = sum_cg_iterations / k;
    x = xalt;

    return true;
}


#endif //GRUN_CG_METHOD_H
