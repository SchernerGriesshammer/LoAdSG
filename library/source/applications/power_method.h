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

                int rank;
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




#endif //GRUN_CG_METHOD_H
