//
// Created by to35jepo on 11/16/22.
//

#include "../source/sparseEXPDE.h"
#include "../source/applications/H1_Norm.h"
#include "../source/functions/functions.h"


//#define testNOW
//#define goolePerf
#include <string>
#include <ctime>
#include <sys/time.h>
#include <iomanip>

void mpi_cout(string s){
    int rank = 0;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank ==0){
        cout << s << endl;
    }
#else
    cout << s  << endl;
#endif


}

int count = 0;
double k = 100.0;
const double CUTCOEFF = 1e6;

double var_coeff( double* x)
{
  
    double r =0.0;
    
    for (int i = 0; i < DimensionSparseGrid; i++) {
        double ci = x[i] - 0.5; // < map from [0,1] to [-s,s]
        r += ci * ci;
    }

    r = sqrt(r);

    

        
    double result = CUTCOEFF;
    
    if (r > 1e-15){
        result = min(1. / r, CUTCOEFF); // < clipping large values of 1/r around the center
    }


    return -1.0*k*result;

}




int main(int argc, char **argv) {
    
    unsigned int seed = std::random_device{}(); // Fallback auf Zufall

    

    int rank = 0;
    int num_tasks = 1;
#ifdef MY_MPI_ON
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif



    mpi_cout(" Dimension " + to_string( DimensionSparseGrid));
    mpi_cout("regular grid, solve: -lap u + k u = f "); 
    mpi_cout("varCoeff, a=100, c=1e+8, k =1");


    int level_start=1;



    //MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    #pragma omp parallel
    {
        #pragma omp single
        {
            int total_threads = omp_get_num_threads();
            mpi_cout("number of mpi_tasks " + to_string(num_tasks));
            mpi_cout("number of omp_threads " + to_string(total_threads));
        }
    }
    
    
    int numberMVprocesses=10;
    int numberLSprocesses=num_tasks - numberMVprocesses;
    mpi_cout("use " + to_string(numberLSprocesses)+ " for LocalStiffnessmatrices ");
    mpi_cout("use " + to_string(num_tasks - numberLSprocesses)+ " for MV ");


    double Linfty_old = 1.0;
    double L2_old = 1.0;
    double H1_old =1.0;



    IndexDimension centerPoint;


    AdaptiveSparseGrid old_grid;
    old_grid.AddRecursiveSonsOfPoint(centerPoint,level_start);

    //mpi_cout("var coeff = 1/(r+1e-5) ");
    //mpi_cout("number of mc points on each cell multiplied with 100");
    for (int level = level_start; level <10; level++)
    {


          AdaptiveSparseGrid grid;
        grid.AddRecursiveSonsOfPoint(centerPoint,level);
        //grid.copy_inner(old_grid);
        





        MultiLevelAdaptiveSparseGrid mgrid(&grid);



        VectorSparseG u(grid), f_vec(grid), f_test(grid),fpu(grid);
        u=0.0;


        double val_u = 0.0;
        double val_f = 0.0;
        double val_rhs = 0.0;

        double a = -100.0;
        double b = 0.5;
        double c=1e8;

      for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++) {
            if(grid.getActiveTable()[i]) {
                IndexDimension I = grid.getIndexOfTable(i);
                double x[DimensionSparseGrid];
                for (int d = 0; d < DimensionSparseGrid; d++) x[d] = I.coordinate(d);

                double f2 = c*g(x)*f(x,b,a)*v(x);
       

                val_u = f2;
              

                double r2 = -1.0*c*laplacian_fgv(x,b,a);

                val_f = r2;
    

                u.setValue(i, val_u);
                f_vec.setValue(i, val_f);
                f_test.setValue(i, val_f);
                


            }

        }



        VectorSparseG cu(grid);

        MatrixVectorHomogen m(grid,mgrid,1,1);

        HelmHoltz helmHoltz(grid);
        Poisson stencilPoisson(grid);
        VectorSparseG prew(grid);
        prew = 0.0;
        

        
        calcPrewByNodal(prew,f_vec);


        

        //VariableCoefficientMass stencilVarCoeff(grid);

        //HelmHoltz stencilVarCoeff(grid);
        //ConstantCoefficientMass stencilVarCoeff(grid);
        //StencilSGpp<double (*)(double *)> stencilVarCoeff(grid, &var_coeff, level+1,MonteCarloFixed);

        StencilMC<double (*)(double *)> stencilVarCoeff(grid,&var_coeff,100,seed,10);









        LocalStiffnessMatricesDynamicDistribution lhs(grid, stencilVarCoeff,numberLSprocesses);
        LocalStiffnessMatricesDynamicDistribution rhs(grid,helmHoltz,numberLSprocesses);



        cu = 0.0;
        prew=0.0;



        calcPrewByNodal(prew,u);

        m.multiplication(prew,cu,lhs);




        prew = 0.0;


        calcPrewByNodal(prew,f_test);
        //calcPrewByNodal(prew,rhs_vec);
        f_test=0.0;
        m.multiplication(prew,f_test,rhs);



        f_vec = f_test+cu;
        //    f_vec=f_test;


        int iterations;
        prew = 0.0;

        double time_p;

        lhs.addStencil(stencilPoisson);
               
        //ADD<Poisson,HelmHoltz> stencil(stencilPoisson,stencilVarCoeff,grid);


        CG::solveHomogen(1e-10, prew, f_vec, &iterations, u, m, lhs, &time_p);


        time_p=0.0;


        grid.WorkOnHangingNodes=false;

        VectorSparseG w(grid),uSolution(grid),d(grid);
        w=0.0;
        calcNodalByPrew(prew, w);
        uSolution = w;


        d = u - uSolution;
        double linfty= L_infty(d);


        L2 L2_(grid,mgrid);
        double l2=L2_.getValue(d);

        H1 h1(grid,mgrid);
        double h1_ = h1.getValue(d);
        h1_ =l2+h1_;


        old_grid.clear();
        f_test=u;


       int rank=0;
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          if(rank ==0){
              cout << level << " " << grid.getDOFS() << "  " << linfty << " " << Linfty_old/linfty << "  " <<  l2 << " " << L2_old/l2 <<   "    " <<  h1_ << " " <<H1_old/h1_<< "  " << to_string(iterations) <<  endl;
          }

  

/*
        double condition1 = condition(100,grid,mgrid,true);
        cout <<level << "  " << condition1 << endl;
*/

        prew=0.0;
        u=0.0;

        if(rank==0){

            const char* raw_id = std::getenv("SLURM_JOB_ID");
            std::string job_id = (raw_id != nullptr) ? raw_id : "local";
            std::string filename = "temp_run_" + job_id + ".dat";
    
            // Append-Modus: Schreibt Level 1 bis 10 nacheinander in dieses File
            std::ofstream file(filename, std::ios::app);
       
    
        if (file.is_open()) {
                 if (level == 1) {
                        file << "\n"; // Die Leerzeile vor einem neuen Block
                    }

                file << job_id << level << " " << grid.getDOFS() << "  " << linfty << " " << to_string(Linfty_old/linfty) << "  " <<  l2 << " " << to_string(L2_old/l2) <<   "    " <<  h1_ << " " <<to_string(H1_old/h1_)<< "  " << iterations << "\n";
            file.close();
            }
        }

              Linfty_old = linfty;
        L2_old = l2;
        H1_old=h1_;
    }
#ifdef MY_MPI_ON
    MPI_Finalize();
#endif
    return 0;


}