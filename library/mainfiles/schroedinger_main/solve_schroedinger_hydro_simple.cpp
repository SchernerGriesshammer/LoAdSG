//
// Created by rumi on 23/09/25.
//

#include "../source/sparseEXPDE.h"
#include "../source/iterator/depthIterator.h"
#include "../source/applications/power_method.h"
//#define testNOW
//#define goolePerf
#include <string>
#include <ctime>
#include <sys/time.h>



double border = 7.5;
//3e-9;
//8e-10;

double eps = 1e-5;
double cg_eps = 1e-10;


const double CUTCOEFF = 1e4;

const int mc_samples = 100;

    // https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
    double bohr_radius = 5.29177210544e-11;
    double b_factor = sqrt(M_PI);


    const double varc_factor = -8.0*border;

    const double rhs_factor = 8.0*border*border;

    const double eig_exp = -0.5;//*rhs_factor;
    const double shiftfactor = 0.9;
    const double shiftvalue = (shiftfactor+1) * (-1) * eig_exp*rhs_factor;



volatile double r_min = 1;
volatile double r_max = 1;




void mpi_cout(string s, bool endline = true){
    int rank = 0;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank ==0){
        cout << s;
        if (endline) cout << endl;
        cout << flush;
    }
#else
        cout << s;
        if (endline) cout << endl;
        cout << flush;
#endif
}




string r_minmax_res_log() {
    string str = "";
    
    int rank = 0;
    #ifdef MY_MPI_ON
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        double globr_min;
        double globr_max;

        MPI_Reduce((const double *)&r_min, &globr_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce((const double *)&r_max, &globr_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank==0) {
            str = to_string(globr_min) + "\t" + to_string(globr_max);
        }
    #endif

    r_min = 1;
    r_max = 1;
    return str;
}

double alpha = 1.0;




double ExponentialCutOff(double r)
{
    if (r < 1e-15){
        return CUTCOEFF;
    } else {
        return (1.0 - std::exp(-CUTCOEFF * r)) / r;
    }

}


double var_coeff( double* coordinates)
{
    double r = 0;
    for (int i = 0; i < DimensionSparseGrid; i++) {
        double ci = 2. * coordinates[i] - 1.; // < map from [0,1] to [-s,s]
        r += ci * ci;
    }

    r = sqrt(r);

        
    double result = CUTCOEFF;
    //double rc = 1.0/CUTCOEFF;
    //if (r < rc) result = ExponentialCutOff(r);
    //    else result=1.0/r;

    if (r > 1e-15){
        result = min(1. / r, CUTCOEFF); // < clipping large values of 1/r around the center
    }else{
        result = CUTCOEFF;
    }

    result = varc_factor * result;
    //result -= 1.9*eig_exp; // + 1.9*1*4*s^2 (eig_exp==1)
    result += shiftvalue;

    return result;

    // log range of sampled values
    //double over_r = (r!=0) ? 1./r : 1;
    //if (alpha==1 && over_r < r_min) r_min = over_r;
    //if (alpha==1 && over_r > r_max) r_max = over_r;



    //return alpha*result;
}




inline double unitSQ2atom(double coordinate) {
    return (2. * coordinate - 1.) * border;
}


double known_eigenfunction(IndexDimension I) {
    // Calculate distance from zero r
    double r = 0.;
    for(int d=0; d<DimensionSparseGrid; d++){
        const double atom_coord = unitSQ2atom(I.coordinate(d));
        r += atom_coord * atom_coord;
    }
    r = sqrt(r);

    double val = 2.0* exp(- 1.0 * r);// / b_factor;


    return val;
}


double laplace_sin(IndexDimension I) {
    double val = 1.;
    for (int d=0; d<DimensionSparseGrid; d++) {
        val *= sin(M_PI * I.coordinate(d));
    }
    return val;
}






int main(int argc, char **argv) {

    unsigned int seed = std::random_device{}(); // Fallback auf Zufall
    seed = 649454244;
    


    

    int rank = 0;
    int num_tasks = 1;


    #ifdef MY_MPI_ON
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    #endif



    mpi_cout(" Dimension " + to_string( DimensionSparseGrid));
    mpi_cout("regular grid, solve: -lap u + (c-" + to_string(1.0+shiftfactor) + "lambda)*u =  -" + to_string(1.0/shiftfactor) + "lamda");


    #pragma omp parallel
        {
    #pragma omp single
            {
                int total_threads = omp_get_num_threads();
                mpi_cout("number of mpi_tasks " + to_string(num_tasks));
                mpi_cout("number of omp_threads " + to_string(total_threads));
            }
        }


    int numberMVprocesses=2;
    int numberLSprocesses=num_tasks - numberMVprocesses;

    mpi_cout("use " + to_string(numberLSprocesses)+ " for LocalStiffnessmatrices ");
    mpi_cout("use " + to_string(num_tasks - numberLSprocesses)+ " for Matrix Vector Multiplication ");


    int level_start=1;
    int level_end=10;
    


    if (rank == 0) {
        cout << "Computing the first eigenvalue of the Schroedinger equation for the Spinless Hydrogen atom." << endl;
        cout << "The expected eigenvalue is " << eig_exp << endl << endl;
        cout << "Dimension " << DimensionSparseGrid << endl;
        cout << "Computing on domain [-" << border << ", " << border << "]" << endl;
        cout << "With epsilon = " << eps << endl;
        cout << "With cg_epsilon = " << cg_eps << endl;
        cout << "Coefficient cutoff = " << CUTCOEFF << endl;
        cout << "Shift factor = " << shiftfactor << endl;
        cout << "Monte Carlo samples and ccenter=10 * mccsamples " << mc_samples << endl << endl;
        cout << " Monte Carlo seed " << seed  << endl;
        cout << "Adaptive CG" << endl;
        cout << "level start " << level_start << endl;
        cout << "level end " << level_end << endl;
        //mpi_cout("level\tDOFS\tSummed\t\tPoisson\t\tVarC\t\t\tDiffADD");

        cout << "Factor k^2 = " << varc_factor << endl << endl;
        cout << "ev without backshift and transformation " << shiftfactor*0.5*rhs_factor << endl; ;

        string legend = "level\tDOFS\t\t\tHelmholz";

        legend += "\t\tpow_iter\tcg_iter";
     
        cout << legend << endl << endl;
        //cout <<"level\tDOFS\tByMult\t\tPower" << endl;
    }



    





    AdaptiveSparseGrid old_grid;
    IndexDimension centerPoint;
    old_grid.AddRecursiveSonsOfPoint(centerPoint,3);


    double ev_old=1.0;
    double err_rel_old= 1.0;

    double Linfty_old = 1.0;
    string err_str = "";
    string err_rel_str = "";
        

    for (int level = level_start; level <level_end+1; level++)
    {
        // Start measuring time
        struct timeval begin_all, end_all;
        gettimeofday(&begin_all, 0);

        IndexDimension centerPoint;
        AdaptiveSparseGrid grid;
        //grid.copy_inner(old_grid);
        grid.AddRecursiveSonsOfPoint(centerPoint,level);


        //grid.completeDirichletGrid();


        MultiLevelAdaptiveSparseGrid mgrid(&grid);
        MatrixVectorHomogen m(grid, mgrid, 1, 1);


        VectorSparseG u(grid);
        VectorSparseG prew_u(grid);
        VectorSparseG u_diff(grid);


        double val_u = 0.0;


        for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++) {
                IndexDimension I = grid.getIndexOfTable(i);
                val_u = known_eigenfunction(I);
                u.setValue(i,val_u);
        }
        calcPrewByNodal(prew_u, u);






        string mults_str = "";
        string powers_str = "";
        




            int pow_interations;
            double cg_interations;
            double time_power;
            double precon = 0;
            VectorSparseG prew_inv(grid);
            prew_inv = 1.0;

        
            Poisson poisson(grid);    
            HelmHoltz helm_holtz(grid);

            StencilMC<double (*)(double *)> stencilVarCoeff(grid,&var_coeff,mc_samples,seed);
            LocalStiffnessMatricesDynamicDistribution lhs(grid, stencilVarCoeff,numberLSprocesses);

            lhs.addStencil(poisson);

            LocalStiffnessMatricesDynamicDistribution rhs(grid, helm_holtz,numberLSprocesses);



  

            double eigenvalue_power;


            prew_inv = 1.0;

            Power::power_inverse_riccarda<LocalStiffnessMatricesDynamicDistribution, LocalStiffnessMatricesDynamicDistribution>(eps,cg_eps, prew_inv, eigenvalue_power, pow_interations, cg_interations, m, lhs, rhs, true, precon, time_power);

        VectorSparseG Ax(grid),Mx(grid),r(grid);
        m.multiplication(prew_inv,Ax,lhs);
        m.multiplication(prew_inv,Mx,rhs);
        r = Ax - eigenvalue_power*Mx;
        double res = product(r,r);





            double ev=eigenvalue_power/(rhs_factor*shiftfactor);
            double ev2=eigenvalue_power;
            powers_str += "\t" + to_string(eigenvalue_power/(rhs_factor*shiftfactor));
            double err= abs(0.5-ev);
            double err_rel = abs(err)/0.5;
        
      

            powers_str += "\t\t" + to_string(pow_interations) + "\t error "+to_string(err)+"\t"+to_string(ev_old/err)+ "\terr_rel " + to_string(err_rel) +"\t"+to_string(err_rel_old/err_rel);
            ev_old=err;
            err_rel_old = err_rel;


            

        
        //mpi_cout(to_string(level) + "\t" + to_string(grid.getDOFS()) + "\t", false);
        //mpi_cout(powers_str,false);
        
// 2. Warten, bis alle Prozesse den Konsolen-Output abgesetzt haben
//MPI_Barrier(MPI_COMM_WORLD);

        if(rank==0){

            cout << level << "\t" << grid.getDOFS() << "\t" << powers_str << endl;
            const char* raw_id = std::getenv("SLURM_JOB_ID");
            std::string job_id = (raw_id != nullptr) ? raw_id : "local";
            std::string filename = "temp_run_" + job_id + ".dat";
    
            // Append-Modus: Schreibt Level 1 bis 10 nacheinander in dieses File
            std::ofstream file(filename, std::ios::app);
       
    
        if (file.is_open()) {
                 if (level == 1) {
                        file << "\n"; // Die Leerzeile vor einem neuen Block
                    }

                file << job_id << "\t" << grid.getDOFS() << "\t" << err<< "\t"  << err_rel  << "\n";
            file.close();
            }
        }
   


    }


    

#ifdef MY_MPI_ON
    MPI_Finalize();
#endif
    return 0;


}
