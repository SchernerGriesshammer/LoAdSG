//
// Created by rumi on 19/01/26.
//

#include "../source/sparseEXPDE.h"
#include "../source/iterator/depthIterator.h"
#include "../source/applications/power_method.h"
//#define testNOW
//#define goolePerf
#include <string>
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <fstream>



double border = 7.5;
double eps = 1e-5;
double cg_eps = 1e-12;

const int count = 1;
double CUTCOEFF = 1e4;
int mc_sample = 400;
double shiftfactor = 5;

const size_t singleDimension = 3;
int num_electrons = 2;

    /// @brief Factor for variable coefficient
    /// @details 2: for factor before Helmholtz operator in atomic units
    /// 4*border: for coordinate transform factor
    const double varc_factor = 2*4*border;

    /// @brief factor for scaling eigenvalue
    /// @details from coordinate transform and Helmholtz factor
    const double rhs_factor = 2*4*border*border;

    /// @brief expected eigenvalue for one core and two electrons in atomic units
    /// 2.9037 2437 7034 1195 9830
    const double eig_exp = -2.90338583;
    const double transf_eig_exp = eig_exp*rhs_factor;
    double shiftvalue;

void set_shiftvalue(double newfactor) {
    shiftvalue = (newfactor+1.0) * (-1.0) * transf_eig_exp;
}




/*
const int num_bins = 100;
vector<double> create_binborders() {
    vector<double> binborders(num_bins+1);
    binborders[num_bins] = INFINITY;
    for (size_t i = 0; i < num_bins; i++)
    {
        binborders[i] = 1. - (log((i+1)/num_bins) / log(1./num_bins));
    }

    return binborders;
}
*/


// logging sample values
double r_max = 1;
double r_min = 1;
/*
vector<double> binborders_norm = create_binborders();
vector<double> binborders(binborders_norm.size());
vector<int> bincontents(binborders.size() - 1);

void set_binborders(double newcutoff) {
    for (size_t i = 0; i < binborders.size(); i++)
    {
        binborders[i] = newcutoff * binborders_norm[i];
    }
}
*/



/// @brief mpi appropriate output
/// @param s output string
/// @param endline whether to end with linebreak
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



/// @brief gather logged minimum and maximum mc sampled values from all mpi threads and write output string
/// @return logging string
string minmax_sample_combine_log() {
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

    // Reset all logging values
    r_min = 1;
    r_max = 1;
    return str;
}


/*
/// @brief gather logged mc sampled histogram values from all mpi threads and write output string
/// @return logging string
string hist_sample_combine_log() {
    string row1 = "";
    string row2 = "";
    
    int rank = 0;
    #ifdef MY_MPI_ON
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int sum[num_bins];
        MPI_Reduce((const int *)&sum, &bincontents[0], bincontents.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank==0) {
            for (size_t i = 0; i < num_bins; i++) {
                row1 += to_string(binborders[i]) + ";";
                row2 += to_string(bincontents[i]) + ";";
            }
            row1 += "\n";
            row2 += "\n\n";
        }
    #endif

    // Reset bin contents
    fill(bincontents.begin(), bincontents.end(), 0);
    return row1 + row2;
}
*/



void inline r_count(double over_r) {
    if (over_r < r_min) r_min = over_r;
    if (!isinf(over_r) && over_r > r_max) r_max = over_r;

    /*
    // binary search to find right bin
    int leftidx = 0;
    int rightidx = binborders.size()-1;
    while(!(binborders[leftidx] <= over_r && binborders[leftidx+1] > over_r)) {
        int middleidx = int((rightidx - leftidx) / 2);
        if (binborders[middleidx] <= over_r) {
            leftidx = middleidx;
        } else {
            rightidx = middleidx;
        }
    }

    // count up respective bin
    bincontents[leftidx] += 1;
    */
}







/// @brief variable coefficient term for electron of given index and a core
/// @param coordinates coordinates for function evaluation
/// @param electron_idx index of the electron in question
/// @return value at coordinates
double electron_core(double* coordinates, size_t electron_idx) {
    double r = 0;
    for (int i = singleDimension*electron_idx; i < singleDimension*(electron_idx+1); i++) {
        double ci = 2. *  coordinates[i] - 1; // < map from [0,1] to [-s,s]
        r += ci * ci;
    }
    r = sqrt(r);

    double result  = CUTCOEFF;
    if (r > 1e-15){
        result = min(1. / r, CUTCOEFF); // < clipping large values of 1/r around the center
    };
    
    // log sampled values
    r_count(1/r);

    return -2. * result;
    // Still to happen in : varc_factor for coordinate transformation
    // trick for positive definiteness
}


/// @brief variable coefficient term for two electrons of given indices
/// @param coordinates coordinates for function evaluation
/// @param e_idx0 index of first electron
/// @param e_idx1 index of second electron
/// @return value at coordinates
double electron_electron(double* coordinates, size_t e_idx0, size_t e_idx1) {
    double r = 0;
    int offset0 = singleDimension*e_idx0;
    int offset1 = singleDimension*e_idx1;
    for(int i = 0; i < singleDimension; i++) {
        double ci = coordinates[offset0 + i] - coordinates[offset1 +  i];
        r += ci*ci;
    }
    r = sqrt(r);

    double result = CUTCOEFF;

    if (r > 1e-15){
        result = min(1. / r, CUTCOEFF); // < clipping large values of 1/r around the center
    };
    
    // log sampled values
    r_count(1/r);

    return 0.5*result;
    // Still to happen: varc_factor for coordinate transformation
    // trick for positive definiteness
}


/// @brief Calculates variable coefficient for two electrons and one core
/// @param coordinates coordinates for function evaluation
/// @return value at coordinates
double var_coeff(double* coordinates)
{
    double result = 0.;

    for (int e = 0; e < num_electrons; e++) {
        // Electron - core interaction term
        result += electron_core(coordinates, e);

        for (int e2 = 0; e2 < e; e2++) {
            // Electron - electron interaction term
            result += electron_electron(coordinates, e, e2);
        }
    }

    // Factor from coordinate transformation
    result *= varc_factor;

    // Trick to make positive definite
    result += shiftvalue;

    return result;
}



int main(int argc, char **argv) {

    unsigned int seed = std::random_device{}(); // Seed für Monte Carlo Integration
    seed = 3219425942;


    int rank = 0;
    int num_tasks = 1;


    #ifdef MY_MPI_ON
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    #endif



    mpi_cout("Grid dimension " + to_string(DimensionSparseGrid));
    mpi_cout("Regular grid.");
    mpi_cout("");

    mpi_cout("Single particle dimension " + to_string(singleDimension));

    if (DimensionSparseGrid % singleDimension != 0)exit(1);


    #pragma omp parallel
        {
    #pragma omp single
            {
                int total_threads = omp_get_num_threads();
                mpi_cout("number of mpi_tasks " + to_string(num_tasks));
                mpi_cout("number of omp_threads " + to_string(total_threads));
            }
        }


    int numberMVprocesses=48;
    int numberLSprocesses= num_tasks - numberMVprocesses;

    mpi_cout("use " + to_string(numberLSprocesses)+ " for LocalStiffnessmatrices ");
    mpi_cout("use " + to_string(numberMVprocesses)+ " for Matrix Vector Multiplication ");




    if (rank == 0) {
        cout << endl;
        cout << "Computing the first eigenvalue of the Schroedinger equation for the HELIUM atom, i.e. with two electrons." << endl;
        cout << "The expected eigenvalue is " << eig_exp << endl;
        cout << "Computing on domain [-" << border << ", " << border << "]" << endl;
        cout << "With epsilon = " << eps << endl;
        cout << "With cg_epsilon = " << cg_eps << endl;
        cout << "Using the shift trick for positive eigenvalue with shift factor " << shiftfactor << endl;
        cout << "Coefficient cutoffcoeff "<< CUTCOEFF << endl;
        cout << "Monte Carlo samples " << mc_sample   << endl;
        cout << "Monte Carlo seed" << seed << endl;

        cout << endl << endl;

        string legend = "level\tDOFS\t\tsmcsampl\tclip\tshift\t\teig\t\teig_diff\t\tcg_iter\t\tpow_iter\ttime";
        cout << legend << endl << endl;


       
    }


    int level_start = 7;
    int level_end = 10;
    double err_old = 1.0;
    double err_rel_old = 1.0;
     for (int level = level_start; level <level_end+1; level++){

    // Regular grid construction
    AdaptiveSparseGrid grid;
    IndexDimension centerPoint;
    grid.AddRecursiveSonsOfPoint(centerPoint,level);

    MultiLevelAdaptiveSparseGrid mgrid(&grid);
    MatrixVectorHomogen m(grid, mgrid, 1, 1);


    set_shiftvalue(shiftfactor);


     // Matrix construction
     Poisson poisson(grid);
     HelmHoltz rhs(grid);
     StencilMC<double (*)(double *)> stencilVarCoeff(grid,&var_coeff,mc_sample,seed);
     LocalStiffnessMatricesDynamicDistribution lhs(grid, stencilVarCoeff,numberLSprocesses);
     lhs.addStencil(poisson);


      // fields for solving
      double eigenvalue_power;
      int pow_interations;
      double cg_interations;
      double time_power;
       double precon = 0;

      VectorSparseG prew_inv(grid);
      prew_inv = 1.0;

      // SOLVING
        Power::power_inverse_riccarda<LocalStiffnessMatricesDynamicDistribution, HelmHoltz>(eps,cg_eps, prew_inv, eigenvalue_power, pow_interations, cg_interations,
                m, lhs, rhs, true, precon, time_power);

                

                // Undoing the shift
                eigenvalue_power /= (rhs_factor * shiftfactor);
                double error = abs(eigenvalue_power + eig_exp);
                double err_rel = error / abs(eig_exp);

         VectorSparseG u_plot(grid);
         calcNodalByPrew(prew_inv,u_plot);

         if (rank ==0)
         {
             std::ofstream outfile("raw_sparse_grid_density.dat");
             outfile << "r1,r2,Dm\n"; // Header für einfaches Einlesen in Python
             for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++)
             {
                 if(grid.getActiveTable()[i])
                 {
                     IndexDimension I = grid.getIndexOfTable(i);
                     double x1,x2,x3;
                     double y1,y2,y3;

                     x1=I.coordinate(0);
                     x2=I.coordinate(1);
                     x3=I.coordinate(2);
                     y1=I.coordinate(3);
                     y1=I.coordinate(4);
                     y1=I.coordinate(5);

                     double u_val = u_plot.getValue(i);

                     // 4. Radiale Paardichte berechnen
                     double D_m = r1 * r1 * r2 * r2 * (u_val * u_val);

                     // 5. Einfach als rohe Punktwolke abspeichern
                     outfile << r1 << "," << r2 << "," << D_m << "\n";



                 }
             }
             outfile.close();
         }
         
             

            calcNodalByPrew(prew_inv,u_plot);



         if(rank == 0){
             cout << level << "\t" << grid.getDOFS() << "\t" << eigenvalue_power << "\t" << error << "\t" << err_old/error << "\t" << err_rel << "\t" << err_rel_old/err_rel << "\t" << pow_interations << endl;
                         const char* raw_id = std::getenv("SLURM_JOB_ID");
            std::string job_id = (raw_id != nullptr) ? raw_id : "local";
            std::string filename = "temp_run_" + job_id + ".dat";
    
            // Append-Modus: Schreibt Level 1 bis 10 nacheinander in dieses File
            std::ofstream file(filename, std::ios::app);
       
    
            if (file.is_open()) {
                 if (level == 1) {
                        file << "\n"; // Die Leerzeile vor einem neuen Block
                    }

                file << job_id << "\t" << grid.getDOFS() << "\t" << error<< "\t"  << err_rel  << "\n";
            file.close();
            }
        }
        err_old = error;
        err_rel_old = err_rel;

    }



             
         

               

     
            
  


#ifdef MY_MPI_ON
    MPI_Finalize();
#endif

    return 0;


}