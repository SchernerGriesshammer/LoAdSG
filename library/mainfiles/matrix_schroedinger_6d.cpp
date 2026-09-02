//
// Created by rumi on 19/01/26.
//

#include "../source/sparseEXPDE.h"

//#define testNOW
//#define goolePerf
#include <string>
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <iomanip>


double border = 7.5;
double eps = 1e-6;
const int lvl = 4;

const int count = 3;

const double CUTCOEFFarray[] = {1e5};
double CUTCOEFF;
const int mcsamplesarray[] = {200};
int mc_sample;
const double shiftfactorarray[] = {5.};


const size_t singleDimension = 3;
const double shiftfactor = 5.0;
int num_electrons = 2;

    /// @brief Factor for variable coefficient
    /// @details 2: for factor before Helmholtz operator in atomic units
    /// 4*border: for coordinate transform factor
    const double varc_factor = 2*4*border;

    /// @brief factor for scaling eigenvalue
    /// @details from coordinate transform and Helmholtz factor
    const double rhs_factor = 2*4*border*border;

    /// @brief expected eigenvalue for one core and two electrons in atomic units
    const double eig_exp = -2.90338583;
    const double transf_eig_exp = eig_exp*rhs_factor;
    double shiftvalue;

void set_shiftvalue(double newfactor) {
    shiftvalue = (newfactor+1) * (-1) * transf_eig_exp;
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
    if (r!=0) {
        result = min(1./r, CUTCOEFF);
    }
    // log sampled values
    r_count(1.0/r);

    return -2.0 * result;
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
    if (r != 0) {
        result = min(1./r, CUTCOEFF);
    }
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


    int numberMVprocesses=16;
    int numberLSprocesses= num_tasks - numberMVprocesses;

    mpi_cout("use " + to_string(numberLSprocesses)+ " for LocalStiffnessmatrices ");
    mpi_cout("use " + to_string(numberMVprocesses)+ " for Matrix Vector Multiplication ");





    if (rank == 0) {
        cout << endl;
        cout << "Computing the first eigenvalue of the Schroedinger equation for the HELIUM atom, i.e. with two electrons." << endl;
        cout << "The expected eigenvalue is " << eig_exp << endl;
        cout << "Computing on domain [-" << border << ", " << border << "]" << endl;
        cout << "With epsilon = " << eps << endl;
        cout << "On level l = " << lvl << endl;
        cout << "Using the shift trick for positive eigenvalue with shift factors [";
            for(double sv : shiftfactorarray) {cout << to_string(sv) << ", ";}
            cout << "]." << endl;
        cout << "Coefficient cutoffs = [";
            for (double cc : CUTCOEFFarray) {cout << to_string(cc) << ", ";}
            cout << "]." << endl;
        cout << "Monte Carlo sampless = [";
            for (int mcs : mcsamplesarray) {cout << to_string(mcs) << ", ";}
            cout << "]." << endl;
        cout << "shiftfactor " << shiftfactor << endl;

        cout << endl << endl;


    }




    // Regular grid construction
    AdaptiveSparseGrid grid;
    IndexDimension centerPoint;
    grid.AddRecursiveSonsOfPoint(centerPoint,lvl);

    MultiLevelAdaptiveSparseGrid mgrid(&grid);
    MatrixVectorHomogen m(grid, mgrid, 1, 1);



    set_shiftvalue(shiftfactor);


                // Matrix construction
    Poisson poisson(grid);
    HelmHoltz rhs(grid);
    StencilMC<double (*)(double *)> stencilVarCoeff(grid,&var_coeff,100);
    LocalStiffnessMatricesDynamicDistribution lhs(grid, stencilVarCoeff,numberLSprocesses);
    lhs.addStencil(poisson);


    // fields for solving
   VectorSparseG x(grid),Ax(grid),Bx(grid);

    int size =grid.getDOFS();

    double mAx[size][size];
    double mBx[size][size];

for (unsigned long k=0; k<grid.getMaximalOccupiedSecondTable();k++)
    {
        IndexDimension p=grid.getIndexOfTable(k);
        x.setValue(p,1.0);

        m.multiplication(x,Ax,lhs);
        for (unsigned long j=0;j< grid.getMaximalOccupiedSecondTable();j++)
        {
            IndexDimension q=grid.getIndexOfTable(j);
            mAx[k][j]=Ax.getValue(q);

        }

        x=0.0;
        Ax=0.0;

    }


    for (unsigned long k=0; k<grid.getMaximalOccupiedSecondTable();k++)
    {
        IndexDimension p=grid.getIndexOfTable(k);
        x.setValue(p,1.0);

        m.multiplication(x,Bx,rhs);
        for (unsigned long j=0;j< grid.getMaximalOccupiedSecondTable();j++)
        {
            IndexDimension q=grid.getIndexOfTable(j);

            mBx[k][j]=Bx.getValue(q);
        }

        x=0.0;
        Bx=0.0;

    }

if (rank==0)
{
    std::ofstream datei("mA4.txt");

    if (datei.is_open()) {
        // Hohe Präzision sicherstellen, damit double-Werte nicht gerundet werden
        datei << std::setprecision(15);

        for (int i = 0; i < grid.getMaximalOccupiedSecondTable(); ++i) {
            for (int j = 0; j < grid.getMaximalOccupiedSecondTable(); ++j) {
                datei << mAx[i][j];
                if (j < grid.getMaximalOccupiedSecondTable() - 1) datei << " "; // Leerzeichen zwischen Spalten
            }
            datei << "\n"; // Zeilenumbruch nach jeder Zeile
        }
        datei.close();
        std::cout << "Matrix  A erfolgreich gespeichert." << std::endl;
    }



    std::ofstream datei_B("mB4.txt");

    if (datei_B.is_open()) {
        // Hohe Präzision sicherstellen, damit double-Werte nicht gerundet werden
        datei_B << std::setprecision(15);

        for (int i = 0; i < grid.getMaximalOccupiedSecondTable(); ++i) {
            for (int j = 0; j < grid.getMaximalOccupiedSecondTable(); ++j) {
                datei_B << mBx[i][j];
                if (j < grid.getMaximalOccupiedSecondTable() - 1) datei_B << " "; // Leerzeichen zwischen Spalten
            }
            datei_B << "\n"; // Zeilenumbruch nach jeder Zeile
        }
        datei_B.close();
        std::cout << "Matrix B erfolgreich gespeichert." << std::endl;
    }
}



#ifdef MY_MPI_ON
    MPI_Finalize();
#endif

    return 0;


}