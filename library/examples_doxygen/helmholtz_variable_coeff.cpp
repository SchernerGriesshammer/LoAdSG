#include "source/sparseEXPDE.h"

double sin_func(const IndexDimension I) {

    double val = 1.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        val*=sin(M_PI*I.coordinate(d));
    }
    return val;

}


double laplace_sin_func(const IndexDimension I) {
    constexpr double dim = static_cast<double>(DimensionSparseGrid);
    return dim*M_PI*M_PI*sin_func(I);
}


void mpi_cout(const string &s){
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


double var_coeff( double* x)
{
    const double sum = (x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2]);
    return 1.0/sum;
}




int main(int argc, char **argv) {
    int rank = 0;
    int num_tasks = 1;


#ifdef MY_MPI_ON
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif



    mpi_cout(" Dimension " + to_string( DimensionSparseGrid));
    mpi_cout("regular grid, solve: -lap u + c*u = f ");


    int level_start=1;




#pragma omp parallel
    {
#pragma omp single
        {
            int total_threads = omp_get_num_threads();
            mpi_cout("number of mpi_tasks " + to_string(num_tasks));
            mpi_cout("number of omp_threads " + to_string(total_threads));
        }
    }


    int numberMVprocesses=1;
    int numberLSprocesses=num_tasks - numberMVprocesses;

    mpi_cout("use " + to_string(numberLSprocesses)+ " for LocalStiffnessmatrices ");
    mpi_cout("use " + to_string(num_tasks - numberLSprocesses)+ " for Matrix Vector Multiplication ");


    double Linfty_old = 1.0;
    double L2_old = 1.0;
    double H1_old =1.0;



    IndexDimension centerPoint;



    mpi_cout(" varCoeff  1/x^2 ");

    for (int level = level_start; level <10; level++) {
        AdaptiveSparseGrid grid;
        grid.AddRecursiveSonsOfPoint(centerPoint,level);

        MultiLevelAdaptiveSparseGrid mgrid(&grid);



        VectorSparseG u(grid), f_vec(grid), f_test(grid),fpu(grid);
        u=0.0;




        for (unsigned long i = 0; i < grid.getMaximalOccupiedSecondTable(); i++) {
            if(grid.getActiveTable()[i]) {
                IndexDimension I = grid.getIndexOfTable(i);
                double x[DimensionSparseGrid];
                for (int d = 0; d < DimensionSparseGrid; d++) x[d] = I.coordinate(d);


                double val_u = sin_func(I);
                double val_f= laplace_sin_func(I);


                u.setValue(i, val_u);
                f_vec.setValue(i, val_f);
                f_test.setValue(i, val_f);
            }






            }



            VectorSparseG cu(grid);

            MatrixVectorHomogen m(grid,mgrid,true,true);

            HelmHoltz helmHoltz(grid);
            Poisson stencilPoisson(grid);
            VectorSparseG prew(grid);


            prew = 0.0;
            calcPrewByNodal(prew,f_vec);



            StencilMC<double (*)(double *)> stencilVarCoeff(grid,&var_coeff,1);
            LocalStiffnessMatricesDynamicDistribution localStiffnessMatrices_VarCoeff(grid, stencilVarCoeff,numberLSprocesses);





            cu = 0.0;
            prew=0.0;

            calcPrewByNodal(prew,u);
            m.multiplication(prew,cu,localStiffnessMatrices_VarCoeff);


            prew = 0.0;
            calcPrewByNodal(prew,f_test);

            f_test=0.0;
            m.multiplication(prew,f_test,helmHoltz);



            f_vec = f_test+cu;


            int iterations;
            prew = 0.0;

            double time_p;

            localStiffnessMatrices_VarCoeff.addStencil(stencilPoisson);



            CG::solveHomogen(1e-10, prew, f_vec, &iterations, u, m, localStiffnessMatrices_VarCoeff, &time_p);


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




            if(rank==0){
                cout <<level << " " << grid.getDOFS() << "  " << linfty << " " << Linfty_old/linfty << "  " <<  l2 << " " << L2_old/l2 <<   "    " <<  h1_ << " " <<H1_old/h1_<< "  " << to_string(iterations) << " " << endl;
            }

            Linfty_old = linfty;
            L2_old = l2;
            H1_old=h1_;


            prew=0.0;
            u=0.0;
        }

#ifdef MY_MPI_ON
        MPI_Finalize();
#endif
        return 0;



}