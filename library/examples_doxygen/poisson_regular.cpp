//! [start]
#include "source/sparseEXPDE.h"


double sin_func(IndexDimension I) {
    double val = 1.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        val*=sin(M_PI*I.coordinate(d));
    }
    return val;
}


double laplace_sin_func(IndexDimension I) {
    double dim = double(DimensionSparseGrid);
    return dim*M_PI*M_PI*sin_func(I);
}

//! [start]

//! [main]
int main(int argc, char **argv) {



    cout << " Dimension " << to_string( DimensionSparseGrid)<< endl;
    cout << " regular grid,constant coefficient, solve: -lap u = f " << endl;




    double Linfty_old = 1.0;
    double L2_old = 1.0;
    double H1_old =1.0;

    IndexDimension centerPoint;


    for (int level = 1; level < 10; level++)
    {


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
              double val_f = laplace_sin_func(I);


              u.setValue(i, val_u);
              f_vec.setValue(i, val_f);
              f_test.setValue(i, val_f);




            }

        }



        VectorSparseG cu(grid);

        MatrixVectorHomogen m(grid,mgrid,true,false);

        HelmHoltz helmHoltz(grid);
        Poisson stencilPoisson(grid);

        VectorSparseG prew(grid);
        prew = 0.0;
        

        
        calcPrewByNodal(prew,f_vec);




        cu = 0.0;





        prew = 0.0;


        calcPrewByNodal(prew,f_test);
        f_test=0.0;
        m.multiplication(prew,f_test,helmHoltz);



        f_vec = f_test+cu;


        int iterations;
        prew = 0.0;

        double time_p;
        CG::solveHomogen(1e-10, prew, f_vec, &iterations, u, m, stencilPoisson, &time_p);





        grid.WorkOnHangingNodes=false;

        VectorSparseG w(grid),uSolution(grid),d(grid);
        w=0.0;
        calcNodalByPrew(prew, w);
        uSolution = w;

        //! [main]

        //! [error]

        //calculate error
        d = u - uSolution;
        double linfty= L_infty(d);


        L2 L2_(grid,mgrid);
        double l2=L2_.getValue(d);

        H1 h1(grid,mgrid);
        double h1_ = h1.getValue(d);
        h1_ =l2+h1_;


        cout << level << " " << grid.getDOFS() << "  " << linfty << " " << Linfty_old/linfty << "  " <<  l2 << " " << L2_old/l2 <<   "    " <<  h1_ << " " <<H1_old/h1_<< "  " << to_string(iterations) <<  endl;


        Linfty_old = linfty;
        L2_old = l2;
        H1_old=h1_;



    }

    return 0;


}
//! [error]