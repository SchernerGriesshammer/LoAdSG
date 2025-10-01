//
// Created by scherner on 26.07.21.
//

#include "testing.h"
#include "../fullgrid/MatrixVectorMultiplication_Fullgrid.h"
#include "../stencils/Stencil.h"
#include "../stencils/MassStencil.h"
#include "../applications/norms.h"


void checkgrid(AdaptiveSparseGrid &sgrid, MultiLevelAdaptiveSparseGrid &mgrid) {

    MatrixVectorInhomogen matrix(sgrid,mgrid);


    int maxlevel = 0;
    for (int d = 0; d < DimensionSparseGrid; d++) {
        if (maxlevel < sgrid.getMaxDepth(d))
            maxlevel = sgrid.getMaxDepth(d);

    }

    // only grid points with coordinate 0<= x <= 1 should be added

    for (unsigned long k = 0; k < sgrid.getMaximalOccupiedSecondTable(); k++) {
        IndexDimension I = sgrid.getIndexOfTable(k);

        for (int d = 0; d < DimensionSparseGrid; d++) {
            int index = I.getIndex(d);
            if (Static_IndexOneD::coordinate(index) < 0 || Static_IndexOneD::coordinate(index) > 1) {
                cout << "error coordinate is not in [0,1]" << endl;
                I.PrintCoord();
            }


        }
    }


    // only one Index per grid point should be stored. (Important for adaptivity)
    for (unsigned long k = 0; k < sgrid.getMaximalOccupiedSecondTable(); k++) {
        IndexDimension I = sgrid.getIndexOfTable(k);
        int count = 0;
        for (unsigned long j = 0; j < sgrid.getMaximalOccupiedSecondTable(); j++) {
            IndexDimension Iv = sgrid.getIndexOfTable(j);
            if (I == Iv) count++;

        }
        if (count > 1) {
            cout << " error more than one Index for a grid point is stored " << endl;
            I.PrintCoord();
        }

    }


    // create fullgrid
    AdaptiveSparseGrid grid;
    IndexDimension centerPoint;
    grid.AddPoint(centerPoint);

    fullgrid(maxlevel, grid);

   // grid.completeGrid();

    VectorSparseG sprew(sgrid);
    VectorSparseG sprew_vergleich(sgrid);
    VectorSparseG sAx(sgrid);
    VectorSparseG sAx_vergleich(sgrid);


    VectorSparseG Ax(grid);
    VectorSparseG prew(grid);


    VectorSparseG diff(sgrid);


    sprew = 1.0;

    for (int d = 0; d < DimensionSparseGrid; d++) {
        Coordinate Xi(sgrid, d);
        sprew = sprew * Xi;
    }



    PoissonStencil poissonStencil;
    matrix.multiplication<PoissonStencil>(sprew, sAx, poissonStencil);

    prew = sprew;

    //matrixmult_prew_fullgrid_inhomog(prew, Ax);
    sgrid.WorkOnHangingNodes = false;
    sAx_vergleich = Ax;


    diff = sAx_vergleich - sAx;

    if (L_infty(diff) > 1e-5) {
        cout << "error: different results for sparse grid and fullgrid" << endl;

        cout << L_infty(diff) << endl;


        cout << "fullgrid " << endl;
        Ax.PrintDouble(maxlevel);
        sAx_vergleich.PrintDouble(maxlevel);

        cout << "sparsegrid " << endl;
        sAx.PrintDouble(maxlevel);
        diff.PrintDouble(maxlevel);

        exit(1);


    }


    VectorSparseG su(sgrid);
    CalcUbyPrew_Neumann(sprew, su);
    CalcPrew_inhomogen(sprew_vergleich, su);

    diff = sprew - sprew_vergleich;
    if (L_infty(diff) > 1e-5) {
        cout << "error: basis transformation error" << endl;
        cout << L_infty(diff) << endl;

    }
}


void checkgrid_dirichlet(AdaptiveSparseGrid &sgrid, MultiLevelAdaptiveSparseGrid &mgrid) {
    cout << "checkgrid " << endl;


    int maxlevel = 0;
    for (int d = 0; d < DimensionSparseGrid; d++) {
        if (maxlevel < sgrid.getMaxDepth(d))maxlevel = sgrid.getMaxDepth(d);

    }

    for (unsigned long k = 0; k < sgrid.getMaximalOccupiedSecondTable(); k++) {
        IndexDimension I = sgrid.getIndexOfTable(k);

        for (int d = 0; d < DimensionSparseGrid; d++) {
            int index = I.getIndex(d);
            if (Static_IndexOneD::coordinate(index) < 0 || Static_IndexOneD::coordinate(index) > 1) {
                cout << "error " << endl;
                I.PrintCoord();
            }


        }
    }
    for (unsigned long k = 0; k < sgrid.getMaximalOccupiedSecondTable(); k++) {
        IndexDimension I = sgrid.getIndexOfTable(k);
        int count = 0;
        for (unsigned long j = 0; j < sgrid.getMaximalOccupiedSecondTable(); j++) {
            IndexDimension Iv = sgrid.getIndexOfTable(j);
            if (I == Iv) count++;

        }
        if (count > 1) {
            cout << " error " << endl;
            I.PrintCoord();
        }

    }


    AdaptiveSparseGrid grid;
    IndexDimension centerPoint;

    grid.AddRecursiveSonsOfPoint(centerPoint, maxlevel+3);
    MultiLevelAdaptiveSparseGrid mgrid2(&grid);

    MatrixVectorHomogen matrix2(grid, mgrid2,0,0);



    VectorSparseG sprew(sgrid);
    VectorSparseG sprew_vergleich(sgrid);
    VectorSparseG sAx(sgrid);

    VectorSparseG sAx_vergleich(sgrid);
    VectorSparseG Ax(grid);
    VectorSparseG prew(grid);

    //sgrid.PrintActiveHanging(3);



    VectorSparseG diff(sgrid);



    sprew = 1.0;



    prew = sprew;

    MatrixVectorHomogen matrix(sgrid, mgrid,0,0);

    //PoissonStencil stencil;
    MassStencil massStencil;
    //StencilTemplate massStencil2(&sgrid);


    Poisson stencil(sgrid);

    PoissonStencil stencil3;

    //matrix.multiplication<StencilTemplate>(sprew, sAx, massStencil2);
    matrix.multiplication(sprew,sAx,stencil);


    //matrix.multiplication2<StencilTemplate>(prew, Ax, stencil);
    //matrixmult_prew_fullgrid(prew, Ax, Stiffness);
    cout << " ----------------------- " << endl;
    //matrix2.multiplication(prew, Ax, stencil3);
    //matrix2.multiplication<MassStencil>(prew, Ax, massStencil);

    sgrid.WorkOnHangingNodes = false;
    grid.WorkOnHangingNodes=false;
    sAx_vergleich = 0.0;
    sAx_vergleich = Ax;
    diff = sAx_vergleich - sAx;
    cout << "ERROR SPARSE-FULL " << L_infty(diff) << endl;
    if (L_infty(diff) > 1e-10) {
        cout << "error in matrix-vector-multiplication check dirichlet_grid " << endl;
        sAx_vergleich.PrintDouble(3);
        sAx.PrintDouble(3);
        exit(1);
    }


    VectorSparseG su(sgrid);
    calcNodalByPrew(sprew, su);
    calcPrewByNodal(sprew_vergleich, su);

    diff = sprew - sprew_vergleich;
    if (L_infty(diff) > 1e-5) {
        cout << "error in prewavelet calculation check grid " << endl;
        cout << L_infty(diff) << endl;

    }


}



