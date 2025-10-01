//
// Created by to35jepo on 3/20/23.
//

#include "MatrixVectorMultiplication_Fullgrid.h"
#include "../BasisTransformations/BasisTransformations_inhomog.h"
#include "../stencils/PoissonStencil.h"
#include "../stencils/Stencil.h"



void restriction_full(VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {
    AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    Depth Tzero(0);
    int tneu = t;


    for (unsigned long i = 0; i < maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);

        if (Index.getDepth(d) <= t) {
            double value = fine.getValue(Index);
            if (Index.isAtLeftBoundary(d))
                value = fine.getValue(Index) + 0.5 * fine.getValue(Index.nextRight(d, tneu + 1));

            if (Index.isAtRightBoundary(d))
                value = fine.getValue(Index) + 0.5 * fine.getValue(Index.nextLeft(d, tneu + 1));

            if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                value = fine.getValue(Index) + 0.5 * fine.getValue(Index.nextRight(d, tneu + 1))
                        + 0.5 * fine.getValue(Index.nextLeft(d, tneu + 1));
            //coarse.setValue(Index, value);
            coarse = value | Index;


        }


    }

}

void restriction_local_boundary(VectorSparseG &fine, VectorSparseG &coarse, Depth Tfine, Depth Tcoarse) {
    coarse = fine;

    Depth Tvergleich(0);;


    for (int d = 0; d < DimensionSparseGrid; d++) {
        for (int t = Tfine.at(d) - 1; t >= int(Tcoarse.at(d)); --t) {

            //restriction3(coarse, t, d);
            restriction_full(coarse, coarse, t, d);


        }
    }
}

void ApplyStencil_Neumann(VectorSparseG &x, VectorSparseG &Ax) {
    PoissonStencil poissonStencil;

    AdaptiveSparseGrid_Base *grid = x.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();


    for (unsigned long i = 0; i < maxocc; i++) {

        IndexDimension Index = grid->getIndexOfTable(i);
        Depth T(0);


        for (int d = 0; d < DimensionSparseGrid; d++) {
            int t = x.getSparseGrid()->getMaxDepth(d);
            T.set(t, d);
        }


        double val = getLocalStencilBoundary<PoissonStencil>(Index, T, x,poissonStencil);

        Ax.setValue(i,val);


    }
}

/*
void matrixmult_prew_fullgrid_inhomog(VectorSparseG &x, VectorSparseG &Ax) {
    AdaptiveSparseGrid *grid = x.getSparseGrid();


    ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outer(list);
    outer.gotoEnd();
    Depth T = outer.getSubgrid()->getT();



    VectorSparseG nodal(*grid);


    VectorSparseG finest_nodal(*grid);
    VectorSparseG Ax_hier(*grid);


    CalcUbyPrew_Neumann(x, nodal);


    ApplyStencil_Neumann(nodal, finest_nodal);


    ListOfDepthOrderedSubgrids::iterator inneriter(list);

    do {

        Depth Tcoarse = outer.getSubgrid()->getT();


        restriction_local_boundary(finest_nodal, Ax_hier, T, Tcoarse);


        ConvertToNeumannPrewavelet(Ax_hier, Ax, Tcoarse);


        Ax_hier = 0.0;
    } while (outer.previous());


}*/
