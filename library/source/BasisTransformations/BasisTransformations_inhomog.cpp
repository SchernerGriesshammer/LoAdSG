//
// Created by scherner on 29.04.21.
//

#include "BasisTransformations_inhomog.h"
#include "../tests/testing.h"

#include "../iterator/RectangularIterator.h"


void CalcNodalBasis_Neumann(VectorSparseG &u, ListOfDepthOrderedSubgrids &orderedSubgrids) {


    ListOfDepthOrderedSubgrids list2(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(list2);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) >= 1) {
                u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
            }
        } while (iterSubgrids.next());
    }
}


void
CalcHierarchicalBasis_Neumann(VectorSparseG &u, const Depth &Tsubgrid, ListOfDepthOrderedSubgrids &orderedSubgrids) {

    ListOfDepthOrderedSubgrids list2(*u.getSparseGrid());

    ListOfDepthOrderedSubgrids::iterator iterSubgrids(list2);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) >= 1 && T <= Tsubgrid) {
                u = u - 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
            }
        } while (iterSubgrids.previous());
    }
}

void CalcHierarchicalBasis_Neumann(VectorSparseG &u, ListOfDepthOrderedSubgrids &orderedSubgrids) {

    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) > 1) {
                u = u - 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
            }
        } while (iterSubgrids.previous());
    }
}

void CalcHierbyPrew_Neumann(VectorSparseG &prew, VectorSparseG &hier, SubgridFixedDepth &subgrid,
                            ListOfDepthOrderedSubgrids &list) {
    Depth T = subgrid.getT();

    //ListOfDepthOrderedSubgrids list2(*hier.getSparseGrid());
    SubgridFixedDepth::iterator iter(*list.getSubgrid(T));




    iter.gotobegin();
    do {

        IndexDimension Index = iter.getPoint();

        double coeff = prew.getValue(Index);
        Depth Tneu(Index);

        for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
            double basis_coeff = 1.0;
            bool next = true;
            IndexDimension J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff,&next);


            if(next) {

                double val = hier.getValue(J) + coeff * basis_coeff;

                hier = val | J;
            }

        }




    } while (iter.next());



    CalcHierarchicalBasis_Neumann(hier, T, list);


}

void CalcUbyPrew_Neumann(VectorSparseG &prew, VectorSparseG &hier, SubgridFixedDepth &subgrid) {


    //ListOfDepthOrderedSubgrids list2(*hier.getSparseGrid());
    //SubgridFixedDepth::iterator iter(*list.getSubgrid(T));
    SubgridFixedDepth::iterator iter(subgrid);
    iter.gotobegin();
    do {

        IndexDimension Index = iter.getPoint();
        Depth Tneu(Index);

        double coeff = prew.getValue(Index);

        for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
            double basis_coeff = 1.0;
            IndexDimension J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff);


            double val = hier.getValue(J) + coeff * basis_coeff;


            //hier.setValue(J, val);
            hier = val | J;

        }


    } while (iter.next());


}

void CalcUbyPrew_Neumann(VectorSparseG &prew, VectorSparseG &u) {
    AdaptiveSparseGrid *grid = u.getSparseGrid();

    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*grid);
    VectorSparseG hier_new(*grid);

    u = 0.0;

    outerIterationsSubgrids.gotoEnd();
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();


        hier_new = 0.0;

        CalcHierbyPrew_Neumann(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), list);


        innerIterationsSubgrids.gotoEnd();
        do {
            Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();

            if (Tlocal <= T) {
                hier = hier + hier_new | *innerIterationsSubgrids.getSubgrid();
            }
        } while (innerIterationsSubgrids.previous());


        u = hier;
        CalcNodalBasis_Neumann(u, orderedSubgrids);

    } while (outerIterationsSubgrids.previous());

}

void SolveSave_Neumann(VectorSparseG *prew, VectorSparseG *u, ZusammenhangsKomponente &komponente) {


    Depth T = komponente.getDepth();
    AdaptiveSparseGrid_Base *grid = prew->getSparseGrid();

    for (int d = 0; d < DimensionSparseGrid; d++) {


        IndexDimension minI = komponente.getMinIndex();
        IndexDimension maxI = komponente.getMaxIndex();


        IndexDimension iterIndex = minI;
        unsigned long k;
        int size = 1;


        if (IndexOneD(minI.getIndex(d)) == IndexOneD(maxI.getIndex(d)) && minI.getDepth(d) < 1) size = 1;
        else
            while (grid->occupied(k, iterIndex.nextRight(d, T.at(d))) && !iterIndex.isAtRightBoundary(d)) {
                if (grid->workonindex(k)) {
                    size++;
                    iterIndex = iterIndex.nextRight(d, T.at(d));
                } else {
                    break;
                }

            }

        double p[size];
        double nodal[size];

        PrewaveletMatrixInhomogen M(size, komponente.getDepth(), prew->getSparseGrid());

        M.create(minI, d);



        maxI.replace(d, IndexOneD(minI.getIndex(d)));
        iterIndex = minI;

        for (RectangularIterator iter(minI, maxI, komponente.getDepth()); iter.goon(); ++iter) {
            iterIndex = iter.getIndex();

            for (int i = 0; i < size; i++) {
                nodal[i] = u->getValue(iterIndex);
                p[i] = 0.0;
                iterIndex = iterIndex.nextRight(d, T.at(d));

            }

            M.solve(p, nodal);

            iterIndex = iter.getIndex();
            for (int i = 0; i < size; i++) {
                if (Depth(iterIndex) == komponente.getDepth())
                    *prew = p[i] | iterIndex;//prew->setValue(iterIndex, p[i]);
                *u = p[i] | iterIndex;//u->setValue(iterIndex, p[i]);
                iterIndex = iterIndex.nextRight(d, T.at(d));
            }


        }


    }
}

void CalcPrew_inhomogen(VectorSparseG &prew, VectorSparseG &u) {

    AdaptiveSparseGrid *grid = u.getSparseGrid();

    VectorSparseG uold(*grid);
    uold = u;

    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    //ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);

    ListOfDepthOrderedSubgrids orderedSubgrids2(*grid);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids2);


    VectorSparseG hier(*grid);
    VectorSparseG hier_new(*grid);

    prew = 0.0;
    hier = u;
    //CalcHierarchicalBasis_Neumann(hier, orderedSubgrids);
    CalcHierarchicalBasis(hier);


    outerIterationsSubgrids.gotoEnd();
    do {

        Depth T = outerIterationsSubgrids.getSubgrid()->getT();


        // hier: Iterator über Zusammenhangskomp. direkt in subgrid
        ZusammenhangsKomponente_Neumann komponente(grid, T);


        unsigned long maxoccup = grid->getMaximalOccupiedSecondTable();
        unsigned long iNow_mapping=-1;

        for (unsigned long iNow = komponente.findAnicht(-1,iNow_mapping); iNow < maxoccup; iNow = komponente.findAnicht(iNow,iNow_mapping)) {
            IndexDimension IndexNow = komponente.StartSearchComponent(iNow);


            komponente.recursiveMarkArbeite(IndexNow, iNow);

            //hier_new = u;


            //hier_new = u | komponente
            hier_new = 0.0;
            IndexDimension minI = komponente.getMinIndex();
            IndexDimension maxI = komponente.getMaxIndex();

            for (RectangularIterator iter2(minI, maxI, T); iter2.goon(); ++iter2) {
                IndexDimension Index = iter2.getIndex();

                double val = u.getValue(Index);
                hier_new = val | Index;
                //hier_new.setValue(Index, u.getValue(Index));
            }

            SolveSave_Neumann(&prew, &hier_new, komponente);


        }

        hier_new = 0.0;

        CalcHierbyPrew_Neumann(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), orderedSubgrids);




        //ziehe hierarchischen Überschuss ab


        innerIterationsSubgrids.gotoBegin();
        do {
            Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();
            if (Tlocal <= T) {
                hier = hier - hier_new | *innerIterationsSubgrids.getSubgrid();
            }
        } while (innerIterationsSubgrids.next());


        u = hier;
        CalcNodalBasis_Neumann(u, orderedSubgrids);


    } while (outerIterationsSubgrids.previous());

    u = uold;

}

void CalcNodalBasis_inhomog(VectorSparseG &u, Depth &Tiefe) {


    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T <= Tiefe) {
                if (T.at(d) >= 1) {
                    u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                    //AddHierUberschuss(u, T, d);
                }
            }
        } while (iterSubgrids.next());
    }

};


void CalcNodalBasis_inhomog(VectorSparseG &u) {


    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();

            if (T.at(d) > 1) {
                u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //AddHierUberschuss(u, T, d);
            }

        } while (iterSubgrids.next());
    }

};


void CalcUbyPrewRestrictions_inhomog(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions) {

    AdaptiveSparseGrid *grid = u.getSparseGrid();

    VectorSparseG test(*grid);
    // grid.getList();
    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*grid);
    VectorSparseG hier_new(*grid);

    u = 0.0;

    Depth Tv(Tiefe);
    for (int d = 0; d < DimensionSparseGrid; d++) {
        if (restrictions[d])Tv.set(Tiefe.at(d) - 1, d);
    }
    outerIterationsSubgrids.gotoEnd();
    unsigned long lonenorm = outerIterationsSubgrids.getDepth().LoneNorm();

    int maxnorm = 0;
    do {
        maxnorm = 0.0;
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();


        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (!restrictions[d] && !(T.at(d) <= Tiefe.at(d))) {

                goto goon;
            }
            if (restrictions[d] && T.at(d) != Tiefe.at(d)) {

                goto goon;
            }
            // if at least two indices differ
            if (!(T.at(d) == Tv.at(d))) maxnorm++;
            //maxnorm = maxnorm + max(T.at(d),Tv.at(d));
        }

        if (maxnorm > 1) {
            goto goon;
        }


        hier_new = 0.0;

        CalcHierbyPrew_Neumann(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), list);


        // hier = hier + hier_new | innerIterationssubgrid.getSubgridleq(Depth T);
        innerIterationsSubgrids.gotoEnd();
        do {
            Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();
            if (Tlocal <= T) {
                hier = hier + hier_new | *innerIterationsSubgrids.getSubgrid();
            }
        } while (innerIterationsSubgrids.previous());

        goon:;

    } while (outerIterationsSubgrids.previous());


    u = hier;

    CalcNodalBasis_inhomog(u, Tiefe);

}


void CalcUbyPrew_inhomog(VectorSparseG &prew, VectorSparseG &u, Depth &Tiefe) {

    AdaptiveSparseGrid *grid = u.getSparseGrid();

    // grid.getList();
    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*grid);
    VectorSparseG hier_new(*grid);

    u = 0.0;

    Depth TV(0);
    TV.set(3, 0);
    TV.set(1, 1);

    outerIterationsSubgrids.gotoEnd();
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();

        if (T == Tiefe) {

            hier_new = 0.0;
            CalcHierbyPrew_Neumann(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), list);





            // hier = hier + hier_new | innerIterationssubgrid.getSubgridleq(Depth T);
            innerIterationsSubgrids.gotoEnd();
            do {
                Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();

                if (Tlocal <= T) {

                    hier = hier + hier_new | *innerIterationsSubgrids.getSubgrid();
                }
            } while (innerIterationsSubgrids.previous());

        }

    } while (outerIterationsSubgrids.previous());


    u = hier;
    /*CalcNodalBasis_Neumann(u, Tiefe);*/

    CalcNodalBasis_inhomog(u, Tiefe);


}