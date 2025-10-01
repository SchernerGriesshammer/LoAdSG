//
// Created by scherner on 22.04.21.
//
#include "BasisTransformations.h"
#include "../tests/testing.h"

#include "../iterator/RectangularIterator.h"
#include "../sgrid/komponente.h"

/**
 * berechnet hierarchische Basis kleiner gleich Tsubgrid
 **/
void CalcHierarchicalBasis(VectorSparseG &vectorSparseG, const Depth &Tsubgrid, ListOfDepthOrderedSubgrids &orderedSubgrids) {


    //ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());

    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();

            if (T.at(d) != 0 && T <= Tsubgrid) {
                //SubtractHierUberschuss(u, T, d);
                vectorSparseG = vectorSparseG - 0.5 * (Links(vectorSparseG) + Rechts(vectorSparseG)) | *iterSubgrids.getSubgrid();
            }

        } while (iterSubgrids.previous());
    }


}

/**
 * berechnet hierarchische Basis kleiner gleich Tsubgrid
 **/
void CalcHierarchicalBasis_NEU(VectorSparseG &vectorSparseG, const Depth &Tsubgrid, std::list<Depth> SortierteTiefen) {





    for (int d = 0; d < DimensionSparseGrid; ++d) {


        for(Depth T:SortierteTiefen){




            if (T.at(d) != 0 && T <= Tsubgrid) {
                SingleDepthHashGrid &depthGrid = vectorSparseG.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                const auto &mapping = depthGrid._mapPosToGridPos;

                for (size_t i = 0; i < mapping.size(); i++) {
                    if (vectorSparseG.getSparseGrid()->getActiveTable()[mapping[i]]) {
                        IndexDimension I = vectorSparseG.getSparseGrid()->getIndexOfTable(mapping[i]);
                        IndexDimension Ileft = I.nextLeft(d);
                        IndexDimension Iright = I.nextRight(d);
                        double value = vectorSparseG.getValue(I) -
                                       0.5 * (vectorSparseG.getValue(Ileft) + vectorSparseG.getValue(Iright));
                        vectorSparseG.setValue(mapping[i], value);
                    }
                }


            }

        }
    }


}

void CalcHierarchicalBasisExact(VectorSparseG &u, const Depth &Tsubgrid) {

    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) != 0 && T == Tsubgrid) {
                u = u - 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //SubtractHierUberschuss(u, T, d);
            }
        } while (iterSubgrids.previous());
    }
}

void GetPrewaveletCoefficients(VectorSparseG *prew, VectorSparseG *u, ZusammenhangsKomponente &komponente) {

    Depth T = komponente.getDepth();
    AdaptiveSparseGrid_Base *grid = u->getSparseGrid();

    for (int d = 0; d < DimensionSparseGrid; d++) {

        IndexDimension minI = komponente.getMinIndex();
        IndexDimension maxI = komponente.getMaxIndex();


        IndexDimension iterIndex = minI;
        unsigned long k;
        int size = 1;



        // Berechne Größe der Matrix
        while (grid->occupied(k, iterIndex.nextRight(d, T.at(d))) &&
               !iterIndex.nextRight(d, T.at(d)).isAtRightBoundary(d) && grid->workonindex(k)) {

            size++;
            iterIndex = iterIndex.nextRight(d, T.at(d));
        }


        double p[size];
        double nodal[size];




        // M.create() in constructor!
        PrewaveletMatrixHomogen M(size, komponente.getDepth(), prew->getSparseGrid());

        M.create(minI, d);






        // "fixiere" Basiskoordinate Bsp: minI = (0,0,0), maxI = (1,1,1) --> maxI = (1,0,1)
        maxI.replace(d, IndexOneD(minI.getIndex(d)));
        iterIndex = minI;

        for (RectangularIterator iter(minI, maxI, komponente.getDepth()); iter.goon(); ++iter) {
            // Idee: OpenMP für jeden Index
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
                //u->setValue(iterIndex, p[i]);
                *u = p[i] | iterIndex;
                iterIndex = iterIndex.nextRight(d, T.at(d));
            }

        }
    }
}

void CalcNodalBasis(VectorSparseG &u) {

    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {

        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);


        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) != 0) {
                u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //AddHierUberschuss(u, T, d);
            }

        } while (iterSubgrids.next());
    }


}


void CalcNodalBasis(VectorSparseG &u, ListOfDepthOrderedSubgrids& orderedSubgrids) {


    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {

        //ShiftOperator Links(d, Left);
        //ShiftOperator Rechts(d, Right);


        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T.at(d) != 0) {
                SingleDepthHashGrid& depthGrid = u.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                const auto& mapping = depthGrid._mapPosToGridPos;
                if(depthGrid.getNumberOfEntries()>0) {
                    //#pragma omp parallel for
                    for (size_t i = 0; i < mapping.size(); i++) {
                        IndexDimension center = depthGrid._map.getIndexOfTable(i);
                        IndexDimension leftI = center.nextLeft(d);
                        IndexDimension rightI = center.nextRight(d);

                        double value = 0.5 * (u.getValue(leftI) + u.getValue(rightI));
                        u.addToValue(mapping[i], value);
                    }
                }



                 //u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //AddHierUberschuss(u, T, d);
            }

        } while (iterSubgrids.next());
    }


}

void CalcNodalBasis_NEU(VectorSparseG &u,std::list<Depth>& SortierteTiefen) {




    for (int d = 0; d < DimensionSparseGrid; ++d) {



        for (auto it = SortierteTiefen.rbegin(); it != SortierteTiefen.rend(); ++it){
            Depth T = *it;
            if (T.at(d) != 0) {
                SingleDepthHashGrid& depthGrid = u.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                const auto& mapping = depthGrid._mapPosToGridPos;
                if(depthGrid.getNumberOfEntries()>0) {
#pragma omp parallel for
                    for (size_t i = 0; i < mapping.size(); i++) {
                        IndexDimension center = depthGrid._map.getIndexOfTable(i);
                        IndexDimension leftI = center.nextLeft(d);
                        IndexDimension rightI = center.nextRight(d);

                        double value = 0.5 * (u.getValue(leftI) + u.getValue(rightI));
                        u.addToValue(mapping[i], value);
                    }
                }


            }

        }
    }


}


void CalcNodalBasis(VectorSparseG &u, const Depth &Tiefe, ListOfDepthOrderedSubgrids &orderedSubgrids) {

    // ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoBegin();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();
            if (T <= Tiefe) {
                if (T.at(d) != 0) {
                    u = u + 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                    //AddHierUberschuss(u, T, d);
                }
            }
        } while (iterSubgrids.next());
    }
}

void CalcNodalBasisExactDepth(VectorSparseG &hier, VectorSparseG &u, const Depth &Tiefe) {

    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);


    VectorSparseG hier_exact(*u.getSparseGrid());
    iterSubgrids.gotoBegin();
    do {
        Depth Tlocal = iterSubgrids.getSubgrid()->getT();
        if (Tlocal == Tiefe)
            hier_exact = hier | *iterSubgrids.getSubgrid();

    } while (iterSubgrids.next());

    u = hier_exact;

    CalcNodalBasis(u, Tiefe, orderedSubgrids);


}

void CalcHierbyPrew(VectorSparseG &prew, VectorSparseG &hier, SubgridFixedDepth &subgrid,
                    ListOfDepthOrderedSubgrids &orderedSubgrids) {


    SubgridFixedDepth::iterator iter(subgrid);

    Depth T = subgrid.getT();

    do {

        IndexDimension Index = iter.getPoint();
        double coeff = prew.getValue(Index);


        for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
            double basis_coeff = 1.0;
            IndexDimension J = Index.nextFive(&mc, Depth(Index), &basis_coeff);

            double val = hier.getValue(J) + coeff * basis_coeff;


            // hier = value | J;
            //hier.setValue(J, val);
            hier = val | J;
        }


    } while (iter.next());


    CalcHierarchicalBasis(hier, T, orderedSubgrids);

}

void CalcHierbyPrew_NEU(VectorSparseG &prew, VectorSparseG &hier, Depth& T, std::list<Depth> SortierteTiefen) {






    SingleDepthHashGrid &depthGrid = prew.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
    const auto &mapping = depthGrid._mapPosToGridPos;

    for (size_t i = 0; i < mapping.size(); i++) {
        if(prew.getSparseGrid()->getActiveTable()[mapping[i]]) {
            IndexDimension Index = prew.getSparseGrid()->getIndexOfTable(mapping[i]);

            double coeff = prew.getValue(mapping[i]);


            for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                double basis_coeff = 1.0;
                IndexDimension J = Index.nextFive(&mc, Depth(Index), &basis_coeff);

                double val = hier.getValue(J) + coeff * basis_coeff;
                hier.setValue(J, val);
            }
        }


    }


    CalcHierarchicalBasis_NEU(hier, T,SortierteTiefen);

}


void CalcHierarchicalBasis(VectorSparseG &u) {

    ListOfDepthOrderedSubgrids orderedSubgrids(*u.getSparseGrid());
    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();

            if (T.at(d) != 0) {
                u = u - 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //SubtractHierUberschuss(u, T, d);
            }
        } while (iterSubgrids.previous());
    }
}

void CalcHierarchicalBasis(VectorSparseG &u, ListOfDepthOrderedSubgrids &orderedSubgrids) {


    ListOfDepthOrderedSubgrids::iterator iterSubgrids(orderedSubgrids);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        ShiftOperator Links(d, Left);
        ShiftOperator Rechts(d, Right);

        iterSubgrids.gotoEnd();
        do {
            Depth T = iterSubgrids.getSubgrid()->getT();

            if (T.at(d) != 0) {
                u = u - 0.5 * (Links(u) + Rechts(u)) | *iterSubgrids.getSubgrid();
                //SubtractHierUberschuss(u, T, d);
            }
        } while (iterSubgrids.previous());
    }
}


void CalcHierarchicalBasis_NEU(VectorSparseG &u,std::list<Depth>& SortierteTiefen) {




    for (int d = 0; d < DimensionSparseGrid; ++d) {

        for (Depth T: SortierteTiefen) {


            if (T.at(d) != 0) {

                SingleDepthHashGrid &depthGrid = u.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                auto mapping = depthGrid._mapPosToGridPos;
                if (depthGrid.getNumberOfEntries() > 0) {
#pragma omp parallel for
                    for (size_t i = 0; i < mapping.size(); i++) {
                        IndexDimension center = depthGrid._map.getIndexOfTable(i);

                            IndexDimension leftI = center.nextLeft(d);
                            IndexDimension rightI = center.nextRight(d);

                            double value = u.getValue(mapping[i]) - 0.5 * (u.getValue(leftI) + u.getValue(rightI));

                            u.setValue(mapping[i], value);
                    }

                }
            }
        }

    }
}


void SolvePrewavelet(VectorSparseG &prew, VectorSparseG &u, const Depth &T) {

    ZusammenhangsKomponente komponente(u.getSparseGrid(), T);
    unsigned long maxoccup = u.getSparseGrid()->getMaximalOccupiedSecondTable();
    unsigned long iNow_mapping=-1;
    for (unsigned long iNow = komponente.findAnicht(-1,iNow_mapping); iNow < maxoccup; iNow = komponente.findAnicht(iNow,iNow_mapping)) {
        // sets minI = maxI
        if (prew.getSparseGrid()->workonindex(iNow)) {
            IndexDimension IndexNow = komponente.StartSearchComponent(iNow);
            // marks points and sets new minI and maxI
            komponente.recursiveMarkArbeite(IndexNow, iNow);

            GetPrewaveletCoefficients(&prew, &u, komponente);
        }

    }
}

void calcPrewByNodal(VectorSparseG &prew, VectorSparseG &u) {
    AdaptiveSparseGrid *grid = u.getSparseGrid();


    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*grid);
    VectorSparseG u_copy(*grid);
    u_copy = u;
    // hier_new braucht man nicht unbedingt
    VectorSparseG hier_new(*grid);
    prew = 0.0;
    hier = u;
    CalcHierarchicalBasis(hier);


    outerIterationsSubgrids.gotoEnd();
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();



        // hier: Iterator über Zusammenhangskomp. direkt in subgrid?
        // anderer Name?


        Depth zero_depth(0);

        if (T > zero_depth) {
            // berechnet prewavelet koeff der tiefe T
            // berechnet u aus hier auf level

            SolvePrewavelet(prew, u, T);


            hier_new = 0.0;
            CalcHierbyPrew(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), orderedSubgrids);

            //ziehe hierarchischen Überschuss ab

            innerIterationsSubgrids.gotoEnd();
            do {
                Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();

                if (Tlocal <= T) {
                    hier = hier - hier_new | *innerIterationsSubgrids.getSubgrid();
                }
            } while (innerIterationsSubgrids.previous());

            u = hier;
            CalcNodalBasis(u);


        }
    } while (outerIterationsSubgrids.previous());
    u = u_copy;

}

void calcPrewByNodal_NEU(VectorSparseG &prew, VectorSparseG &u, std::list<Depth>& SortierteTiefen) {

    AdaptiveSparseGrid *grid = u.getSparseGrid();


    VectorSparseG hier(*grid);
    VectorSparseG u_copy(*grid);
    u_copy = u;
    // hier_new braucht man nicht unbedingt
    VectorSparseG hier_new(*grid);
    prew = 0.0;
    hier = u;
    CalcHierarchicalBasis_NEU(hier,SortierteTiefen);


    for(Depth T:SortierteTiefen){

        Depth zero_depth(0);

        if (T > zero_depth) {

            // berechnet prewavelet koeff der tiefe T
            // berechnet u aus hier auf level

            SolvePrewavelet(prew, u, T);


            hier_new = 0.0;

            CalcHierbyPrew_NEU(prew, hier_new,T,SortierteTiefen);

            for(Depth Tlocal:SortierteTiefen){


                if (Tlocal <= T) {
                    SingleDepthHashGrid &depthGrid = prew.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
                    const auto &mapping = depthGrid._mapPosToGridPos;
#pragma omp parallel for
                    for (size_t i = 0; i < mapping.size(); i++) {
                        if (prew.getSparseGrid()->getActiveTable()[mapping[i]]) {
                            double val = hier.getValue(mapping[i]) - hier_new.getValue(mapping[i]);
                            hier.setValue(mapping[i], val);
                        }
                    }
                }
            }

            u = hier;
            CalcNodalBasis_NEU(u,SortierteTiefen);



        }
    }


    u = u_copy;
}


void calcPrewByNodal(VectorSparseG &prew, VectorSparseG &u, ListOfDepthOrderedSubgrids& orderedSubgrids) {


    AdaptiveSparseGrid *grid = u.getSparseGrid();




    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*grid);
    VectorSparseG u_copy(*grid);
    u_copy = u;
    // hier_new braucht man nicht unbedingt
    VectorSparseG hier_new(*grid);
    prew = 0.0;
    hier = u;
    CalcHierarchicalBasis(hier,orderedSubgrids);
    int count =0;
    outerIterationsSubgrids.gotoEnd();
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();


        // hier: Iterator über Zusammenhangskomp. direkt in subgrid?
        // anderer Name?


        Depth zero_depth(0);

        if (T > zero_depth) {
            // berechnet prewavelet koeff der tiefe T
            // berechnet u aus hier auf level

            SolvePrewavelet(prew, u, T);


            hier_new = 0.0;
            CalcHierbyPrew(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), orderedSubgrids);

            //ziehe hierarchischen Überschuss ab

            innerIterationsSubgrids.gotoEnd();
            do {
                Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();

                if (Tlocal <= T) {
                    hier = hier - hier_new | *innerIterationsSubgrids.getSubgrid();
                }
            } while (innerIterationsSubgrids.previous());

            u = hier;
            CalcNodalBasis(u,orderedSubgrids);
            count++;

        }
    } while (outerIterationsSubgrids.previous());
    u = u_copy;


}

void calcNodalByPrew(VectorSparseG &prew, VectorSparseG &u) {
    AdaptiveSparseGrid *grid = u.getSparseGrid();

    // grid.getList();
    ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    VectorSparseG hier(*prew.getSparseGrid());
    VectorSparseG hier_new(*grid);

    u = 0.0;

    outerIterationsSubgrids.gotoEnd();
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();

        hier_new = 0.0;
        CalcHierbyPrew(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), orderedSubgrids);


        // hier = hier + hier_new | innerIterationssubgrid.getSubgridleq(Depth T);
        innerIterationsSubgrids.gotoEnd();
        do {
            Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();
            if (Tlocal <= T) {
                hier = hier + hier_new | *innerIterationsSubgrids.getSubgrid();
            }
        } while (innerIterationsSubgrids.previous());


    } while (outerIterationsSubgrids.previous());

    u = hier;
    CalcNodalBasis(u);

}


void CalcUbyPrewRestrictions(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions,
                             ListOfDepthOrderedSubgrids &orderedSubgrids) {

    AdaptiveSparseGrid *grid = u.getSparseGrid();

    // grid.getList();
    //ListOfDepthOrderedSubgrids orderedSubgrids(*grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    ListOfDepthOrderedSubgrids::iterator innerIterationsSubgrids(orderedSubgrids);

    //VectorSparseG hier(*grid);
    VectorSparseG hier_new(*grid);

    u = 0.0;

    Depth Tv(Tiefe);
    for (int d = 0; d < DimensionSparseGrid; d++) {
        if (restrictions[d])Tv.set(Tiefe.at(d) - 1, d);
    }


    outerIterationsSubgrids.gotoEnd();
    unsigned int lonenorm = outerIterationsSubgrids.getDepth().LoneNorm();

    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();

        int maxnorm = 0;
        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (!restrictions[d] && !(T.at(d) <= Tiefe.at(d)))goto goon;
            if (restrictions[d] && T.at(d) != Tiefe.at(d))goto goon;

            if (!(T.at(d) == Tv.at(d)))maxnorm = maxnorm + max(T.at(d), Tv.at(d));
        }

        if (maxnorm >= lonenorm)goto goon;
        hier_new = 0.0;
        CalcHierbyPrew(prew, hier_new, *outerIterationsSubgrids.getSubgrid(), orderedSubgrids);




        // hier = hier + hier_new | innerIterationssubgrid.getSubgridleq(Depth T);
        innerIterationsSubgrids.gotoEnd();
        do {
            Depth Tlocal = innerIterationsSubgrids.getSubgrid()->getT();

            if (Tlocal <= T) {

                u = u + hier_new | *innerIterationsSubgrids.getSubgrid();
            }
        } while (innerIterationsSubgrids.previous());

        goon:;

    } while (outerIterationsSubgrids.previous());


    //u = hier;
    CalcNodalBasis(u, Tiefe, orderedSubgrids);


}

