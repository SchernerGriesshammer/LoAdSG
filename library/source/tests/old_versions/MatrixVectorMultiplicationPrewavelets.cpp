//
// Created by scherner on 22.04.21.
//

#include "MatrixVectorMultiplicationPrewavelets.h"
#include "../../mympi.h"
#include "../testing.h"
#include "../../stencils/Stencil.h"


bool *CasesIterator::getcase(int snum_) {
    static bool cases[DimensionSparseGrid];


    for (int i = 0; i < DimensionSparseGrid; ++i) cases[i] = false;

    for (int i = 0; i < DimensionSparseGrid; ++i) {
        // the following condition checks
        // if the ith bit of snum in binary form
        // is 1 or not
        if ((snum_ & (1 << i)) != 0) {

            cases[i] = true;
        }
    }
    return cases;
};
///////////////////////////////////////////////////////////////////// Matrix-Vector Dirichlet

void MatrixVectorPrewavelet(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, MultiLevelVector &gM) {

    Ax = 0.0;
    AdaptiveSparseGrid *grid = prew.getSparseGrid();
    VectorSparseG u(*grid);
    VectorSparseG Ax_neu(*grid);
    VectorSparseG Ax_hier(*grid);
    VectorSparseG g(*grid);
    VectorSparseG z(*grid);

    MultiLevelVector nodal(*gM.getSparseGrid());



    for (CasesIterator iter; iter.goon(); ++iter) {
        Ax_neu = 0.0;
        bool *restrictions = iter.getcase();
        ListOfDepthOrderedSubgrids list(*grid, restrictions);


        //setze alle punkte von Ax_neu gleich 0
        grid->WorkOnHangingNodes = true;
        Ax_neu = 0.0;
        //Ax_hier = 0.0;

        g = 0.0;
        z = 0.0;
        gM = 0.0;

        //PP_Alg(prew, u, Ax_neu, type, restrictions, Ax_hier,g,z);
        MatrixVector_Case(prew, u, Ax_neu, type, restrictions, g, z, gM, list, nodal);


        grid->WorkOnHangingNodes = false;

        Ax = Ax + Ax_neu;


    }
}

// Alte Version MatrixVector_Case
/*
void PP_Alg(VectorSparseG &prew, VectorSparseG &u, VectorSparseG &Ax, StencilType type, bool *restrictions,VectorSparseG& Ax_hier, VectorSparseG& g, VectorSparseG& z) {


    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();

    // Sortiere Tiefen nach Restriktionen bzw. Prolongationen
    // Iteriere rekursiv über d = 1,...,DimensionSparseGrid
    // Restriktion in Richtung d --> fine to coarse (Ende hier bereits bei t = 2, da nicht von t = 1 auf t = 0 restringiert wird)
    // Prolongation in Richtung d --> coarse to fine
    // Beispiel: restrictions = (0,1)
    // maximale Tiefe = (3,3)
    // (1,3) -> (1,2) ->
    // (2,3) -> (2,2) ->
    // (3,3) -> (3,2)
    ListOfDepthOrderedSubgrids list(*prew.getSparseGrid(), restrictions);

    Depth Told = *list.getSortierteTiefen()->begin();



    for (Depth T:*list.getSortierteTiefen()) {


        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (!restrictions[d] && (Told.at(d) != T.at(d))) {
                // alle Tiefen die vorher bereits verwendet wurden, werden nicht mehr gebraucht
                // weil sie in CalcUbyPrew neu verwendet werden.
                // Das geht sicher klüger, dann braucht man aber auch hier MultilevelVectoren und Prolongationsoperatoren

                grid->WorkOnHangingNodes = true;
                g = 0.0;
                z = 0.0;

            }
        }



        prew.getSparseGrid()->WorkOnHangingNodes = true;
        u = 0.0;
        // Berechne auf Tiefe T alle Tiefen für die gilt:
        // if restrictions[d] == 0: Tiefe >= T
        //    restrictions[d] == 1: Tiefe < T

        CalcUbyPrewRestrictions(prew, u, T, restrictions);



        //ApplyStencil on HangingNodes
        // Wende Stencil auf Knoten der Tiefe T oder kleiner an
        z = 0.0;

        ApplyStencil(u, z, T, type);


        // Addiere "
        z = z + g;
        g = z;




        //restrict on hangingnodes
        Depth Tcoarse(T);
        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {

                Tcoarse.set(T.at(d) - 1, d);
                z = g;
                g = 0.0;
                restriction(z, g, T.at(d) - 1, d);

            }
        }






        //convert on activenodes
        prew.getSparseGrid()->WorkOnHangingNodes = false;
        ConvertToPrewavelet(g, Ax, Tcoarse,*list.getSubgrid(Tcoarse));

        Told = T;
    }

}

*/

void MatrixVector_Case(VectorSparseG &prew, VectorSparseG &u, VectorSparseG &Ax, StencilType type, bool *restrictions,
                       VectorSparseG &g, VectorSparseG &z, MultiLevelVector &gM, ListOfDepthOrderedSubgrids &list,
                       MultiLevelVector &nodal) {


    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();
    ListOfDepthOrderedSubgrids list2(*grid);


    // Sortiere Tiefen nach Restriktionen bzw. Prolongationen
    // Iteriere rekursiv über d = 1,...,DimensionSparseGrid
    // Restriktion in Richtung d --> fine to coarse (Ende hier bereits bei t = 2, da nicht von t = 1 auf t = 0 restringiert wird)
    // Prolongation in Richtung d --> coarse to fine
    // Beispiel: restrictions = (0,1)
    // maximale Tiefe = (3,3)
    // (1,3) -> (1,2) ->
    // (2,3) -> (2,2) ->
    // (3,3) -> (3,2)

    //ListOfDepthOrderedSubgrids list(*prew.getSparseGrid(), restrictions);
    list2.sortDepths(restrictions);
    for (Depth T: *list2.getSortierteTiefen()) {



        prew.getSparseGrid()->WorkOnHangingNodes = true;
        u = 0.0;

        CalcUbyPrewRestrictions(prew, u, T, restrictions, list2);



        //ApplyStencil on HangingNodes
        z = 0.0;
        ApplyStencil(u, z, T, type);



        // g = gM | T --> g.setMultiLevelValues(gM,T);
        g = 0.0;
        g.setMultiLevelValues(gM, T);


        z = z + g;



        // g = z nur falls keine Restriktionen durchgeführt werden
        g = z;



        //restrict on hangingnodes

        Depth Tcoarse(T);

        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {

                Tcoarse.set(T.at(d) - 1, d);

                g = 0.0;
                restriction(z, g, Tcoarse.at(d), d);
                z = 0.0;
                z.setMultiLevelValues(gM, Tcoarse);
                z = z + g;
                gM.setMultiLevelValues(g, Tcoarse);
            }
        }


        //convert on activenodes

        prew.getSparseGrid()->WorkOnHangingNodes = false;
          ConvertToPrewavelet(g, Ax, Tcoarse, *list2.getSubgrid(Tcoarse));



    }
}

void ApplyStencil(VectorSparseG &x, VectorSparseG &Ax, Depth &T, StencilType type) {
    AdaptiveSparseGrid_Base *grid = x.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    Depth Tzero(0);


    for (unsigned long i = 0; i <= maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);

        Depth Tlocal(Index);
        if (Tlocal > Tzero && Tlocal <= T) {
            double val = CalcStencilValue(Index, T, x, type);
            //Ax = val | Index;
            Ax.setValue(i,val);
        }
    }
}

double CalcStencilValue(IndexDimension Index, Depth &T, VectorSparseG &u, StencilType type) {
    PoissonStencil stencil;
    stencil.initialize(T);
    double val = 0.0;
    for (MultiDimCompass mc; mc.goon(); ++mc) {

        IndexDimension NextIndex;
        double val_new = 0.0;
        double val_test = 0.0;


        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (type == Mass) d = DimensionSparseGrid - 1;

            NextIndex = GetStencilValueIndex(Index, mc, &val_test, T, d, type);

            val_new = val_new + val_test;
        }

//       NextIndex = Index.nextThree(&mc,T.returnTiefen());
        if(mc.goon()) {
            //val_new = stencil.returnValue(Index, mc);
            val = val + val_new * u.getValue(NextIndex);
        }
    }
    return val;
}

IndexDimension GetStencilValueIndex(IndexDimension Index, MultiDimCompass mc, double *stencilvalue, Depth T, int dir,
                                    StencilType type) {
    double val = 1.0;

    IndexDimension J = Index;

    for (int d = 0; d < DimensionSparseGrid; ++d) {

        int exp = -1 * T.at(d);
        //if(T.at(d)==0) exp = -1;
        double h = pow(2, exp);


        if (d == dir && type == Stiffness) val = 1 / h * val;
        if (d != dir || type == Mass) val = (1.0 / 6.0) * h * val;

        int t = T.at(d);
        //if(T.at(d)==0) t = 1;
        Richtung r = mc.getRichtung(d);
        if (r == Links) {
            if (d == dir && type == Stiffness) val = -1.0 * val;
            J = J.nextLeft(d, t);

        }


        if (r == Rechts) {

            if (d == dir && type == Stiffness) val = -1 * val;
            J = J.nextRight(d, t);

        }

        if (r == Mitte) {
            if (d == dir && type == Stiffness)val = 2 * val;
            else if (d != dir || type == Mass) val = 4 * val;
        }


    }


    *stencilvalue = val;
    return J;
}

void ConvertToPrewavelet(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T, SubgridFixedDepth &subgrid) {

    unsigned long k;

    SubgridFixedDepth::iterator iter(subgrid);
    do {
        IndexDimension Index = iter.getPoint();
        unsigned long i = iter.geti();
        if(Ax_hier.getSparseGrid()->getActiveTable()[i]) {
            double coeff = 0.0;
            for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                double basis_coeff = 0.0;
                IndexDimension J = Index.nextFive(&mc, T, &basis_coeff);
                if(mc.goon()&&Ax_hier.getSparseGrid()->occupied(k,J)){
                    double val = Ax_hier.getValue(J);
                    double d = val*basis_coeff;
                    coeff = coeff + d;

                }
            }


            Ax.setValue(i, coeff);
        }
    } while (iter.next());

}

void ConvertToPrewavelet2(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T) {
    SingleDepthHashGrid& depthGrid = Ax_hier.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
    const auto& mapping = depthGrid._mapPosToGridPos;
    // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
    for (size_t i = 0; i < mapping.size(); i++)
    {
        if(Ax_hier.getSparseGrid()->getActiveTable()[mapping[i]]) {
            IndexDimension I = depthGrid._map.getIndexOfTable(i);
            double coeff = 0.0;
            IndexDimension J;

            for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                double basis_coeff = 1.0;
                J = I.nextFiveP(&mc, T, &basis_coeff);
                if (mc.goon())
                    coeff = coeff + Ax_hier.getValue(J) * basis_coeff;
            }

            Ax.setValue(mapping[i], coeff);
        }
    }
}


void ConvertDualNodalToPrewavelet(VectorSparseG &nodal_neumann, VectorSparseG &prewavelet_dirichlet, Depth &T) {
    SingleDepthHashGrid& depthGrid = nodal_neumann.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
    const auto& mapping = depthGrid._mapPosToGridPos;
    // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
    for (size_t i = 0; i < mapping.size(); i++)
    {
        if(nodal_neumann.getSparseGrid()->getActiveTable()[mapping[i]]) {
            IndexDimension I = depthGrid._map.getIndexOfTable(i);
            unsigned long k;
            if(prewavelet_dirichlet.getSparseGrid()->occupied(k,I))
                if(prewavelet_dirichlet.getSparseGrid()->getActiveTable()[k]){
                    double coeff = 0.0;
                    IndexDimension J;

                    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                        double basis_coeff = 1.0;
                        J = I.nextFiveP(&mc, T, &basis_coeff);
                        if (mc.goon())
                            coeff = coeff + nodal_neumann.getValue(J) * basis_coeff;
                    }

                    prewavelet_dirichlet.setValue(I, coeff);
                }
        }
    }
}


void restriction(VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {
    AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    Depth Tzero(0);
    for (unsigned long i = 0; i < maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);
        Depth Tlocal(Index);

        if (Index.getDepth(d) <= t && Tlocal > Tzero) {
            double value = 0.0;

            value = fine.getValue(i) + 0.5 * fine.getValue(Index.nextLeft(d, t + 1)) +
                    0.5 * fine.getValue(Index.nextRight(d, t + 1));



            coarse.setValue(i,value);
        }


    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////// Matrix-Vector Neumann
void MatrixVectorPrewavelet_inhomog(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, MultiLevelVector &gM) {


    Ax = 0.0;
    AdaptiveSparseGrid *grid = prew.getSparseGrid();

    VectorSparseG u(*grid);
    VectorSparseG Ax_neu(*grid);
    Process *mpi = grid->mpi;
    int numberofnodes = 1;
    // MPI_Comm_size(MPI_COMM_WORLD, &numberofnodes);

    int k = 1;


    for (CasesIterator iter; iter.goon(); ++iter) {
        Ax_neu = 0.0;
        if (mpi->getMyRank() == k % numberofnodes || true) {
            //cout << "do case "  << k << " on rank " << k % numberofnodes << endl;

            bool *restrictions = iter.getcase();

            grid->WorkOnHangingNodes = true;

            u = 0.0;

            //PP_Alg_inhomog(prew, u, Ax_neu, type, restrictions);
            gM = 0.0;

            test_neumann(prew, Ax_neu, type, restrictions, gM);
            grid->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;
        }

        k++;
    }
    //Ax_neu = Ax;
    //Ax.ReduceSum(0);
    //Ax_neu.Broadcast(0);

/*
    Ax_neu = 0.0;
    ConvertDualSpace(Ax,Ax_neu);
    Ax = 0.0;


    grid->WorkOnHangingNodes = false;
    Ax= Ax_neu;
*/

}

void test_neumann(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, bool *restrictions, MultiLevelVector &gM) {


    AdaptiveSparseGrid *grid = prew.getSparseGrid();

    ListOfDepthOrderedSubgrids list(*prew.getSparseGrid());
    list.SortiereTiefenBoundary(restrictions);


    VectorSparseG g(*prew.getSparseGrid());
    VectorSparseG z(*prew.getSparseGrid());

    VectorSparseG u(*prew.getSparseGrid());


    for (Depth T: *list.getSortierteTiefen()) {

        prew.getSparseGrid()->WorkOnHangingNodes = true;

        u = 0.0;
        CalcUbyPrewRestrictions_inhomog(prew, u, T, restrictions);




        //ApplyStencil on HangingNodes
        z = 0.0;
        ApplyStencil_inhomog(u, z, T, type);



        // g = gM | T --> g.setMultiLevelValues(gM,T);
        g = 0.0;
        g.setMultiLevelValues(gM, T);




        //  z = z+g;
        z = z + g;
        g = z;



        //restrict on hangingnodes

        Depth Tcoarse(T);
        VectorSparseG test(*grid);


        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {

                Tcoarse.set(T.at(d) - 1, d);

                g = 0.0;
                restriction_neumann(z, g, Tcoarse.at(d), d);
                z = 0.0;
                z.setMultiLevelValues(gM, Tcoarse);
                z = z + g;
                gM.setMultiLevelValues(g, Tcoarse);
            }
        }


        //convert on activenodes
        prew.getSparseGrid()->WorkOnHangingNodes = false;
        ConvertToNeumannPrewavelet(g, Ax, Tcoarse);


    }


}


double CalcStencilValue_Boundary(IndexDimension Index, Depth &T, VectorSparseG &u, StencilType type, Stencil stencil) {



    // Stencil stencil(type,T);

    double val = 0.0;
    IndexDimension NextIndex;
    for (MultiDimCompass mc; mc.goon(); ++mc) {

        double val_new = 0.0;
        double val_test = 0.0;


       for (int d = 0; d < DimensionSparseGrid; d++) {

            if (type == Mass) d = DimensionSparseGrid - 1;

            val_test = 1.0;
            NextIndex = GetStencilValueIndex_Boundary(Index, mc, &val_test, T, d, type);


            val_new = val_new + val_test;
        }

        val = val + val_new * u.getValue(NextIndex);


    }

    return val;
}
// möglichst inline
// per referenz


IndexDimension
GetStencilValueIndex_Boundary(IndexDimension Index, MultiDimCompass mc, double *stencilvalue, Depth T, int dir,
                              StencilType type) {
    double val = 1.0;

    IndexDimension J = Index;
    for (int d = 0; d < DimensionSparseGrid; ++d) {

        int exp = T.at(d);
        //if (T.at(d) == 0) exp = -1;

        // 2 ^ Tiefe -> nicht mit pow
        // h nur einmal berechnen und übergeben
        double h = pow(2, exp);
        h = 1/h;

        if (d == dir && type == Stiffness) val = 1 / h * val;
        if (d != dir || type == Mass) val = (1.0 / 6.0) * h * val;

        int t = T.at(d);
        //if (T.at(d) == 0) t = 1;
        Richtung r = mc.getRichtung(d);
        if (r == Links) {
            if (d == dir && type == Stiffness) {
                val = -1.0 * val;
            }
            if (J.isAtLeftBoundary(d)){
                val = 0.0;
                *stencilvalue = val;
                return J;
            }
            J = J.nextLeft(d, t);
        }
        if (r == Rechts) {
            if (d == dir && type == Stiffness) val = -1 * val;
            if (J.isAtRightBoundary(d)){
                val = 0.0;
                *stencilvalue = val;
                return J;
            }
            J = J.nextRight(d, t);
        }

        if (r == Mitte) {
            if (d == dir && type == Stiffness) {
                if ((!Index.isAtLeftBoundary(d)) && (!Index.isAtRightBoundary(d)))
                    val = 2.0 * val;
            } else if (d != dir || type == Mass) {
                if ((!Index.isAtLeftBoundary(d)) && (!Index.isAtRightBoundary(d)))
                    val = 4.0 * val;
                else {
                    val = 2.0 * val;
                }

            }

        }


    }

    *stencilvalue = val;
    return J;
}

double
GetStencilValue_Boundary(MultiDimCompass mc, Depth T, int dir,
                         StencilType type) {

    double val = 1.0;

    for (int d = 0; d < DimensionSparseGrid; ++d) {

        int exp = -1 * T.at(d);
        //if (T.at(d) == 0) exp = -1;
        double h = pow(2, exp);


        if (d == dir && type == Stiffness) val = 1 / h * val;
        if (d != dir || type == Mass) val = (1.0 / 6.0) * h * val;

        int t = T.at(d);
        //if (T.at(d) == 0) t = 1;
        Richtung r = mc.getRichtung(d);
        if (r == Links || r == Rechts) {
            if (d == dir && type == Stiffness) {
                val = -1.0 * val;
            }
        }


        if (r == Mitte) {
            if (d == dir && type == Stiffness) {

                val = 2.0 * val;
            } else if (d != dir || type == Mass) {
                val = 4.0 * val;

            }

        }

    }


    return val;


}


bool getNextIndex(IndexDimension Index, MultiDimCompass mc, Depth T, IndexDimension &NextIndex) {
    IndexDimension J = Index;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        Richtung r = mc.getRichtung(d);
        int t = T.at(d);
        if (r == Links) {
                if(Index.isAtLeftBoundary(d))
                    return false;

                J = J.nextLeft(d, t);
        }
        if (r == Rechts) {
                if(Index.isAtRightBoundary(d))
                    return false;

                J = J.nextRight(d, t);
        }

    }
    NextIndex = J;
    return true;
}


void ApplyStencil_inhomog(VectorSparseG &x, VectorSparseG &Ax, Depth &T, StencilType type) {

    AdaptiveSparseGrid_Base *grid = x.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    PoissonStencil poissonStencil;


    for (unsigned long i = 0; i < maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);

        Depth Tlocal(Index);
        if (Tlocal <= T) {

            double val = getLocalStencilBoundary<PoissonStencil>(Index, T, x,poissonStencil);
            //Ax = val | Index;
            Ax.setValue(i,val);
        }
    }
}

void restriction_neumann(VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {

    AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();

    for (unsigned long i = 0; i < maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);

        if (Index.getDepth(d) <= t) {

            double value = fine.getValue(Index);
            double valueL = 0.0;
            double valueR = 0.0;
            if (t >= 0) {
                if (!(Index.isAtRightBoundary(d))) valueL = 0.5 * fine.getValue(Index.nextRight(d, t + 1));
                if (!(Index.isAtLeftBoundary(d))) valueR = 0.5 * fine.getValue(Index.nextLeft(d, t + 1));
            }


            value = value + valueL + valueR;
            //coarse.setValue(Index, value);
            coarse = value | Index;
        }


    }


}

void ConvertToNeumannPrewavelet(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T) {


    AdaptiveSparseGrid_Base *grid = Ax.getSparseGrid();


    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    for (unsigned long i = 0; i < maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);

        Depth Tneu(Index);
        if (Tneu == T && grid->getActiveTable()[i]) {

                double coeff = 0.0;
                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    IndexDimension J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff);

                    if(mc.goon())
                    coeff = coeff + Ax_hier.getValue(J) * basis_coeff;

            }
            Ax.setValue(i,coeff);
        }
    }





}

void
ConvertToNeumannPrewavelet2(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T, ListOfDepthOrderedSubgrids &list) {

    AdaptiveSparseGrid *grid = Ax.getSparseGrid();
    VectorSparseG Ax_hierCP(*grid);
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();

    SubgridFixedDepth::iterator iter(*list.getSubgrid(T));
    iter.gotobegin();
    do {


        IndexDimension Index = iter.getPoint();

        Depth Tneu(Index);
        if (Tneu == T) {
            if (grid->checkMultiDimFiveStencil(Index)) {
                double coeff = 0.0;
                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    IndexDimension J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff);

                    coeff = coeff + Ax_hier.getValue(J) * basis_coeff;
                }


                //Ax.setValue(Index,coeff);
                Ax = coeff | Index;
            }
        }
    } while (iter.next());

}
/////////////////////////////////////////////////////



//////////////////////// Volle Gitter

void ApplyStencil(VectorSparseG &x, VectorSparseG &Ax, StencilType type) {
    AdaptiveSparseGrid_Base *grid = x.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();

    for (unsigned long i = 0; i <= maxocc; i++) {
        IndexDimension Index = grid->getIndexOfTable(i);
        Depth T(0);


        bool doit = true;
        for (int d = 0; d < DimensionSparseGrid; d++) {

            if (Index.getDepth(d) == 0) doit = false;
            else {
                int t = x.getSparseGrid()->getMaxDepth(d, Index);

                T.set(t, d);
            }

        }

        if (doit) {
            double val = CalcStencilValue(Index, T, x, type);
            //Ax.setValue(Index,val);
            Ax = val | Index;

        }
    }
}







void
prolongation3(VectorSparseG &fine, VectorSparseG &coarse, int t, int d, Depth T, ListOfDepthOrderedSubgrids &list) {
    AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    Depth Tzero(0);
    int tneu = t;

    IndexDimension minI;
    IndexDimension maxI;


    SubgridFixedDepth::iterator iter(*list.getSubgrid(T));


    iter.gotobegin();
    do {
        IndexDimension Index = iter.getPoint();


        if (Index.getDepth(d) == t) {
            double value = fine.getValue(Index);
            if (Index.isAtLeftBoundary(d))
                value = fine.getValue(Index);

            if (Index.isAtRightBoundary(d))
                value = fine.getValue(Index);

            if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                value = fine.getValue(Index) + 0.5 * fine.getValue(Index.nextRight(d))
                        + 0.5 * fine.getValue(Index.nextLeft(d));
            //coarse.setValue(Index, value);
            coarse = value | Index;


        }


    } while (iter.next());

}


void restriction_local(VectorSparseG &fine, VectorSparseG &coarse, Depth Tfine, Depth Tcoarse) {

    //coarse = fine;
    AdaptiveSparseGrid* test = fine.getSparseGrid();
    VectorSparseG fine_new(*test);
    VectorSparseG coarse_new(*test);
    fine_new = fine;

    for (int d = 0; d < DimensionSparseGrid; d++) {
        for (int t = Tfine.at(d); t >= int(Tcoarse.at(d)); --t) {

            coarse_new = 0.0;
            restriction(fine_new, coarse_new, t, d);
            fine_new = coarse_new;

        }
    }

    coarse = coarse_new;
}






void matrixmult_prew_fullgrid(VectorSparseG &x, VectorSparseG &Ax, StencilType type) {
    AdaptiveSparseGrid *grid = x.getSparseGrid();


    ListOfDepthOrderedSubgrids list(*grid);
    ListOfDepthOrderedSubgrids::iterator outer(list);
    outer.gotoEnd();
    Depth T = outer.getSubgrid()->getT();


    VectorSparseG nodal(*grid);
    VectorSparseG Ax_test(*grid);


    VectorSparseG finest_nodal(*grid);
    VectorSparseG Ax_hier(*grid);

    calcNodalByPrew(x, nodal);


    ApplyStencil(nodal, finest_nodal, type);



    Depth Tzero(0);
    do {

        Depth Tcoarse = outer.getSubgrid()->getT();
        if (Tcoarse > Tzero) {


            restriction_local(finest_nodal, Ax_hier, T, Tcoarse);


            Ax_test = 0.0;
            ConvertToPrewavelet(Ax_hier, Ax_test, Tcoarse, *list.getSubgrid(Tcoarse));


            Ax = Ax + Ax_test;
            Ax_hier = 0.0;
        }
    } while (outer.previous());


}




Stencil::Stencil(Depth T_, StencilType type_) : T(T_), type(type_) {

    for (MultiDimCompass mc; mc.goon(); ++mc) {

        double val_new = 0.0;
        double val_test = 0.0;


        for (int d = 0; d < DimensionSparseGrid; d++) {

            if (type == Mass) d = DimensionSparseGrid - 1;

            val_test = GetStencilValue_Boundary(mc, T, d, type);

            val_new = val_new + val_test;


        }


        values[mc.getShiftNumber()] = val_new;


    }
}

double Stencil::returnValue(MultiDimCompass &mc, IndexDimension Index) {


    double value = 0.0;
    if (Index.isNotAtBoundary()) {

        //value = values[mc.getShiftNumber()];
        double val_test = 0.0;
        for (int d = 0; d < DimensionSparseGrid; d++) {

            if (type == Mass) d = DimensionSparseGrid-1;
            val_test = 0.0;
            GetStencilValueIndex_Boundary(Index, mc, &val_test, T, d, type);


            value = value + val_test;
        }


    } else {
        value = 0.0;
        double val_test = 0.0;
        for (int d = 0; d < DimensionSparseGrid; d++) {

            if (type == Mass) d = DimensionSparseGrid - 1;

            GetStencilValueIndex_Boundary(Index, mc, &val_test, T, d, type);


            value = value + val_test;
        }


    }
    return value;


}
