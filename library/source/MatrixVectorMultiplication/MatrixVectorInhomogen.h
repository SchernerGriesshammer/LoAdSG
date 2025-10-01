//
// Created by to35jepo on 12/7/22.
//

#ifndef SGRUN_MATRIXVECTORINHOMOGEN_H_H
#define SGRUN_MATRIXVECTORINHOMOGEN_H_H
#include "../extemp/vector.h"
#include "../tests/old_versions/MatrixVectorMultiplicationPrewavelets.h"
#include "../sgrid/multiDepthHashGrid.h"
#include "../BasisTransformations/BasisTransformations.h"
#include "../stencils/PoissonStencil.h"
#include "../iterator/depthIterator.h"


class MatrixVectorInhomogen { ;
public:
    MatrixVectorInhomogen(AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid) : gM(mgrid), g(grid), z(grid), u(grid),
                                                                                           nodal(mgrid),
                                                                                           u_new(grid),
                                                                                           u_old(grid), Ax_neu(grid), list(grid), P(mgrid){




    };


    template<class Problem>
    void multiplication(VectorSparseG &prew, VectorSparseG &Ax, Problem &matrixProblem);


    /**
     *
     * @tparam Problem
     * @param prew
     * @param Ax
     * @param matrixProblem
     *
     * Berechnet rechte Seite der schwachen Formulierung und transformiert den Dualraum Neumann --> Dirichlet
     */
    template<class Problem>
    void multiplicationRHSHomogen(VectorSparseG &prew, VectorSparseG &Ax, Problem &matrixProblem);




    // hier wird if(type==irgendwas)
    //void multiplication(VectorSparseG &prew, VectorSparseG &Ax, StencilType type);




private:

    MultiLevelVector gM;
    MultiLevelVector nodal;
    MultiLevelVector P;
    VectorSparseG Ax_neu;
    VectorSparseG g;
    VectorSparseG z;
    VectorSparseG u;
    VectorSparseG u_new;
    VectorSparseG u_old;

    ListOfDepthOrderedSubgrids list;

    template<class Problem>
    void CaseFunction(VectorSparseG &prew, bool *restrictions, Problem &matrixProblem);

    template<class Problem>
    void CaseFunctionTransformation(VectorSparseG &prew, bool *restrictions, Problem &matrixProblem);

    inline void calcU(VectorSparseG &prew, Depth &Tiefe, bool *restrictions);

    inline void calcNodal(VectorSparseG &prew, Depth &Tiefe);

    template<class Problem>
    void applyProblemGlobal(VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {
        AdaptiveSparseGrid_Base *grid = input.getSparseGrid();
        int maxocc = grid->getMaximalOccupiedSecondTable();
        for (unsigned long i = 0; i < maxocc; i++) {
            IndexDimension Index = grid->getIndexOfTable(i);
            Depth Tlocal(Index);
            if (Tlocal <= T) {

                double val = applyProblemLocal(Index, input, matrixProblem, T);
                //output = val | Index;
                output.setValue(i, val);
            }
        }
    };

    template<class Problem>
    double applyProblemLocal(IndexDimension &Index, VectorSparseG &input, Problem &matrixProblem, Depth &T) {

        double val = 0.0;
        IndexDimension NextTest;
        bool exist_index;
        for (MultiDimCompass mc; mc.goon(); ++mc) {
            NextTest = Index.nextThree_boundary(&mc, T.returnTiefen());
            //exist_index = getNextIndex(Index, mc, T, NextTest);
            if(mc.goon()) {
                double value = 0.0;
                value = matrixProblem.returnValue(Index, mc);
                val = val + value * input.getValue(NextTest);
            }
        }
        return val;
    };

    inline void prolongation1D(VectorSparseG &coarse, VectorSparseG &fine, int &t, int &d) {
        AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
        unsigned long maxocc = grid->getMaximalOccupiedSecondTable();


        for (unsigned long i = 0; i < maxocc; i++) {
            IndexDimension Index = grid->getIndexOfTable(i);

            if (Index.getDepth(d) == t) {
                double value = fine.getValue(i);
/*
                    if (Index.isAtLeftBoundary(d))
                        value = fine.getValue(Index);

                    if (Index.isAtRightBoundary(d))
                        value = fine.getValue(Index);
*/

                if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d))) {
                    value = value + 0.5 * (fine.getValue(Index.nextRight(d)) + fine.getValue(Index.nextLeft(d)));

                }
                //coarse.setValue(Index, value);
                coarse = value | Index;


            }


        }


    }

    void prolongation1D_inplace(VectorSparseG &vec, int &t, int &d) {
        AdaptiveSparseGrid_Base *grid = vec.getSparseGrid();


        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            if (Tlocal.at(d) == t) {
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                do {
                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    double value = vec.getValue(i);
                    if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                        value = value + 0.5 * (vec.getValue(Index.nextRight(d)) + vec.getValue(Index.nextLeft(d)));
                    //coarse = value | Index;
                    vec.setValue(i, value);

                } while (inneriter.next());
            }

        } while (iter.next());
    }



    inline void prolongation(VectorSparseG &coarse, VectorSparseG &fine, Depth &Tcoarse, Depth &Tfine) {
        fine = coarse;
        int tcoarse;
        int tfine;

        VectorSparseG testvec(*coarse.getSparseGrid());
        testvec = coarse;
        VectorSparseG testvec2(*coarse.getSparseGrid());

        for (int d = 0; d < DimensionSparseGrid; d++) {

            tcoarse = Tcoarse.at(d) + 1;
            tfine = Tfine.at(d);
            for (int t = tcoarse; t <= tfine; ++t) {
                //prolongation1D(fine, fine, t, d);
                //prolongation1D_inplace(fine,t,d);
                prolongation1D_inplace_singleHash(list,fine,t,d);
     /*           testvec2 = testvec-fine;
                if(L_infty(testvec2)>1e-10) {
                    for(unsigned long k=0; k < testvec.getSparseGrid()->getMaximalOccupiedSecondTable(); k++){
                        if(abs(testvec2.getValue(k))>1e-10) {
                            testvec2.getSparseGrid()->getIndexOfTable(k).Print();
                            cout << testvec2.getValue(k) << endl;
                            cout <<  fine.getValue(k) << endl;
                            cout << testvec.getValue(k) << endl;
                            exit(EXIT_FAILURE);
                        }
                    }

                }*/
            }
        }
    };

    static void prolongation1D_inplace_singleHash(ListOfDepthOrderedSubgrids& list, VectorSparseG &vec,const int &t,const int &d) {
        AdaptiveSparseGrid_Base *grid = vec.getSparseGrid();


        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {

            Depth Tlocal = iter.getDepth();


            if (Tlocal.at(d) == t) {


                vector<SingleDepthHashGrid*> fineDepthGrids = grid->getMultiDepthHashGrid()->getGridsForDepthInDirection(Tlocal,d);
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                do {
                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    double value = vec.getValue(i);
                    if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d))){
                        auto right = Index.nextRight(d);
                        auto left = Index.nextLeft(d);
                        // auto* gridRight = fineDepthGrids[right.getDepth(d)];
                        // auto* gridLeft = fineDepthGrids[left.getDepth(d)];
                        // auto valRight = vec.getValue(right,gridRight);
                        // auto valLeft = vec.getValue(left,gridLeft);
                        value = value + 0.5 * (vec.getValue(right,fineDepthGrids[right.getDepth(d)]) + vec.getValue(left,fineDepthGrids[left.getDepth(d)]));

                        // value = value + 0.5 * (vec.getValue(Index.nextRight(d)) + vec.getValue(Index.nextLeft(d)));
                        // auto p= vec.getValue(right,fineDepthGrids[right.getDepth(d)]);
                    }
                    //coarse = value | Index;
                    vec.setValue(i, value);
                    // cout << Depth(Index.nextRight(d)).at(d)<<", \t" <<Depth(Index.nextLeft(d)).at(d) << endl;

                    // cout << "Left: " ;
                    //     Depth(Index.nextLeft(d)).Print();
                    // cout << "Right: " ;
                    //     Depth(Index.nextRight(d)).Print();
                } while (inneriter.next());


            }

            // cout << "iter" << endl;
        } while (iter.next());
        // cout << "end" << endl;


    }

    inline void restriction1D(VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {

        AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
        unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
        double value;
        double valueL;
        double valueR;
        for (unsigned long i = 0; i < maxocc; i++) {
            IndexDimension Index = grid->getIndexOfTable(i);

            if (Index.getDepth(d) <= t) {

                value = fine.getValue(i);
                valueL = 0.0;
                valueR = 0.0;
                if (t >= 0) {
                    if (!(Index.isAtRightBoundary(d))) valueL = 0.5 * fine.getValue(Index.nextRight(d, t + 1));
                    if (!(Index.isAtLeftBoundary(d))) valueR = 0.5 * fine.getValue(Index.nextLeft(d, t + 1));
                }


                value = value + valueL + valueR;
                //coarse.setValue(Index, value);
                coarse.setValue(i, value);
                //coarse = value | Index;
            }


        }

    };


    inline void restriction1D_inverted_inplace( const VectorSparseG &fine,const VectorSparseG &coarse, int t, int d) {

        double value;

        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            // SingleDepthHashGrid& coarseDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
            if(Tlocal.at(d)<=t){

                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                // SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    unsigned long i = inneriter.geti();
                    double value = fine.getValue(i);
                    coarse.setValue(i,value);

                } while (inneriter.next());

            }

            if (Tlocal.at(d) == t + 1) {

                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                // SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    value = 0.5 * fine.getValue(i);

                    if (t >= 0) {

                        IndexDimension rightIndex = Index.nextRight(d);
                        if (!Index.isAtRightBoundary(d)){
                            // vec.getValue(rightIndex,value,singleDepthHashGrid)
                           coarse.addToValue(rightIndex,value);
                            // value += vec.getValue(rightIndex);
                        }


                        IndexDimension leftIndex = Index.nextLeft(d);
                        if (!Index.isAtLeftBoundary(d)){
                            // vec.getValue(rightIndex,value,singleDepthHashGrid)
                            coarse.addToValue(leftIndex,value);
                            // value += vec.getValue(rightIndex);
                        }
                    }
                } while (inneriter.next());
            }
        } while (iter.next());

    };


    inline void restriction1D_inverted_inplace( const VectorSparseG &vec, int t, int d) {

        double value;

        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            if (Tlocal.at(d) == t + 1) {

                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                // SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    value = 0.5 * vec.getValue(i);

                    if (t >= 0) {

                        IndexDimension rightIndex = Index.nextRight(d);
                        if (!Index.isAtRightBoundary(d)){
                            // vec.getValue(rightIndex,value,singleDepthHashGrid)
                            vec.addToValue(rightIndex,value);
                            // value += vec.getValue(rightIndex);
                        }


                        IndexDimension leftIndex = Index.nextLeft(d);
                        if (!Index.isAtLeftBoundary(d)){
                            // vec.getValue(rightIndex,value,singleDepthHashGrid)
                            vec.addToValue(leftIndex,value);
                            // value += vec.getValue(rightIndex);
                        }
                    }
                } while (inneriter.next());
            }
        } while (iter.next());

    };


    inline void DualSpaceTransformationNeumann(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T) {
        SingleDepthHashGrid& depthGrid = Ax_hier.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
        const auto& mapping = depthGrid._mapPosToGridPos;

        for (size_t i = 0; i < mapping.size(); i++)
        {

            IndexDimension I = depthGrid._map.getIndexOfTable(i);
            double coeff = 0.0;
            IndexDimension J;
            if(Ax_hier.getSparseGrid()->getActiveTable()[mapping[i]]) {
                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    J = I.nextFive_Neumann(&mc, T, &basis_coeff);
                    if (mc.goon())
                        coeff = coeff + Ax_hier.getValue(J) * basis_coeff;
                }

                Ax.setValue(mapping[i], coeff);
            }
        }

    }

    inline void DualSpaceTransformationDirichlet(VectorSparseG &g, VectorSparseG &Ax, Depth &T) {
        SingleDepthHashGrid& depthGrid = g.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
        const auto& mapping = depthGrid._mapPosToGridPos;

        for (size_t i = 0; i < mapping.size(); i++)
        {

            IndexDimension I = depthGrid._map.getIndexOfTable(i);
            Richtung richtung[DimensionSparseGrid];
            if(I.isRandnah(richtung)){

            }
            double coeff = 0.0;
            IndexDimension J;
            if(g.getSparseGrid()->getActiveTable()[mapping[i]]) {
                for (MultiDimCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    J = I.nextThree_boundary(&mc, T);
                    if (mc.goon())
                        coeff = coeff + g.getValue(J) * basis_coeff;
                }

                Ax.setValue(mapping[i], coeff);
            }
        }

    }

    class iterator{
    public:
        iterator(Richtung* richtung_){
            maxShift = 1;
            for(int d=0; d<DimensionSparseGrid; d++){
                richtung[d]=richtung_[d];
                if(richtung[d]!=Mitte)
                    maxShift *=2;

            }
        }
        bool goon(){
                return false;
        }
    private:
        Richtung richtung[DimensionSparseGrid];
        int maxShift;
        int shiftnumber;
    };

};

///////////////////////////////////////////////////////////////////

template<class Problem>
void MatrixVectorInhomogen::multiplication(VectorSparseG &prew, VectorSparseG &Ax, Problem &matrixProblem) {

    Ax = 0.0;
    nodal = 0.0;

    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();

    int k = 1;
    grid->WorkOnHangingNodes = true;

    ListOfDepthOrderedSubgrids::iterator iter(list);
    iter.gotoBegin();
    do {
        Depth T = iter.getDepth();
        u = 0.0;
        calcNodal(prew, T);
        nodal.setMultiLevelValues(u, T);
    } while (iter.next());

    grid->WorkOnHangingNodes = false;


    for (CasesIterator iter; iter.goon(); ++iter) {
        Ax_neu = 0.0;
        bool *restrictions = iter.getcase();

        grid->WorkOnHangingNodes = true;

        Ax_neu = 0.0;
        gM = 0.0;


        CaseFunction(prew, restrictions, matrixProblem);



        grid->WorkOnHangingNodes = false;
        Ax = Ax + Ax_neu;

        k++;
    }

}


template<class Problem>
void MatrixVectorInhomogen::multiplicationRHSHomogen(VectorSparseG &prew, VectorSparseG &Ax, Problem &matrixProblem) {

    Ax = 0.0;
    nodal = 0.0;
    AdaptiveSparseGrid *grid = prew.getSparseGrid();

    int k = 1;
    grid->WorkOnHangingNodes = true;


    ListOfDepthOrderedSubgrids::iterator iter(list);
    iter.gotoBegin();
    do {
        Depth T = iter.getDepth();
        u = 0.0;
        calcNodal(prew, T);
        nodal.setMultiLevelValues2(u, T);
    } while (iter.next());

    grid->WorkOnHangingNodes = false;


    for (CasesIterator iter; iter.goon(); ++iter) {
        Ax_neu = 0.0;
        bool *restrictions = iter.getcase();

        grid->WorkOnHangingNodes = true;

        u = 0.0;
        g = 0.0;
        z = 0.0;
        Ax_neu = 0.0;
        gM = 0.0;


        CaseFunctionTransformation(prew, restrictions, matrixProblem);

        grid->WorkOnHangingNodes = false;
        Ax = Ax + Ax_neu;

        k++;
    }

}

template<class Problem>
void MatrixVectorInhomogen::CaseFunction(VectorSparseG &prew, bool *restrictions, Problem &matrixProblem) {

    list.SortiereTiefenBoundary(restrictions);
    ListOfDepthOrderedSubgrids::iterator iter(list);

    std::list<Depth>* sortedDepths = list.getSortierteTiefen(); // TODO maybe not deference to prevent useless copy of all elements

    if (sortedDepths->size() == 0) return;

    P+=nodal;

    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();

    for(int d=0; d < DimensionSparseGrid; d++) {
        if(restrictions[d]==0)
            for (Depth T: *sortedDepths) {
                int t = T.at(d)+1;
                Depth T_fine = T;
                T_fine.set(t,d);

                if(list.isIncluded(T_fine)) {
                    u_new = 0.0;
                    u_new.setMultiLevelValues2(P, T);

                    prolongation1D_inplace(u_new, t, d);
                    P.addMultiLevelValues(u_new, T_fine);
                }

            }
    }


    for (Depth T: *sortedDepths) {
        grid->WorkOnHangingNodes = true;
        u = 0.0;


        u.setMultiLevelValues2(P,T);



        z = 0.0;

        matrixProblem.initialize(T);
        applyProblemGlobal(u, z, matrixProblem, T);

        g = 0.0;
        g.setMultiLevelValues(gM, T);

        z = z + g;
        g = z;






        Depth Tcoarse(T);
        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {


                Tcoarse.set(T.at(d) - 1, d);
                if (list.isIncluded(Tcoarse)) {
                    g = 0.0;


                    restriction1D_inverted_inplace(z, g, Tcoarse.at(d), d);


                    z = 0.0;
                    z.setMultiLevelValues(gM, Tcoarse);
                    z = z + g;
                    gM.setMultiLevelValues(g, Tcoarse);
                }


                }


        }


        if (list.isIncluded(Tcoarse)) {
            grid->WorkOnHangingNodes = false;

            DualSpaceTransformationNeumann(g, Ax_neu, Tcoarse);
        }

    }

};


template<class Problem>
void MatrixVectorInhomogen::CaseFunctionTransformation(VectorSparseG &prew, bool *restrictions, Problem &matrixProblem) {

    list.SortiereTiefenBoundary(restrictions);
    ListOfDepthOrderedSubgrids::iterator iter(list);

    std::list<Depth>* sortedDepths = list.getSortierteTiefen(); // TODO maybe not deference to prevent useless copy of all elements

    if (sortedDepths->size() == 0) return;

    P+=nodal;

    AdaptiveSparseGrid *grid = prew.getSparseGrid();
    for (Depth T: *sortedDepths) {

        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d] == 0) {
                int t = T.at(d)+1;
                Depth T_fine = T;
                T_fine.set(t,d);
            }
        }
    }

    for(int d=0; d < DimensionSparseGrid; d++) {
        if(restrictions[d]==0)
            for (Depth T: *sortedDepths) {
                int t = T.at(d)+1;
                Depth T_fine = T;
                T_fine.set(t,d);

                if(list.isIncluded(T_fine)) {
                    u_new = 0.0;
                    u_new.setMultiLevelValues2(P, T);

                    prolongation1D_inplace(u_new, t, d);
                    P.addMultiLevelValues(u_new, T_fine);
                }

            }
    }


    for (Depth T: *sortedDepths) {
        grid->WorkOnHangingNodes = true;
        u = 0.0;


        u.setMultiLevelValues2(P,T);



        z = 0.0;

        matrixProblem.initialize(T);
        applyProblemGlobal(u, z, matrixProblem, T);

        z.addMultiLevelValues2(gM,T);






        Depth Tcoarse(T);
        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {

                Tcoarse.set(T.at(d) - 1, d);
                if (list.isIncluded(Tcoarse)) {

                    restriction1D_inverted_inplace(z,Tcoarse.at(d), d);
                    g=z;
                    z.addMultiLevelValues2(gM,Tcoarse);
                    gM.setMultiLevelValues2(g, Tcoarse);


                }

            }
        }


        if (list.isIncluded(Tcoarse)) {
            grid->WorkOnHangingNodes = false;
            g = 0.0;
            DualSpaceTransformationNeumann(z, g, Tcoarse);
            DualSpaceTransformationDirichlet(g,Ax_neu,Tcoarse);
        }

    }

};


void MatrixVectorInhomogen::calcNodal(VectorSparseG &prew, Depth &Tiefe) {
    SubgridFixedDepth::iterator iter(*list.getSubgrid(Tiefe));
    iter.gotobegin();
    do {

        IndexDimension Index = iter.getPoint();
        double coeff = prew.getValue(Index);
        // Depth Tneu(Index);
        double basis_coeff = 1.0;
        for (MultiDimFiveCompass mc; mc.hasNext(); ++mc) {

            basis_coeff = 1.0;
            IndexDimension J = Index.nextFive_Neumann(&mc, Tiefe.returnTiefen(), &basis_coeff);

            double val = u.getValue(J) + coeff * basis_coeff;

            //hier.setValue(J, val);
            //u = val | J;
            u.setValue(J, val);


        }


    } while (iter.next());
}

#endif //SGRUN_MATRIXVECTORINHOMOGEN_H_H
