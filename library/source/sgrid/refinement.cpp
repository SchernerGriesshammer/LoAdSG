//
// Created by scherner on 27.07.21.
//

#include "refinement.h"
#include "../stencils/Stencil.h"
#include "../applications/cg_method.h"
//#include "../applications/norms.h"


bool AdaptiveRefinement(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps) {

    VectorSparseG hier(sgrid);
    hier = nodal;
    CalcHierarchicalBasis(hier);

    int rank;




    bool changed = false;
    unsigned long end = sgrid.getMaximalOccupiedSecondTable();
    for (unsigned long i = 0; i < end; i++) {
        if (sgrid.getSecondTable()[i] != 0) {
            IndexDimension I = sgrid.getIndexOfTable(i);
            Depth T(I);
            if (abs(hier.getValue(i)) > eps && sgrid.getActiveTable()[i]) {
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (I.getIndex(d) != 0) {
                        if (sgrid.AddPoint(I.leftSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.leftSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                    if (I.getIndex(d) != 1) {
                        if (sgrid.AddPoint(I.rightSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.rightSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                }

            }
        }
    }

    return changed;
}

bool AdaptiveRefinement(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps, int level) {

    VectorSparseG hier(sgrid);
    hier = nodal;
    CalcHierarchicalBasis(hier);


    int rank;




    bool changed = false;
    unsigned long end = sgrid.getMaximalOccupiedSecondTable();
    for (unsigned long i = 0; i < end; i++) {
        if (sgrid.getSecondTable()[i] != 0) {
            IndexDimension I = sgrid.getIndexOfTable(i);
            Depth T(I);
            if (abs(hier.getValue(i)) > eps && sgrid.getActiveTable()[i]) {
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (I.getIndex(d) != 0) {
                        if (sgrid.AddPoint(I.nextLeft(d,level)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k,I.nextLeft(d,level));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                    if (I.getIndex(d) != 1) {
                        if (sgrid.AddPoint(I.nextRight(d,level)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.nextRight(d,level));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                }

            }
        }
    }

    return changed;
}

bool AdaptiveRefinement(AdaptiveSparseGrid &oldgrid,AdaptiveSparseGrid &newgrid,  VectorSparseG &nodal, double eps) {



    VectorSparseG hier(*nodal.getSparseGrid());
    hier = nodal;

    CalcHierarchicalBasis(hier);
    bool changed = false;
    for (unsigned long i = 0; i < nodal.getSparseGrid()->getMaximalOccupiedSecondTable(); i++) {


        if (abs(hier.getValue(i)) > eps && nodal.getSparseGrid()->getActiveTable()[i]) {
                IndexDimension I = nodal.getSparseGrid()->getIndexOfTable(i);
                newgrid.AddPoint(I);
                unsigned long k;

                Depth T(I);


                for (int d = 0; d < DimensionSparseGrid; d++) {
                       if(T.at(d)>0) {
                               if (!nodal.getSparseGrid()->occupied(k, I.leftSon(d)) &&
                                   !nodal.getSparseGrid()->occupied(k, I.rightSon(d))) {
                                   newgrid.AddPoint(I.leftSon(d));
                                   newgrid.AddPoint(I.rightSon(d));
                                   changed = true;
                               }
                       }
                }

        }

        if (abs(hier.getValue(i)) <= eps) hier.setValue(i,0.0);

    }

    nodal=hier;
    CalcNodalBasis(nodal);

    return changed;
}



bool AdaptiveRefinementL2(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps) {

    VectorSparseG hier(sgrid);
    hier = nodal;
    CalcHierarchicalBasis(hier);

    bool changed = false;
    for (unsigned long i = 0; i < sgrid.getMaximalOccupiedSecondTable(); i++) {
        if (sgrid.getSecondTable()[i] != 0 &&sgrid.getActiveTable()[i]) {
            IndexDimension I = sgrid.getIndexOfTable(i);
            Depth T(I);


            double value = 1.0;
            double h = 0.0;
            for (int d = 0; d < DimensionSparseGrid; d++) {
                int exp = T.at(d);
                if(exp == 0) {
                    h = 1.0;
                    value *=2*h/6;
                }
                else {

                    int value2 = 1 << exp;
                    h = 1 / double(value2);

                }
                value *= 4*h/6;
            }
            double l2 = hier.getValue(i)*value;


            if (abs(l2) > eps && sgrid.getActiveTable()[i]) {
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (I.getIndex(d) != 0) {
                        if (sgrid.AddPoint(I.leftSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.leftSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                    if (I.getIndex(d) != 1) {
                        if (sgrid.AddPoint(I.rightSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.rightSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                }

            }
        }
    }

    return changed;
}

bool AdaptiveRefinementEnergy(AdaptiveSparseGrid& newgrid, VectorSparseG &prew, double eps) {

    Poisson poissonStencil(*prew.getSparseGrid());
    VectorSparseG zPoisson(prew.getSparseGrid());
    zPoisson=prew;
    Preconditioning<Poisson> pPoisson(zPoisson,poissonStencil);

    HelmHoltz helmHoltz(*prew.getSparseGrid());
    VectorSparseG zHelm(prew.getSparseGrid());
    zHelm = prew;
    Preconditioning<HelmHoltz> pHelm(zHelm,helmHoltz);
    bool changed = false;



    for (unsigned long i = 0; i < prew.getSparseGrid()->getMaximalOccupiedSecondTable(); i++) {

        if(prew.getSparseGrid()->getActiveTable()[i]){
            double norm = prew.getValue(i)*prew.getValue(i)*(pHelm.getValue(i)+pPoisson.getValue(i));




            if (norm > eps) {
                IndexDimension I = prew.getSparseGrid()->getIndexOfTable(i);
                newgrid.AddPoint(I);
                unsigned long k;

                Depth T(I);


                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if(T.at(d)>0) {
                        if (!prew.getSparseGrid()->occupied(k, I.leftSon(d)) &&
                            !prew.getSparseGrid()->occupied(k, I.rightSon(d))) {
                            newgrid.AddPoint(I.leftSon(d));
                            newgrid.AddPoint(I.rightSon(d));
                            changed = true;
                        }
                    }
                }

            }

        }


    }




    return changed;
}

double integral2(double a, double b, double mp, double mq, double cp, double cq){
    double term1 = (1.0/3.0*b*b*b*mp*mq)-(1.0/3.0*a*a*a*mp*mq);
    double term2 = (0.5*b*b*mp*cq)-(0.5*a*a*mp*cq);
    double term3 = (0.5*b*b*mq*cp)-(0.5*a*a*mq*cp);
    double term4 = (b-a)*cp*cq;

    return term1+term2+term3+term4;

}
double CalcIntegral(IndexDimension p, IndexDimension q, Depth T, int d){

    double integral_value=0.0;


    //Steigung m
    double mLeft;
    double mRight;


    // Shift
    double cpLeft;
    double cqLeft;
    double cpRight;
    double cqRight;

    // Meshsize
    double h;


    //Coordinates
    double coordp;
    double coordq;

    // support
    double suppleftP;
    double supprightP;

    double suppleftQ;
    double supprightQ;



        mLeft=pow(2,T.at(d));


        h= 1.0/mLeft;

        mRight= -1.0*mLeft;
        coordp=p.coordinate(d);
        cpLeft=-(coordp-h)/h;
        cpRight=+(coordp+h)/h;



        suppleftP=coordp-h;
        supprightP=coordp+h;







        coordq=q.coordinate(d);
        cqLeft=-(coordq-h)/h;
        cqRight=+(coordq+h)/h;


        suppleftQ=coordq-h;
        supprightQ=coordq+h;



    double leftP_leftQ;
    double leftP_rightQ;
    double rightP_leftQ;
    double rightP_rightQ;


    double sum_inner=0.0;
    double prod=1.0;


    bool docase1 = false;
    bool docase2 = false;
    bool docase3 = false;
    bool docase4 = false;
    if (coordq==coordp) {
        docase1 = true;
        docase4 = true;
    } else if(coordp<coordq) docase3 = true;
    else if(coordq>coordp) docase2 = true;


    sum_inner = 0.0;
    leftP_leftQ = 0.0;
    leftP_rightQ = 0.0;
    rightP_leftQ = 0.0;
    rightP_rightQ = 0.0;

    // int d/dx vp * d/dx vp mit x==dinner

    //1st case: integral over leftp cap leftq
    double intLp = suppleftP;
    double intRp = coordp;

    double intLq = suppleftQ;
    double intRq = coordq;


    if (docase1) {


        if (intRp < intLq) {
            leftP_leftQ = 0.0;
        } else if (intRq < intLp) {
            leftP_leftQ = 0.0;
        } else {
            leftP_leftQ = integral2(max(intLq, intLp), min(intRp, intRq), mLeft, mLeft,
                                   cpLeft, cqLeft);
        }
    }





    //2nd case: integral over leftp cap rightq
    intLp = suppleftP;
    intRp = coordp;

    intLq = coordq;
    intRq = supprightQ;

    if (docase2) {
        if (intRp < intLq) {
            leftP_rightQ = 0.0;
        } else if (intRq < intLp) {
            leftP_rightQ = 0.0;
        } else {
            leftP_rightQ = integral2(max(intLq, intLp), min(intRp, intRq), mLeft,
                                    mRight, cpLeft, cqRight);
        }
    }

    //3rd case: integral over rightp cap leftq
    intLp = coordp;
    intRp = supprightP;

    intLq = suppleftQ;
    intRq = coordq;

    if (docase3) {

        if (intRp < intLq) {
            rightP_leftQ = 0.0;
        } else if (intRq < intLp) {
            rightP_leftQ = 0.0;
        } else {
            rightP_leftQ = integral2(max(intLq, intLp), min(intRp, intRq), mRight,
                                    mLeft, cpRight, cqLeft);
        }
    }

    //4th case: integral over rightp cap rightq
    intLp = coordp;
    intRp = supprightP;

    intLq = coordq;
    intRq = supprightQ;

    if (docase4) {

        if (intRp < intLq) {
            rightP_rightQ = 0.0;
        } else if (intRq < intLp) {
            rightP_rightQ = 0.0;
        } else {
            rightP_rightQ = integral2(max(intLq, intLp), min(intRp, intRq), mRight,
                                     mRight, cpRight, cqRight);
        }
    }


    sum_inner = leftP_leftQ + leftP_rightQ + rightP_leftQ + rightP_rightQ;






    return  sum_inner;


}

double CalcL2ofPrew(IndexDimension I){

    Depth T(I);

    double prod= 1.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        double sum = 0.0;
        int exp = T.at(d);
        if(exp == 0) {
            sum+= CalcIntegral(I,I,T,d);
        } else {
            if (I.isLinksRandNah(d)) {
                IndexDimension IR = I.nextRight(d, T.at(d));
                IndexDimension IRR = IR.nextRight(d, T.at(d));
                sum += 0.9 * 0.9 * CalcIntegral(I, I, T, d);
                sum += 2.0 * 0.9 * -0.6 * CalcIntegral(I, IR, T, d);
                sum += -0.6 * -0.6 * CalcIntegral(IR, IR, T, d);
                sum += 2.0 * -0.6 * 0.1 * CalcIntegral(IR, IRR, T, d);
            } else if (I.isRechtsRandNah(d)) {
                IndexDimension IL = I.nextLeft(d, T.at(d));
                IndexDimension ILL = IL.nextLeft(d, T.at(d));
                sum += 0.9 * 0.9 * CalcIntegral(I, I, T, d);
                sum += 2.0 * 0.9 * -0.6 * CalcIntegral(I, IL, T, d);
                sum += -0.6 * -0.6 * CalcIntegral(IL, IL, T, d);
                sum += 2.0 * -0.6 * 0.1 * CalcIntegral(IL, ILL, T, d);
            } else {
                IndexDimension IR = I.nextRight(d, T.at(d));
                IndexDimension IRR = IR.nextRight(d, T.at(d));
                IndexDimension IL = I.nextLeft(d, T.at(d));
                IndexDimension ILL = IL.nextLeft(d, T.at(d));
                sum += 0.1 * 0.1 * CalcIntegral(ILL, ILL, T, d);
                sum +=2.0* 0.1 * -0.6 * CalcIntegral(ILL, IL, T, d);
                sum += 0.6 * 0.6 * CalcIntegral(IL, IL, T, d);
                sum += 2.0-0.6 * CalcIntegral(IL, I, T, d);
                sum += CalcIntegral(I, I, T, d);
                sum += 2.0*-0.6* CalcIntegral(I, IR, T, d);
                sum += 0.6 * 0.6 * CalcIntegral(IR, IR, T, d);
                sum += 2.0 * -0.6 * 0.1 * CalcIntegral(IR, IRR, T, d);
                sum += 0.1 * 0.1 * CalcIntegral(IRR, IRR, T, d);


            }

        }
        prod*=sum;
    }
    return prod;
}

/**
 * Noch nicht getestet!!!!
 *
 * @param sgrid
 * @param prew
 * @param eps
 * @return
 */
bool AdaptiveRefinementL2_Prew(AdaptiveSparseGrid &sgrid, VectorSparseG &prew, double eps) {

    double l2;

    bool changed = false;
    for (unsigned long i = 0; i < sgrid.getMaximalOccupiedSecondTable(); i++) {
        if (sgrid.getSecondTable()[i] != 0 &&sgrid.getActiveTable()[i]) {
            IndexDimension I = sgrid.getIndexOfTable(i);
            Depth T(I);
            l2 =sqrt(CalcL2ofPrew(I));
            l2 = l2*prew.getValue(i);

            if (abs(l2) > eps && sgrid.getActiveTable()[i]) {
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (I.getIndex(d) != 0) {
                        if (sgrid.AddPoint(I.leftSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.leftSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                    if (I.getIndex(d) != 1) {
                        if (sgrid.AddPoint(I.rightSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid.occupied(k, I.rightSon(d));
                        if (sgrid.getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid.getActiveTable()[k] = true;
                        }
                    }
                }

            }
        }
    }

    return changed;
}


bool Refinement::apply(VectorSparseG &nodal, double eps) {

    AdaptiveSparseGrid* sgrid = nodal.getSparseGrid();
    VectorSparseG hier(*sgrid);
    hier = nodal;

    calcPrewByNodal(hier,nodal);



    bool changed = false;
    for (unsigned long i = 0; i < sgrid->getMaximalOccupiedSecondTable(); i++) {
        if (sgrid->getSecondTable()[i] != 0 && refined.getValue(i)<1.0) {
            IndexDimension I = sgrid->getIndexOfTable(i);
            Depth T(I);
            if (abs(hier.getValue(I)) > eps && sgrid->getActiveTable()[i] && T > 0) {
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (I.getIndex(d) != 0) {
                        if (sgrid->AddPoint(I.leftSon(d))){
                            changed = true;
                        }

                        unsigned long k;
                        sgrid->occupied(k, I.leftSon(d));
                        if (sgrid->getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid->getActiveTable()[k] = true;

                        }
                    }
                    if (I.getIndex(d) != 1) {
                        if (sgrid->AddPoint(I.rightSon(d)))
                            changed = true;
                        unsigned long k;
                        sgrid->occupied(k, I.rightSon(d));
                        if (sgrid->getActiveTable()[k] == 0) {
                            changed = true;
                            sgrid->getActiveTable()[k] = true;
                        }
                    }
                }
                refined.setValue(i,1.0);
            }
        }
    }

    return changed;

}
