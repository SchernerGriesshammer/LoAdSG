//
// Created by to35jepo on 3/13/23.
//

#include "RHS.h"
void CalcHierarchicalBasisForRHS(VectorSparseG &u) {
    DepthList depthList(*u.getSparseGrid());
    bool prolongations[DimensionSparseGrid]={0};
    depthList.sortDepthsNeumann(prolongations);
    std::list<Depth>* sortedDepths = depthList.getSortierteTiefen();
    sortedDepths->reverse();
    for(int d=0; d<DimensionSparseGrid;d++){
        for(Depth T: *sortedDepths) {
            if(!(T>>0)) {
                if(T.at(d)>0){

                    SingleDepthHashGrid &depthGrid =u.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                    const auto &mapping = depthGrid._mapPosToGridPos;
                    if (depthGrid.getNumberOfEntries() > 0) {
                        for (size_t i = 0; i < mapping.size(); i++){

                            IndexDimension INodal = depthGrid._map.getIndexOfTable(i);
                            IndexDimension LeftI = INodal.nextLeft(d);
                            IndexDimension RightI = INodal.nextRight(d);

                            unsigned long Il,Ir;
                            u.getSparseGrid()->occupied(Il,LeftI);
                            u.getSparseGrid()->occupied(Ir,RightI);

                            double val = u.getValue(mapping[i]) - 0.5 * (u.getValue(Il) + u.getValue(Ir));
                            u.setValue(mapping[i],val);
                        }
                    }


                }
            }

        }
    }
}

double integral_gradient(double a, double b, double mp, double mq){

    return (b-a)*mp*mq;
}
double integral(double a, double b, double mp, double mq, double cp, double cq){
    double term1 = (1.0/3.0*b*b*b*mp*mq)-(1.0/3.0*a*a*a*mp*mq);
    double term2 = (0.5*b*b*mp*cq)-(0.5*a*a*mp*cq);
    double term3 = (0.5*b*b*mq*cp)-(0.5*a*a*mq*cp);
    double term4 = (b-a)*cp*cq;

    return term1+term2+term3+term4;

}
double CalcIntegral(IndexDimension p, IndexDimension q, Depth TP, Depth TQ){



    double mpLeft[DimensionSparseGrid];
    double mqLeft[DimensionSparseGrid];
    double mpRight[DimensionSparseGrid];
    double mqRight[DimensionSparseGrid];

    double cpLeft[DimensionSparseGrid];
    double cqLeft[DimensionSparseGrid];
    double cpRight[DimensionSparseGrid];
    double cqRight[DimensionSparseGrid];

    double hp[DimensionSparseGrid];
    double hq[DimensionSparseGrid];

    double coordp[DimensionSparseGrid];
    double coordq[DimensionSparseGrid];

    double suppleftP[DimensionSparseGrid];
    double supprightP[DimensionSparseGrid];

    double suppleftQ[DimensionSparseGrid];
    double supprightQ[DimensionSparseGrid];

    for(int d=0; d<DimensionSparseGrid; d++){

        mpLeft[d]=POW2(TP.at(d));


        hp[d]= 1.0/mpLeft[d];

        mpRight[d]= -1.0*mpLeft[d];
        coordp[d]=p.coordinate(d);
        cpLeft[d]=-(coordp[d]-hp[d])/hp[d];
        cpRight[d]=+(coordp[d]+hp[d])/hp[d];



        suppleftP[d]=coordp[d]-hp[d];
        supprightP[d]=coordp[d]+hp[d];






        mqLeft[d]=POW2(TQ.at(d));
        hq[d]= 1.0/mqLeft[d];
        mqRight[d]= -1.0*mqLeft[d];
        coordq[d]=q.coordinate(d);
        cqLeft[d]=-(coordq[d]-hq[d])/hq[d];
        cqRight[d]=+(coordq[d]+hq[d])/hq[d];


        suppleftQ[d]=coordq[d]-hq[d];
        supprightQ[d]=coordq[d]+hq[d];

    }

    double leftP_leftQ;
    double leftP_rightQ;
    double rightP_leftQ;
    double rightP_rightQ;

    double sum_outer=0.0;
    double sum_inner;
    double prod;


    for(int douter=0; douter<DimensionSparseGrid; douter++) {

        {

            prod = 1.0;
            for (int dinner = 0; dinner < DimensionSparseGrid; dinner++) {




                bool docase1=true;
                bool docase2=true;
                bool docase3=true;
                bool docase4=true;
                if(q.getIndex(dinner)==0){
                    docase1=false;
                    docase3=false;
                }
                if(q.getIndex(dinner)==1){
                    docase2= false;
                    docase4= false;
                }

                sum_inner = 0.0;
                leftP_leftQ = 0.0;
                leftP_rightQ = 0.0;
                rightP_leftQ = 0.0;
                rightP_rightQ = 0.0;

                // int d/dx vp * d/dx vp mit x==dinner

                //1st case: integral over leftp cap leftq
                double intLp = suppleftP[dinner];
                double intRp = coordp[dinner];

                double intLq = suppleftQ[dinner];
                double intRq = coordq[dinner];


                if(docase1) {


                    if (intRp < intLq) {
                        leftP_leftQ = 0.0;
                    } else if (intRq < intLp) {
                        leftP_leftQ = 0.0;
                    } else {
                        if (dinner == douter) {
                            leftP_leftQ = integral_gradient(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner],
                                                            mqLeft[dinner]);
                        } else {
                            leftP_leftQ = integral(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner], mqLeft[dinner],
                                                   cpLeft[dinner], cqLeft[dinner]);
                        }
                    }



                }

                //2nd case: integral over leftp cap rightq
                intLp = suppleftP[dinner];
                intRp = coordp[dinner];

                intLq = coordq[dinner];
                intRq = supprightQ[dinner];

                if(docase2) {

                    if (intRp < intLq) {
                        leftP_rightQ = 0.0;
                    } else if (intRq < intLp) {
                        leftP_rightQ = 0.0;
                    } else {
                        if (dinner == douter) {
                            leftP_rightQ = integral_gradient(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner],
                                                             mqRight[dinner]);

                        } else {
                            leftP_rightQ = integral(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner],
                                                    mqRight[dinner], cpLeft[dinner], cqRight[dinner]);

                        }
                    }

                }

                //3rd case: integral over rightp cap leftq
                intLp = coordp[dinner];
                intRp = supprightP[dinner];

                intLq = suppleftQ[dinner];
                intRq = coordq[dinner];

                if(docase3) {

                    if (intRp < intLq) {
                        rightP_leftQ = 0.0;
                    } else if (intRq < intLp) {
                        rightP_leftQ = 0.0;
                    } else {
                        if (dinner == douter) {
                            rightP_leftQ = integral_gradient(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                                             mqLeft[dinner]);
                        } else {
                            rightP_leftQ = integral(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                                    mqLeft[dinner], cpRight[dinner], cqLeft[dinner]);
                        }
                    }
                }

                //4th case: integral over rightp cap rightq
                intLp = coordp[dinner];
                intRp = supprightP[dinner];

                intLq = coordq[dinner];
                intRq = supprightQ[dinner];

                if(docase4) {

                    if (intRp < intLq) {
                        rightP_rightQ = 0.0;
                    } else if (intRq < intLp) {
                        rightP_rightQ = 0.0;
                    } else {
                        if (dinner == douter) {
                            rightP_rightQ = integral_gradient(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                                              mqRight[dinner]);
                        } else {
                            rightP_rightQ = integral(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                                     mqRight[dinner], cpRight[dinner], cqRight[dinner]);
                        }
                    }
                }




                sum_inner = leftP_leftQ + leftP_rightQ + rightP_leftQ + rightP_rightQ;


                prod *= sum_inner;


            }
            sum_outer += prod;
        }


    }

    return  sum_outer;

}




double integral_coeff(double a, double b, double mp, double mq, double cp, double cq){
    double term1a = (1.0/3.0*b*b*b*mp*mq)-(1.0/3.0*a*a*a*mp*mq);
    double term2a = (0.5*b*b*mp*cq)-(0.5*a*a*mp*cq);
    double term3a = (0.5*b*b*mq*cp)-(0.5*a*a*mq*cp);
    double term4a = (b-a)*cp*cq;

    double term1b = (1.0/5.0*b*b*b*b*b*mp*mq)-(1.0/5.0*a*a*a*a*a*mp*mq);
    double term2b = (1.0/4.0*b*b*b*b*mp*cq)-(0.25*a*a*a*a*mp*cq);
    double term3b = (0.25*b*b*b*b*mq*cp)-(0.25*a*a*a*a*mq*cp);
    double term4b = (1.0/3.0*b*b*b*cp*cq)-(1.0/3.0*a*a*a*cp*cq);

    double reta = term1a+term2a+term3a+term4a;
    double retb= term1b+term2b+term3b+term4b;

    return reta-retb;

}

double max(double a, double b){
    if (a > b) return a;
    return b;
}
double min(double a, double b){
    if (a < b) return a;
    return b;
}



double CalcIntegral_mass(IndexDimension p, IndexDimension q, Depth TP, Depth TQ){



    double mpLeft[DimensionSparseGrid];
    double mqLeft[DimensionSparseGrid];
    double mpRight[DimensionSparseGrid];
    double mqRight[DimensionSparseGrid];

    double cpLeft[DimensionSparseGrid];
    double cqLeft[DimensionSparseGrid];
    double cpRight[DimensionSparseGrid];
    double cqRight[DimensionSparseGrid];

    double hp[DimensionSparseGrid];
    double hq[DimensionSparseGrid];

    double coordp[DimensionSparseGrid];
    double coordq[DimensionSparseGrid];

    double suppleftP[DimensionSparseGrid];
    double supprightP[DimensionSparseGrid];

    double suppleftQ[DimensionSparseGrid];
    double supprightQ[DimensionSparseGrid];

    for(int d=0; d<DimensionSparseGrid; d++){

        mpLeft[d]=POW2(TP.at(d));


        hp[d]= 1.0/mpLeft[d];

        mpRight[d]= -1.0*mpLeft[d];
        coordp[d]=p.coordinate(d);
        cpLeft[d]=-(coordp[d]-hp[d])/hp[d];
        cpRight[d]=+(coordp[d]+hp[d])/hp[d];



        suppleftP[d]=coordp[d]-hp[d];
        supprightP[d]=coordp[d]+hp[d];






        mqLeft[d]=POW2(TQ.at(d));
        hq[d]= 1.0/mqLeft[d];
        mqRight[d]= -1.0*mqLeft[d];
        coordq[d]=q.coordinate(d);
        cqLeft[d]=-(coordq[d]-hq[d])/hq[d];
        cqRight[d]=+(coordq[d]+hq[d])/hq[d];


        suppleftQ[d]=coordq[d]-hq[d];
        supprightQ[d]=coordq[d]+hq[d];

    }

    double leftP_leftQ;
    double leftP_rightQ;
    double rightP_leftQ;
    double rightP_rightQ;


    double sum_inner=0.0;
    double prod=1.0;


    for (int dinner = 0; dinner < DimensionSparseGrid; dinner++) {


        bool docase1 = true;
        bool docase2 = true;
        bool docase3 = true;
        bool docase4 = true;
        if (q.getIndex(dinner) == 0) {
            docase1 = false;
            docase3 = false;
        }
        if (q.getIndex(dinner) == 1) {
            docase2 = false;
            docase4 = false;
        }

        sum_inner = 0.0;
        leftP_leftQ = 0.0;
        leftP_rightQ = 0.0;
        rightP_leftQ = 0.0;
        rightP_rightQ = 0.0;

        // int d/dx vp * d/dx vp mit x==dinner

        //1st case: integral over leftp cap leftq
        double intLp = suppleftP[dinner];
        double intRp = coordp[dinner];

        double intLq = suppleftQ[dinner];
        double intRq = coordq[dinner];


        if (docase1) {


            if (intRp < intLq) {
                leftP_leftQ = 0.0;
            } else if (intRq < intLp) {
                leftP_leftQ = 0.0;
            } else {
                leftP_leftQ = integral(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner], mqLeft[dinner],
                                       cpLeft[dinner], cqLeft[dinner]);
            }
        }





        //2nd case: integral over leftp cap rightq
        intLp = suppleftP[dinner];
        intRp = coordp[dinner];

        intLq = coordq[dinner];
        intRq = supprightQ[dinner];

        if (docase2) {
            if (intRp < intLq) {
                leftP_rightQ = 0.0;
            } else if (intRq < intLp) {
                leftP_rightQ = 0.0;
            } else {
                leftP_rightQ = integral(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner],
                                        mqRight[dinner], cpLeft[dinner], cqRight[dinner]);
            }
        }

        //3rd case: integral over rightp cap leftq
        intLp = coordp[dinner];
        intRp = supprightP[dinner];

        intLq = suppleftQ[dinner];
        intRq = coordq[dinner];

        if (docase3) {

            if (intRp < intLq) {
                rightP_leftQ = 0.0;
            } else if (intRq < intLp) {
                rightP_leftQ = 0.0;
            } else {
                rightP_leftQ = integral(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                        mqLeft[dinner], cpRight[dinner], cqLeft[dinner]);
            }
        }

        //4th case: integral over rightp cap rightq
        intLp = coordp[dinner];
        intRp = supprightP[dinner];

        intLq = coordq[dinner];
        intRq = supprightQ[dinner];

        if (docase4) {

            if (intRp < intLq) {
                rightP_rightQ = 0.0;
            } else if (intRq < intLp) {
                rightP_rightQ = 0.0;
            } else {
                rightP_rightQ = integral(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                         mqRight[dinner], cpRight[dinner], cqRight[dinner]);
            }
        }


        sum_inner = leftP_leftQ + leftP_rightQ + rightP_leftQ + rightP_rightQ;


        prod *= sum_inner;
    }





    return  prod;

}

double CalcIntegral_coeff(IndexDimension p, IndexDimension q, Depth TP, Depth TQ){



    double mpLeft[DimensionSparseGrid];
    double mqLeft[DimensionSparseGrid];
    double mpRight[DimensionSparseGrid];
    double mqRight[DimensionSparseGrid];

    double cpLeft[DimensionSparseGrid];
    double cqLeft[DimensionSparseGrid];
    double cpRight[DimensionSparseGrid];
    double cqRight[DimensionSparseGrid];

    double hp[DimensionSparseGrid];
    double hq[DimensionSparseGrid];

    double coordp[DimensionSparseGrid];
    double coordq[DimensionSparseGrid];

    double suppleftP[DimensionSparseGrid];
    double supprightP[DimensionSparseGrid];

    double suppleftQ[DimensionSparseGrid];
    double supprightQ[DimensionSparseGrid];

    for(int d=0; d<DimensionSparseGrid; d++){

        mpLeft[d]=POW2(TP.at(d));


        hp[d]= 1.0/mpLeft[d];

        mpRight[d]= -1.0*mpLeft[d];
        coordp[d]=p.coordinate(d);
        cpLeft[d]=-(coordp[d]-hp[d])/hp[d];
        cpRight[d]=+(coordp[d]+hp[d])/hp[d];



        suppleftP[d]=coordp[d]-hp[d];
        supprightP[d]=coordp[d]+hp[d];






        mqLeft[d]=POW2(TQ.at(d));
        hq[d]= 1.0/mqLeft[d];
        mqRight[d]= -1.0*mqLeft[d];
        coordq[d]=q.coordinate(d);
        cqLeft[d]=-(coordq[d]-hq[d])/hq[d];
        cqRight[d]=+(coordq[d]+hq[d])/hq[d];


        suppleftQ[d]=coordq[d]-hq[d];
        supprightQ[d]=coordq[d]+hq[d];

    }

    double leftP_leftQ;
    double leftP_rightQ;
    double rightP_leftQ;
    double rightP_rightQ;


    double sum_inner=0.0;
    double prod=1.0;


    for (int dinner = 0; dinner < DimensionSparseGrid; dinner++) {


        bool docase1 = true;
        bool docase2 = true;
        bool docase3 = true;
        bool docase4 = true;
        if (q.getIndex(dinner) == 0) {
            docase1 = false;
            docase3 = false;
        }
        if (q.getIndex(dinner) == 1) {
            docase2 = false;
            docase4 = false;
        }

        sum_inner = 0.0;
        leftP_leftQ = 0.0;
        leftP_rightQ = 0.0;
        rightP_leftQ = 0.0;
        rightP_rightQ = 0.0;

        // int d/dx vp * d/dx vp mit x==dinner

        //1st case: integral over leftp cap leftq
        double intLp = suppleftP[dinner];
        double intRp = coordp[dinner];

        double intLq = suppleftQ[dinner];
        double intRq = coordq[dinner];


        if (docase1) {


            if (intRp < intLq) {
                leftP_leftQ = 0.0;
            } else if (intRq < intLp) {
                leftP_leftQ = 0.0;
            } else {
                leftP_leftQ = integral_coeff(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner], mqLeft[dinner],
                                             cpLeft[dinner], cqLeft[dinner]);
            }
        }





        //2nd case: integral over leftp cap rightq
        intLp = suppleftP[dinner];
        intRp = coordp[dinner];

        intLq = coordq[dinner];
        intRq = supprightQ[dinner];

        if (docase2) {
            if (intRp < intLq) {
                leftP_rightQ = 0.0;
            } else if (intRq < intLp) {
                leftP_rightQ = 0.0;
            } else {
                leftP_rightQ = integral_coeff(max(intLq, intLp), min(intRp, intRq), mpLeft[dinner],
                                              mqRight[dinner], cpLeft[dinner], cqRight[dinner]);
            }
        }

        //3rd case: integral over rightp cap leftq
        intLp = coordp[dinner];
        intRp = supprightP[dinner];

        intLq = suppleftQ[dinner];
        intRq = coordq[dinner];

        if (docase3) {

            if (intRp < intLq) {
                rightP_leftQ = 0.0;
            } else if (intRq < intLp) {
                rightP_leftQ = 0.0;
            } else {
                rightP_leftQ = integral_coeff(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                              mqLeft[dinner], cpRight[dinner], cqLeft[dinner]);
            }
        }

        //4th case: integral over rightp cap rightq
        intLp = coordp[dinner];
        intRp = supprightP[dinner];

        intLq = coordq[dinner];
        intRq = supprightQ[dinner];

        if (docase4) {

            if (intRp < intLq) {
                rightP_rightQ = 0.0;
            } else if (intRq < intLp) {
                rightP_rightQ = 0.0;
            } else {
                rightP_rightQ = integral_coeff(max(intLq, intLp), min(intRp, intRq), mpRight[dinner],
                                               mqRight[dinner], cpRight[dinner], cqRight[dinner]);
            }
        }


        sum_inner = leftP_leftQ + leftP_rightQ + rightP_leftQ + rightP_rightQ;


        prod *= sum_inner;
    }





    return  prod;

}


void InHomoBoundaryRHSHierarchical::multiply(VectorSparseG &hierarchical, VectorSparseG &dirichlet) {

        DepthList depthList(*dirichlet.getSparseGrid());
    DepthList depthList1(*hierarchical.getSparseGrid());
        for(auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
            nodal = 0.0;
            Depth T = *it;



            auto iter = GetNextSmallerDepthsIterator(T);
            do {
                Depth Tlocal = *iter;
                if(Tlocal>>0) {

                    SingleDepthHashGrid &depthGrid = dirichlet.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                            Tlocal);
                    const auto &mapping = depthGrid._mapPosToGridPos;
                    // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
                    if (depthGrid.getNumberOfEntries() > 0) {



                        #pragma omp parallel for schedule(runtime)
                        for (size_t i = 0; i < mapping.size(); i++) {

                            IndexDimension INodal = depthGrid._map.getIndexOfTable(i);

                            double sum = 0.0;
                            for (auto it_inner = depthList1.begin_all();
                                 it_inner != depthList1.end_all(); ++it_inner) {

                                Depth THier = *it_inner;
                                if (!(THier >> 0)) {

                                    SingleDepthHashGrid &depthGrid_hier = hierarchical.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                                            THier);
                                    const auto &mapping_hier = depthGrid_hier._mapPosToGridPos;
                                    const auto &map_hier = depthGrid_hier._map;
                                    auto end = mapping_hier.size();


                                    if (depthGrid_hier.getNumberOfEntries() > 0) {

                                            for (size_t j = 0; j < end; j++) {
                                                IndexDimension IHier = map_hier.getIndexOfTable(j);


                                                double val = CalcIntegral(INodal, IHier, T, THier);


                                                sum += val * hierarchical.getValue(mapping_hier[j]);
                                            }



                                    }
                                }
                            }


                            nodal.setValue(mapping[i], sum);



                        }
                    }


                }
            }while(iter.next());

            ConvertToPrewavelet2(nodal, dirichlet, T);
            nodal = 0.0;
        }


}

void InHomoBoundaryRHSHierarchical::multiply_mass(VectorSparseG &hierarchical, VectorSparseG &dirichlet) {

    DepthList depthList(*dirichlet.getSparseGrid());
    DepthList depthList1(*hierarchical.getSparseGrid());
    for(auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
        nodal = 0.0;
        Depth T = *it;



        auto iter = GetNextSmallerDepthsIteratorInner(T);
        do {
            Depth Tlocal = *iter;
            if(Tlocal>>0) {

                SingleDepthHashGrid &depthGrid = dirichlet.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                        Tlocal);
                const auto &mapping = depthGrid._mapPosToGridPos;

                if (depthGrid.getNumberOfEntries() > 0) {
#pragma omp parallel for schedule(runtime)
                    for (size_t i = 0; i < mapping.size(); i++) {

                        IndexDimension INodal = depthGrid._map.getIndexOfTable(i);

                        double sum = 0.0;
                        for (auto it_inner = depthList1.begin_all();
                             it_inner != depthList1.end_all(); ++it_inner) {

                            Depth THier = *it_inner;
                            if (!(THier >> 0)) {

                                SingleDepthHashGrid &depthGrid_hier = hierarchical.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                                        THier);
                                const auto &mapping_hier = depthGrid_hier._mapPosToGridPos;
                                const auto &map_hier = depthGrid_hier._map;
                                auto end = mapping_hier.size();


                                if (depthGrid_hier.getNumberOfEntries() > 0) {
//#pragma omp parallel reduction(+:sum)
                                    {
//#pragma omp for schedule(runtime)
                                        for (size_t j = 0; j < end; j++) {
                                            IndexDimension IHier = map_hier.getIndexOfTable(j);


                                            double val = CalcIntegral_mass(INodal, IHier, T, THier);


                                            sum += val * hierarchical.getValue(mapping_hier[j]);
                                        }


                                    }
                                }
                            }
                        }




                        nodal.setValue(mapping[i], sum);



                    }
                }


            }
        }while(iter.next());

        ConvertToPrewavelet2(nodal, dirichlet, T);
        nodal = 0.0;
    }
}

void InHomoBoundaryRHSHierarchical::multiply_mass_coeff(VectorSparseG &hierarchical, VectorSparseG &dirichlet) {



    DepthList depthList(*dirichlet.getSparseGrid());
    DepthList depthList1(*hierarchical.getSparseGrid());

    for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
        Depth T = *it;
        auto iter = GetNextSmallerDepthsIterator(T);
        do {
            Depth Tlocal = *iter;
            if (Tlocal >> 0) {

                SingleDepthHashGrid &depthGrid = dirichlet.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                        Tlocal);
                const auto &mapping = depthGrid._mapPosToGridPos;

                if (depthGrid.getNumberOfEntries() > 0) {


#pragma omp parallel for schedule(runtime)
                        for (size_t i = 0; i < mapping.size(); i++) {

                            IndexDimension INodal = depthGrid._map.getIndexOfTable(i);

                            //double sum = 0.0;
                            double sum_omp = 0.0;
                            double sum = 0.0;


                            for (auto it_inner = depthList1.begin_all();
                                 it_inner != depthList1.end_all(); ++it_inner) {

                                Depth THier = *it_inner;
                                if (!(THier >> 0)) {

                                    SingleDepthHashGrid &depthGrid_hier = hierarchical.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                                            THier);
                                    const auto &mapping_hier = depthGrid_hier._mapPosToGridPos;
                                    const auto &map_hier = depthGrid_hier._map;
                                    auto end = mapping_hier.size();


                                        if (depthGrid_hier.getNumberOfEntries() > 0) {
//#pragma omp parallel reduction(+:sum_omp)
                                          {
//#pragma omp for schedule(runtime)
                                            for (size_t j = 0; j < end; j++) {
                                                //IndexDimension IHier = depthGrid_hier._map.getIndexOfTable(j);
                                                IndexDimension IHier = map_hier.getIndexOfTable(j);


                                                double val = CalcIntegral_coeff(INodal, IHier, T, THier);

                                                sum_omp += val *
                                                           hierarchical.getValue(mapping_hier[j]);


                                            }

                                        }
                                    }

                                }
                            }


                            nodal.setValue(mapping[i], sum_omp);


                    }
                }
            }


        } while (iter.next());

        ConvertToPrewavelet2(nodal, dirichlet, T);
        nodal = 0.0;
    }














  /*  DepthList depthList(*dirichlet.getSparseGrid());

    for(auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
        Depth T = *it;




        auto iter = GetNextSmallerDepthsIterator(T);
        do {
            Depth Tlocal = *iter;
            if(Tlocal>>0) {

                SingleDepthHashGrid &depthGrid = dirichlet.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                        Tlocal);
                const auto &mapping = depthGrid._mapPosToGridPos;

                if (depthGrid.getNumberOfEntries() > 0) {
                    for (size_t i = 0; i < mapping.size(); i++) {

                        IndexDimension INodal = depthGrid._map.getIndexOfTable(i);

                        double sum = 0.0;
                        for (unsigned long k = 0;
                             k < hierarchical.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {
                            IndexDimension IHier = hierarchical.getSparseGrid()->getIndexOfTable(k);
                            Depth THier(IHier);

                            if (!(THier >> 0)) {



                                double val = CalcIntegral_coeff(INodal, IHier, T, THier);




                                sum += val * hierarchical.getValue(k);


                            }
                        }



                        nodal.setValue(mapping[i], sum);



                    }
                }


            }
        }while(iter.next());

        ConvertToPrewavelet2(nodal, dirichlet, T);
        nodal = 0.0;
    }*/
}
