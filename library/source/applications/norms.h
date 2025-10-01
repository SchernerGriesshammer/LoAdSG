//
// Created by to35jepo on 5/11/23.
//

#ifndef RUN_NORMS_H
#define RUN_NORMS_H


#include "../extemp/extempAlg.h"

template <class A>
double L_infty ( const ExprSparseG<A>& a ) {
    const A& ao ( a );

    AdaptiveSparseGrid_Base* sparseGrid = ao.getSparseGrid();

    unsigned long endIndex   = ao.length;
    //double* dataTable        = sparseGrid->dataTable;
    dataInteger* secondTable = sparseGrid->secondTable;
    //unsigned int numberOfData = sparseGrid->numberOfData;


    double maximum = 0.0;
    if(ao.getDescription().isIndexNeeded()) {
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(secondTable[i]!=0) {

                //double x = ao.getValue(&dataTable[i * numberOfData],sparseGrid->getIndexOfTable(i));
                double x = ao.getValue(i,sparseGrid->getIndexOfTable(i));

                if(x < 0.0) x = x * (-1.0);
                if(x>maximum) maximum = x;
            }
        }
        return maximum;
    }
    else {
        IndexDimension Idummy;
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(secondTable[i]!=0) {


                //int size1 = *(&dataTable[i*numberOfData] + 1) - dataTable[i*numberOfData];

                //cout << "number " << ao.getNumber() << endl;
                //cout <<"size of dataTable[i*numberOfData]   " << size1 << endl;
                //cout <<"size of dataTable[i*numberOfData]   " << dataTable[i*numberOfData][1]<< endl;
                //double x = ao.getValue(&dataTable[i * numberOfData],Idummy);
                double x = ao.getValue(i,Idummy);

                if(x < 0.0) x = x * (-1.0);
                if(x>maximum) maximum = x;
            }
        }
        return maximum;
    }
}


template <class A>
double L_2 ( const ExprSparseG<A>& a ) {
    const A& ao ( a );

    AdaptiveSparseGrid_Base* sparseGrid = ao.getSparseGrid();

    unsigned long endIndex   = sparseGrid->getMaximalOccupiedSecondTable();




    double sum = 0.0;
    int dofs = 0;
    if(ao.getDescription().isIndexNeeded()) {
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(sparseGrid->getActiveTable()[i]) {


                //double x = ao.getValue(&dataTable[i * numberOfData],sparseGrid->getIndexOfTable(i));
                double x = ao.getValue(i, sparseGrid->getIndexOfTable(i));

                dofs++;
                sum += (x * x);
            }

        }
        return ::sqrt(sum/dofs);
    }
    else {
        IndexDimension Idummy;
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(sparseGrid->getActiveTable()[i]) {

                double x = ao.getValue(i, Idummy);

                dofs++;
                sum += (x * x);
            }

        }
        return ::sqrt(sum/dofs);
    }
}


#endif //RUN_NORMS_H
