//
// Created by to35jepo on 7/30/24.
//

#ifndef RUN_H1_NORM_H
#define RUN_H1_NORM_H


#include "../sgrid/sparseGrid.h"
#include "../sgrid/multilevelSparseGrid.h"

class H1{
public:
    H1(AdaptiveSparseGrid& grid,MultiLevelAdaptiveSparseGrid& mgrid_):sparseGrid(grid),mgrid(mgrid_){};


    double getValue(VectorSparseG& u);


private:

   AdaptiveSparseGrid& sparseGrid;
   MultiLevelAdaptiveSparseGrid& mgrid;


};

class L2{
public:
    L2(AdaptiveSparseGrid& grid,MultiLevelAdaptiveSparseGrid& mgrid_):sparseGrid(grid),mgrid(mgrid_){};


    double getValue(VectorSparseG& u);


private:

    AdaptiveSparseGrid& sparseGrid;
    MultiLevelAdaptiveSparseGrid& mgrid;


};



#endif //RUN_H1_NORM_H
