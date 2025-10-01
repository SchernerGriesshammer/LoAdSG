//
// Created by scherner on 27.07.21.
//

#ifndef GRUN_REFINEMENT_H
#define GRUN_REFINEMENT_H

#include "../BasisTransformations/BasisTransformations.h"


class Refinement{
public:
    Refinement(AdaptiveSparseGrid& grid):refined(grid){};

    bool apply(VectorSparseG&nodal,double eps);
private:

    VectorSparseG refined;
};

bool AdaptiveRefinement(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps);


bool AdaptiveRefinement(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps,int level);

bool AdaptiveRefinement(AdaptiveSparseGrid& oldgrid,AdaptiveSparseGrid& newgrid, VectorSparseG &nodal, double eps);


bool AdaptiveRefinementL2(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps);

bool AdaptiveRefinementEnergy(AdaptiveSparseGrid &sgrid, VectorSparseG &nodal, double eps);

#endif //GRUN_REFINEMENT_H
