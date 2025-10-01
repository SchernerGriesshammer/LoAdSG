//
// Created by scherner on 26.07.21.
//

#ifndef GRUN_TESTING_H
#define GRUN_TESTING_H

#include "../sgrid/sparseGrid.h"
#include "../fullgrid/FullGrid.h"
#include "../extemp/coordinate.h"
#include "old_versions/MatrixVectorMultiplicationPrewavelets.h"


#include "../MatrixVectorMultiplication/MatrixVectorInhomogen.h"

#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"


void checkgrid(AdaptiveSparseGrid &sgrid, MultiLevelAdaptiveSparseGrid &mgrid);

void checkgrid_dirichlet(AdaptiveSparseGrid &sgrid, MultiLevelAdaptiveSparseGrid &mgrid);


#endif //GRUN_TESTING_H
