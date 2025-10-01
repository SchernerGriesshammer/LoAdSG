//
// Created by to35jepo on 9/19/23.
//

#ifndef RUN_CONDITION_H
#define RUN_CONDITION_H


#include "../sgrid/sparseGrid.h"
#include "../sgrid/multilevelSparseGrid.h"

double potenzmethode(int k, AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, bool conditioned);
double potenzmethode_min(int k,AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid,bool conditioned, double lambda);

/**
 *
 * @param k Number of Iterations
 * @param grid Sparse Grid
 * @param mgrid MultiLevelGrid
 * @param emin Eigenvector
 * @return Eigenvalue
 */
double potenzmethode_invers(int k, AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, VectorSparseG& emin);

double condition(int k, AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid&mgrid, bool conditioned);

#endif //RUN_CONDITION_H
