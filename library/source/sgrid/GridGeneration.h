//
// Created by to35jepo on 4/12/23.
//

#ifndef SGRUN_GRIDGENERATION_H
#define SGRUN_GRIDGENERATION_H

#include "../indices/index.h"
#include "../iterator/depthIterator.h"
#include "sparseGrid.h"

void support(const IndexDimension& P, IndexDimension& Imin, IndexDimension& Imax);
bool overlap(const IndexDimension& P, const IndexDimension& Q, IndexDimension& Imin, IndexDimension& Imax);
Depth max(Depth& TP, Depth& TQ);

bool refill(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q);
bool refill_NEU(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q, std::vector<IndexDimension>& pointsToAdd);
void addPoints(AdaptiveSparseGrid& dGrid);
void addPoints2(AdaptiveSparseGrid& dGrid);


void addPoints2_NEU(AdaptiveSparseGrid& dGrid, DepthList list);
#endif //SGRUN_GRIDGENERATION_H
