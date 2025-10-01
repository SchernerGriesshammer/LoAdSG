/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
#include "coordinate.h"

Coordinate::Coordinate(AdaptiveSparseGrid_Base& sg, int dimension_) {
  assertDimension(dimension_);
  dimension = dimension_;
  sparseGrid = &sg;
};

