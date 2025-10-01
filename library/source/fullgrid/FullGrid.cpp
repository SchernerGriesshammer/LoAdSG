//
// Created by scherner on 26.07.21.
//

#include "FullGrid.h"
#include "../iterator/RectangularIterator.h"

void fullgrid(int depth, AdaptiveSparseGrid &adaptiveSGrid) {
    Depth T(depth);

    IndexDimension minI;
    IndexDimension maxI;

    for (int i = 0; i < DimensionSparseGrid; i++) {
        minI.replace(i, 0);
        maxI.replace(i, 1);

    }

   for (int i = 0; i < DimensionSparseGrid; i++) {
        minI = minI.nextRight(i, depth);
        maxI = maxI.nextLeft(i, depth);

    }


    adaptiveSGrid.AddPoint(minI);
    for (RectangularIterator iter(minI, maxI, T); iter.goon(); ++iter) {
        IndexDimension Index = iter.getIndex();
        adaptiveSGrid.AddPoint(Index);
    }

};

