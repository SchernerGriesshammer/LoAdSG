
#include "SingleDepthCellStructure.h"

SingleDepthHashCellStructure::SingleDepthHashCellStructure(Depth depth, AdaptiveSparseGrid_Base &grid): _grid(grid), _depth(depth), _map(depth){};

bool SingleDepthHashCellStructure::addPoint(const CellDimension &key) {
    if (!_map.addPoint(key)) return false;
    return true;
}
