#ifndef SINGLEDEPTHCELLSTRUCTURE_H
#define SINGLEDEPTHCELLSTRUCTURE_H


#include "../sgrid/sparseGrid.h"
#include "celldimension.h"

class SingleDepthHashCellStructure{

#ifndef BENCHMARKING
private:
#else
    public:
#endif




    AdaptiveSparseGrid_Base &_grid;
    const Depth _depth;


public:
    SimpleMultiHash<CellDimension,Depth> _map;

    SingleDepthHashCellStructure(Depth depth, AdaptiveSparseGrid_Base &grid);//:_grid(grid),_depth(depth){}


    bool addPoint(const CellDimension &key);


    Depth getDepth() const {return _depth;}

    bool isDepth(const Depth & depth) const {return _depth == depth;}

    size_t getNumberOfEntries() const { return _map.size();}

};

#endif