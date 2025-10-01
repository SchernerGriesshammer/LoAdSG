#include "CellStructure.h"

bool CellData::addPoint(const CellDimension &key){
    Depth D;
    for(int d=0; d<DimensionSparseGrid; d++)
        D.set(key.getDepth(d),d);

    SingleDepthHashCellStructure& sgrid = getGridForDepth(D);
    return sgrid.addPoint(key);;
}

SingleDepthHashCellStructure &CellData::getGridForDepth(const Depth &D) {
    unsigned long hashVal = hash(D);
    _map.hash_function();

    auto its = _map.equal_range(hashVal);

    for (auto it = its.first; it != its.second; ++it) {
        if(it->second.isDepth(D)){
            return it->second;
        }
    }
    // SingleDepthHashGrid* grid = new SingleDepthHashGrid(D,this->_grid);

    SingleDepthHashCellStructure grid(D,this->gridBase);
    auto iter = _map.insert({hashVal,grid});

    return iter->second;
}

unsigned long CellData::hash(Depth D) {
    unsigned long value = D.at(0);
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        value = value + D.at(d) * PrimeNumbers::getPrimeForHash(d);
    }
    return value;
}


vector<SingleDepthHashCellStructure*> CellData::getAllCells(){
    vector<SingleDepthHashCellStructure*> grids;
    grids.reserve(_map.size());
    auto iter = _map.begin();
    while (iter!=_map.end()){
        SingleDepthHashCellStructure* grid = &(iter->second);
        if(grid->_map.size()>0)  grids.push_back(grid);
        iter++;
    }
    return grids;
}

int CellData::countCells() {

   int count=0;
    auto iter = _map.begin();
    while (iter!=_map.end()){
        SingleDepthHashCellStructure* grid = &(iter->second);
        count += int(grid->_map.size());
        iter++;
    }
    return count;
}
