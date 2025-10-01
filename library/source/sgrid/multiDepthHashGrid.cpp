
#include <vector>
#include <unordered_map>
#include "multiDepthHashGrid.h"
#include "sparseGrid.h"

unsigned long MultiDepthHashGrid::hash(Depth D){
    unsigned long value = D.at(0);
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        value = value + D.at(d) * PrimeNumbers::getPrimeForHash(d);
    }
    return value;
}

    unsigned long hashOriginal(const IndexDimension& index,const unsigned long tableSize)  {

        Depth T(index);
        unsigned long value = index.getIndex(0);
        for (int d = 1; d < DimensionSparseGrid; ++d) {
            value = value + index.getIndex(d) * PrimeNumbers::getPrimeForHash(d);
        }

        value = value + T.at(0) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid);
        for (int d = 1; d < DimensionSparseGrid; d++)
            value = value + T.at(d) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid + d);

        return value % tableSize;

    }





SingleDepthHashGrid::SingleDepthHashGrid(Depth depth, AdaptiveSparseGrid_Base &grid): _grid(grid), _depth(depth), _map(depth){};

bool MultiDepthHashGrid::addPoint(const IndexDimension &key, unsigned long pos){
    Depth D(key);
    SingleDepthHashGrid& sgrid = getGridForDepth(D);
    return sgrid.addPoint(key,pos);;
}

bool MultiDepthHashGrid::addPoint(const IndexDimension &key, unsigned long pos, Depth D){
    SingleDepthHashGrid& sgrid = getGridForDepth(D);
    return sgrid.addPoint(key,pos);;
}


std::vector<unsigned long> MultiDepthHashGrid::reorder(bool collision){
    size_t size = _grid.getMaximalOccupiedSecondTable();
    std::vector<unsigned long> resultFromTo(size);
    if(!collision){
    size_t nextPos = 0; 
    auto iter = _map.begin();
    while (iter!=_map.end())
    {
        SingleDepthHashGrid &singleGrid = iter->second;
        std::vector<unsigned long>& mapPosToGridPos =singleGrid._mapPosToGridPos;
        for (size_t i =0; i < mapPosToGridPos.size(); i++)
        {
            resultFromTo[mapPosToGridPos[i]] = nextPos;
            mapPosToGridPos[i] = nextPos;
            nextPos++;
        }
        
        

        ++iter;
    }
    }
    else
    { 
        std::vector<unsigned long> sortReference(size);
        for (size_t i = 0; i < size; i++)
        {
            sortReference[i]=i; 
        }
        std::sort(sortReference.begin(), sortReference.end(), [&](unsigned long i,unsigned long j) { return (Depth(_grid.getIndexOfTable(i))<Depth(_grid.getIndexOfTable(j))); });
        resultFromTo = sortReference;

        auto iter = _map.begin();
        while (iter!=_map.end())
        {
            SingleDepthHashGrid &singleGrid = iter->second;
            std::vector<unsigned long>& mapPosToGridPos =singleGrid._mapPosToGridPos;
            for (size_t i =0; i < mapPosToGridPos.size(); i++)
            {
                auto nextPos = resultFromTo[mapPosToGridPos[i]];
                mapPosToGridPos[i] = nextPos;
                nextPos++;
            }
            ++iter;
        }
    }

    
    //  iter = _map.begin();
    // while (iter!=_map.end())
    // {
    //     SingleDepthHashGrid &singleGrid = iter->second;
    //     auto singleMapIter = singleGrid._map.begin();
    //     while(singleMapIter != singleGrid._map.end()){
    //         cout<< singleMapIter->second <<endl;
    //         ++singleMapIter;
    //     }
    //     ++iter;
    // }
    // cout << "Reorder: " <<endl;
    // for (size_t i = 0; i <size; i++)
    // {
    //     cout << resultFromTo[i] << ", ";
    // }
    // cout << endl;
    //     cout << endl;

    unsigned long * gridTable = _grid.getSecondTable();
    unsigned long * gridTableTmp = new unsigned long[size]; 
    indexInteger* indicesTable =_grid.indicesSecondTable;
    indexInteger * indicesTableTmp = new indexInteger[size*DimensionSparseGrid]; 
    for (size_t i = 0; i <size; i++)
    {
        gridTableTmp[i] = gridTable[i]; //TODO memcpy
    }
    for (size_t i = 0; i <size*DimensionSparseGrid; i++)  
    {
        indicesTableTmp[i] = indicesTable[i];  //TODO memcpy
    }
    for (size_t i = 0; i <size; i++)
    {
        size_t oldPos = i;
        size_t newPos = resultFromTo[oldPos];
        // gridTable[newPos] = gridTableTmp[oldPos]; 
        unsigned long newVal = gridTableTmp[oldPos];
        if (newVal > 1)
        {
            newVal = resultFromTo[newVal-2]+2;
        }
        
        gridTable[newPos] = newVal;
        for (size_t d = 0; d < DimensionSparseGrid; d++)
        {
            indicesTable[newPos*DimensionSparseGrid+d] = indicesTableTmp[oldPos*DimensionSparseGrid+d];
        }
        
    }
    unsigned long*  primeTable =_grid.getPrimeTable();
    // cout << "Original: " <<endl;
    // for (size_t i = 0; i <4; i++)
    // {
    //     cout << primeTable[i] << ", ";
    // }
    //     cout << endl;
    // for (size_t i = 0; i <size; i++)
    // {
    //     cout << gridTableTmp[i] << ", ";
    // }
    //     cout << endl;
    //     cout << endl;
    for (size_t i = 0; i < _grid.getPrimeTableLength(); i++)
    {
        if(primeTable[i]==0) continue;
        primeTable[i] = resultFromTo[primeTable[i]-1]+1;
    }
    // cout << "Reordered: " <<endl;
    // for (size_t i = 0; i <4; i++)
    // {
    //     cout << primeTable[i] << ", ";
    // }
    //     cout << endl;
    // for (size_t i = 0; i <size; i++)
    // {
    //     cout << gridTable[i] << ", ";
    // }
    //     cout << endl;
    delete[] gridTableTmp;
    delete[] indicesTableTmp;
    return resultFromTo;
    
}


SingleDepthHashGrid* MultiDepthHashGrid::tryGetGridForDepth(const Depth &D){
    unsigned long hashVal = hash(D);

    auto its = _map.equal_range(hashVal);

    for (auto it = its.first; it != its.second; ++it) {
        if(it->second.isDepth(D)){
             return &it->second;
        }
    }
    return nullptr;
}

SingleDepthHashGrid& MultiDepthHashGrid::getGridForDepth(const Depth &D){
    unsigned long hashVal = hash(D);
    _map.hash_function();

    auto its = _map.equal_range(hashVal);

    for (auto it = its.first; it != its.second; ++it) {
        if(it->second.isDepth(D)){
             return it->second;
        }
    }
    // SingleDepthHashGrid* grid = new SingleDepthHashGrid(D,this->_grid);
    
    SingleDepthHashGrid grid(D,this->_grid);
    auto iter = _map.insert({hashVal,grid});
    
    return iter->second;
}

vector<SingleDepthHashGrid*>  MultiDepthHashGrid::getGridsForDepthInDirection(const Depth &D, int d){
    Depth cD = D;

    int maxd = cD.at(d);

    vector<SingleDepthHashGrid*> A;
    A.reserve(maxd);
    for (size_t i = 0; i < maxd; i++)
    {
        cD.set(i,d);

        SingleDepthHashGrid* grid = &getGridForDepth(cD);
        //SingleDepthHashGrid* grid =tryGetGridForDepth(cD);
        A.push_back(grid);
    }
    return A;

}
vector<SingleDepthHashGrid*> MultiDepthHashGrid::getAllGrids(){
    vector<SingleDepthHashGrid*> grids;
    grids.reserve(_map.size());
    auto iter = _map.begin();
    while (iter!=_map.end()){
        SingleDepthHashGrid* grid = &(iter->second);
        grids.push_back(grid);
        iter++;
    }
    return grids;
}


bool SingleDepthHashGrid::addPoint(const IndexDimension &key, unsigned long pos){
    if(!_map.addPoint(key)) return false;
    _mapPosToGridPos.push_back(pos);
    return true;
    // unsigned long hashVal = _grid.hash(key);
    // // _map.insert(std::unordered_multimap<unsigned,unsigned>::value_type(hashVal,currIndx)); 

    // unsigned long _tmpIndex;
    // if(! occupied(_tmpIndex,key)){
    //     _map.insert({hashVal,pos}); 
    // }
    
    // return true;
}

// std::unordered_multimap<size_t,size_t>::const_iterator SingleLevelHashGrid::_occupied( unsigned long &indexOfData, const IndexDimension &I) const{
//     size_t hash = _grid.hash(I);
//     auto its = _map.equal_range(hash);
//     for (auto it = its.first; it != its.second; ++it) {
//         if(I == _grid.getIndexOfTable(it->second)){
//              indexOfData = it->second;
//              return it;
//         }
//     }
//     return _map.end();
// }


bool SingleDepthHashGrid::occupied( unsigned long &indexOfData, const IndexDimension &I) const{
    bool result  = _map.occupied(indexOfData,I);
    if(!result) return false;
    indexOfData = _mapPosToGridPos[indexOfData];
    return result;
    // unsigned long hash = _grid.hash(I);
    // auto its = _map.equal_range(hash);
    // for (auto it = its.first; it != its.second; ++it) {
    //     if(I == _grid.getIndexOfTable(it->second)){
    //          indexOfData = it->second;
    //          return true;
    //     }
    // }
    // return false;
}

bool SingleDepthHashGrid::occupied(unsigned long &indexOfData, const IndexDimension &I, Depth &T) const {
    return false;
}
