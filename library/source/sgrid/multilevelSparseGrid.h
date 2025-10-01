//
// Created by scherner on 20.05.21.
//

#ifndef GRUN_MULTILEVELSPARSEGRID_H
#define GRUN_MULTILEVELSPARSEGRID_H

#include "sparseGrid.h"
#include <algorithm>
#include <list>


class MultiLevelAdaptiveSparseGrid{
public:


    MultiLevelAdaptiveSparseGrid(AdaptiveSparseGrid *grid){


        primeTableLength = PrimeNumbers::getNextPrime(10000000);
        secondTableLength = 10000000;



        primeTable = new dataInteger[primeTableLength];
        for(unsigned long i=0;i<primeTableLength;++i) primeTable[i] = 0;

        secondTable = new dataInteger[secondTableLength];
        for(unsigned long i=0;i<secondTableLength;++i) secondTable[i] = 0;

        isActiveNodeTable= new bool[secondTableLength];
        for (unsigned long i = 0; i < secondTableLength; ++i) isActiveNodeTable[i] = false;


        indicesSecondTable = new indexInteger[secondTableLength *
                                              DimensionSparseGrid];///>for every point we need to store three indices i.e. i_x and i_y
        for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSecondTable[i] = 0;

        minimalEmptySecondTable = 0;
        maximalOccupiedSecondTable = 0;


        depthTable = new int[secondTableLength * DimensionSparseGrid];
        for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) depthTable[i] = -1;
        createbysparsegrid(grid);
    }

    ~MultiLevelAdaptiveSparseGrid() {
        delete[] primeTable;
        delete[] secondTable;
        delete[] isActiveNodeTable;

        delete[] indicesSecondTable;

        delete[] depthTable;
    };

    void createbysparsegrid(AdaptiveSparseGrid_Base *grid);


    void PrintActiveHangingDepth(int level, Depth T);

    void setDepthInTable(unsigned long k, Depth T) {

        for (int d = 0; d < DimensionSparseGrid; d++)
            depthTable[d + k * DimensionSparseGrid] = T.at(d);

    };

    inline bool occupied(unsigned long &indexOfData, IndexDimension &Index, Depth &T) {


        unsigned long indArray = hash(Index, T);

        unsigned long i = primeTable[indArray];

        if (i == 0) { // in prime table is something free;
            return false;
        } else { // anaylse what is going on
            i = i - 1; // shift since data are stored with shift


            if (Index == getIndexOfTable(i) && T == getDepthOfTable(i)) {
                indexOfData = i;
                return true;
            } else { // search in second table;
                unsigned long iNext = secondTable[i];
                while (iNext > 1) {
                    i = iNext - 2;
                    if (Index == getIndexOfTable(i) && T == getDepthOfTable(i)) {
                        indexOfData = i;
                        return true;
                    }
                    iNext = secondTable[i];
                }
                return false;
            }
        }

    };


    inline unsigned long hash(IndexDimension &Index, Depth &Tfine) {
        unsigned long value = Index.getIndex(0);
        for (int d = 1; d < DimensionSparseGrid; ++d) {
            value = value + Index.getIndex(d) * PrimeNumbers::getPrimeForHash(d);
        }

        for (int d = 0; d < DimensionSparseGrid; d++) {
            value = value + Tfine.at(d) * PrimeNumbers::getPrimeForHash(d + DimensionSparseGrid);
        }
        return value % primeTableLength;
    };

    inline Depth getDepthOfTable(unsigned long i);

    unsigned long getMaximalOccupiedSecondTable() { return maximalOccupiedSecondTable;};
    unsigned long getLengthSecondTable() { return secondTableLength; }

    int getDOFS();
private:
    unsigned long primeTableLength;
    unsigned long secondTableLength;
    unsigned long minimalEmptySecondTable;
    unsigned long maximalOccupiedSecondTable;
    dataInteger *primeTable;      ///> 0 means empty; v>0 means v-1 is array index
    dataInteger *secondTable;     ///> 0 means empty; 1 means occupied, but no next;  v>1 means v-2 is next array index
    bool *isActiveNodeTable;///> true means AdaptiveSparseGrid_Base::getIndexOfTable (i) is active node.

    indexInteger *indicesSecondTable; ///> Contains coded indices. Length: secondTableLength*DimensionSparseGrid
    int *depthTable;
    std::list<Depth> liste;

    bool AddPointDepth(IndexDimension &Index, Depth &Tfine);

    void setIndexAndDepthInTable(const IndexDimension Index, Depth Tfine, unsigned long iSetz);

    inline IndexDimension getIndexOfTable(unsigned long i);


    unsigned long getFreeSpaceNumberInSecondTable();

};


Depth MultiLevelAdaptiveSparseGrid::getDepthOfTable(unsigned long i) {
    Depth T(0);
    for (int d = 0; d < DimensionSparseGrid; d++) {
        T.set(depthTable[d + i * DimensionSparseGrid], d);
    }
    return T;
}


inline IndexDimension MultiLevelAdaptiveSparseGrid::getIndexOfTable(unsigned long i) {
    IndexDimension back;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        back.replace(d, indicesSecondTable[d + i * DimensionSparseGrid]);
    }
    return back;
}
#endif //GRUN_MULTILEVELSPARSEGRID_H
