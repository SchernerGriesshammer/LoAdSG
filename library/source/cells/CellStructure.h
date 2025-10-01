#ifndef CELLSTRUCTURE_H
#define CELLSTRUCTURE_H

#include "SingleDepthCellStructure.h"
#include "../iterator/depthIterator.h"
#include "iterator_neu.h"

class CellData{


private:

    std::unordered_multimap<unsigned long,SingleDepthHashCellStructure> _map;
    AdaptiveSparseGrid_Base& gridBase;





public:
    bool addPoint(const CellDimension &key);

    CellData(AdaptiveSparseGrid_Base &grid):gridBase(grid){

        CellDimension bottomLeft, topRight;


        DepthList depthList(grid);

        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
            Depth T = *it;

            for(int d=0; d<DimensionSparseGrid;d++){
                bottomLeft.replace(d, Static_CellOneD::nextRight(0,T.at(d)+1));
                topRight.replace(d,Static_CellOneD::nextLeft(1,T.at(d)+1));
            }
            for(RectangularIteratorCells iteratorCells(bottomLeft,topRight); iteratorCells.goon(); ++iteratorCells){
                CellDimension cellDimension = iteratorCells.getCell();
                if(cellingrid(cellDimension,grid))
                    addPoint(cellDimension);
            }





        }

    }

    bool cellingrid(CellDimension cell, AdaptiveSparseGrid_Base& grid){


        CellIndexIterator cellIndexIterator(&cell);
        for(;cellIndexIterator.goon(); ++cellIndexIterator){
            IndexDimension Index = cellIndexIterator.getIndex();
            unsigned long k;
            if(Index.isNotAtBoundary() && !(grid.occupied(k,Index))) return false;
        }
        return  true;
    }
    std::unordered_multimap<unsigned long,SingleDepthHashCellStructure>* getMap(){return &_map;};

    SingleDepthHashCellStructure& getGridForDepth(const Depth &D);
    vector<SingleDepthHashCellStructure*> getAllCells();

    int countCells();


    unsigned long hash(Depth D);


};






#endif
