#ifndef MULTIDEPTHHASHGRID_H
#define MULTIDEPTHHASHGRID_H

#include "depth.h"
#include "../primes/prime.h"

#include <vector>
#include <unordered_map>
#include "simpleMultiHash.h"
#include <algorithm>


class AdaptiveSparseGrid_Base;
class MultiDepthHashGrid;
class MultiLevelVector;


class SingleDepthHashGrid{

#ifndef BENCHMARKING
    private:
#else
public:
#endif
        friend MultiDepthHashGrid;
        friend MultiLevelVector;
        // std::unordered_multimap<unsigned long,unsigned long> _map;// = std::unordered_multimap<unsigned,unsigned>(100); // maps hash -> position in sparseGrid second table
        AdaptiveSparseGrid_Base &_grid;
        const Depth _depth;
        // std::unordered_multimap<size_t,size_t>::const_iterator _occupied(unsigned long &indexOfData, const IndexDimension &I);
    public:
        SimpleMultiHash<IndexDimension,Depth> _map;

        vector<unsigned long> _mapPosToGridPos;

        SingleDepthHashGrid(Depth depth, AdaptiveSparseGrid_Base &grid);//:_grid(grid),_depth(depth){}



        bool addPoint(const IndexDimension &key, unsigned long pos);

        // bool addPoint(const unsigned long hash, unsigned long pos);

        bool occupied(unsigned long &indexOfData, const IndexDimension &I) const;

        bool occupied(unsigned long &indexOfData, const IndexDimension &I, Depth& T) const;

        Depth getDepth() const {return _depth;}

        bool isDepth(const Depth & depth) const {return _depth == depth;}

        size_t getNumberOfEntries() const { return _map.size();}

        friend std::ostream& operator<< (std::ostream& stream, const SingleDepthHashGrid& grid){
            stream << "SingleDepthHashGrid: Depth(";
            for (size_t i = 0; i < DimensionSparseGrid-1; i++)
            {
                stream << grid._depth.at(i) << ", ";
            }
            stream << grid._depth.at(DimensionSparseGrid-1) << ") with " << grid.getNumberOfEntries() << " entries";
            return stream;
        }
};


class MultiDepthHashGrid{
#ifndef BENCHMARKING
    private:
#else
public:
#endif
        std::unordered_multimap<unsigned long,SingleDepthHashGrid> _map;// = std::unordered_multimap<unsigned,unsigned>(100); // maps hash -> position in sparseGrid second table
        AdaptiveSparseGrid_Base &_grid;
        // std::unordered_multimap<size_t,size_t>::const_iterator _occupied(unsigned long &indexOfData, const IndexDimension &I);
        friend AdaptiveSparseGrid_Base;
        /**
         * @brief Adds the point to the multiDepthgrid. This method should never be called as is!! Only AdaptiveSparseGrid_Base should call this method upon adding a new point to make sure both grids remain in sync.
         * That is why this method is private 
         * 
         * @param key 
         * @param pos 
         * @return true 
         * @return false 
         */

    public:
        bool addPoint(const IndexDimension &key, unsigned long pos);
        bool addPoint(const IndexDimension &key, unsigned long pos, Depth D);
    /**
     * @brief Reorders the corresponding sparse grid in order to improve cache hit performance.
     * WARNING: SparseVectors already created will not be reordered and will result in wrong calculation results.
     * 
     * @param collision If true will reorded by collision. If false will reorder by Depth
     * @return std::vector<unsigned long> the new order where [oldIdx] = newIdx
     */
        std::vector<unsigned long> reorder(bool collision = false);
        MultiDepthHashGrid(AdaptiveSparseGrid_Base &grid): _grid(grid){}
        /**
         * @brief Get the Grid of the given depth. Will create a new grid if it does not already exist.
         * 
         * @param D 
         * @return SingleDepthHashGrid& 
         */
        SingleDepthHashGrid& getGridForDepth(const Depth &D);
        /**
         * @brief Tries to get the grid of the given depth. Will return nullptr if grid could not be found
         * 
         * @param D 
         * @return nullptr when nothing is found
         */
        SingleDepthHashGrid* tryGetGridForDepth(const Depth &D);
        /**
         * @brief Get all the Grids of the given depth that have the same Depth as D except for the dth element.
         * And where the dth element is smaller or equal the given d
         * 
         * @param D 
         * @return SingleDepthHashGrid& 
         */
        vector<SingleDepthHashGrid*> getGridsForDepthInDirection(const Depth &D, int d);
        /**
         * @brief Wrapper for getGridForDepth(const Depth &d)
         * 
         * @param I 
         * @return SingleDepthHashGrid& 
         */
        SingleDepthHashGrid& getGridForDepth(const IndexDimension &I);
        // bool addPoint( unsigned long hash, unsigned long pos);
        // bool occupied(unsigned long &indexOfData, const IndexDimension &I) const;
        unsigned long hash(Depth D);

        vector<SingleDepthHashGrid*> getAllGrids();
        // unsigned long hash(Depth D){
        //     unsigned long value = D.at(0);
        //     for (int d = 1; d < DimensionSparseGrid; ++d) {
        //         value = value + D.at(d) * PrimeNumbers::getPrimeForHash(d);
        //     }
        //     return value;
        // }
        void printSizes(){
            for (auto it = _map.begin(); it != _map.end(); ++it) {
                SingleDepthHashGrid& res = it->second;
                // res.getDepth().Print();
                // cout << "with " << res.getNumberOfEntries() << " entries."<<endl;
                cout << res << endl;
            }
        }
};
/*
bool MultiDepthHashGrid::addPoint(const IndexDimension &key, unsigned long pos){
    Depth D(key);
    SingleDepthHashGrid& sgrid = getGridForDepth(D);
    return sgrid.addPoint(key,pos);;
}

SingleDepthHashGrid& MultiDepthHashGrid::getGridForDepth(const Depth &D){
    unsigned long hashVal = hash(D);

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
*/
#endif  // MULTIDEPTHHASHGRID_H