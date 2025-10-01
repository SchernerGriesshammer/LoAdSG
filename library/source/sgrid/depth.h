/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
#ifndef DEPTH_H
#define DEPTH_H

#include "../myAssert.h"


#include <iostream>
#include <fstream>
#include <list>

#include "../indices/index.h"




/***
 * Einfache KLasse fuer Tiefe 
 ***/
class Depth {
public:


    inline Depth(){
        for (int d = 0; d < DimensionSparseGrid; d++)
            depth[d] = 0;
    };

    inline explicit Depth(int depthall) {

        for (int d = 0; d < DimensionSparseGrid; d++)
            depth[d] = depthall;
    };

    inline explicit Depth(IndexDimension indexD);

    inline Depth(const Depth &T_) {
        for (int d = 0; d < DimensionSparseGrid; d++) depth[d] = T_.depth[d];
    };

    inline int *returnTiefen() {
        return depth;
    };

    /**
     * @param t Tiefe
     * @param d Dimension 
     **/
    inline void set(int t, int d);

    inline int at(int d) const;


    inline unsigned int maxNorm();

    inline unsigned int LoneNorm();

    inline bool operator==(const Depth &TR) const;

    inline bool operator<=(const Depth &TR);

    inline bool operator<(const Depth &TR);

    inline bool operator>(const Depth &TR);

    inline bool operator>(int t);

    inline bool operator>>(int t);

    inline bool operator==(int t);

    inline void operator++(){
        for(int d=0;d<DimensionSparseGrid; d++){
            depth[d]++;
        }
    }


    inline void Print();

private:
    int depth[DimensionSparseGrid];
};







///////////////////////////////////////////////////////////////
// implementation of inline functions
///////////////////////////////////////////////////////////////


inline Depth::Depth(IndexDimension indexD) {
    for (int d = 0; d < DimensionSparseGrid; ++d)
        depth[d] = indexD.getDepth(d);
}

inline int Depth::at(int d) const{
    assertDimension(d);
    return depth[d];
}


inline void Depth::set(int t, int d) {
    assertDimension(d);
    depth[d] = t;
}

inline bool Depth::operator==(const Depth &TR) const{
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] != TR.depth[d]) return false;
    }
    return true;
}

inline bool Depth::operator<=(const Depth &TR) {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] > TR.depth[d]) return false;
    }
    return true;
}

inline bool Depth::operator<(const Depth &TR) {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] >= TR.depth[d]) return false;
    }
    return true;
}

inline bool Depth::operator>(const Depth &TR) {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] <= TR.depth[d]) return false;
    }
    return true;
}

inline bool Depth::operator>(int t) {


    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] <= t) return false;
    }
    return true;
}

inline bool Depth::operator>>(int t) {

    bool checkall = true;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (depth[d] <= t) checkall = false;

    }
    return checkall;
}

inline bool Depth::operator==(int t) {
    for (int d = 0; d < DimensionSparseGrid; d++) {
        if (depth[d] != t)return false;
    }
    return true;
}

inline unsigned int Depth::maxNorm() {
    unsigned int maxN = depth[0];
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        if (maxN < depth[d]) maxN = depth[d];
    }
    return maxN;
}

inline unsigned int Depth::LoneNorm() {
  unsigned int loneN = depth[0];
  for(int d=1;d<DimensionSparseGrid;++d) {
      loneN = loneN + depth[d];
  }
  return loneN;
}

inline void Depth::Print() {
    std::cout << "Depth ";
    std::cout << " (" << depth[0];
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        std::cout << ", " << depth[d];
    }
    std::cout << ") " ;

}

////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////


    
#endif  // DEPTH_H

