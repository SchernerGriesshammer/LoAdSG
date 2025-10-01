#ifndef CELLDIMENSION_H
#define CELLDIMENSION_H

#include "../myAssert.h"
#include "../indices/index.h"
#include "MultiDepthHashCellStructure.h"

#include <iostream>


class Static_CellOneD {

public:
    inline static double coordinate(indexInteger index);


    inline static indexInteger nextRight(indexInteger index, int level);

    inline static indexInteger nextLeft(indexInteger index, int level);

protected:


    inline static int depth(indexInteger index);

    inline static indexInteger father(const indexInteger index);

    inline static indexInteger leftSon(const indexInteger index);

    inline static indexInteger rightSon(const indexInteger index);

    inline static indexInteger nextLeft(indexInteger index);   ///> next left index on same level
    inline static indexInteger nextRight(indexInteger index); ///> next right index on same level

    ///> next left  index on level >= depth
    ///> next right index on level >= depth



private:
    inline static int internDepth(indexInteger index);
};


class CellOneD : private Static_CellOneD {
public:




    inline CellOneD(indexInteger ind) : index(ind) {};

    inline int depth();


    //! next left index on same level
    /*!
     * Example:
     * * IndexOneD::index = 2;
     * * IndexOneD::nextLeft()=0
     *
     */
    inline CellOneD nextLeft();
    //! next right index on same level

    inline CellOneD nextRight();
    //!  next left index on level \f$\geq\f$ depth
    /*!
     * Example:
     * * IndexOneD::index=2;
     * * IndexOneD::nextLeft(2)=4;
     */

    inline CellOneD nextRightLevel(int level);


    inline double coordinate();

    indexInteger getInteger() const { return index; }



private:
    indexInteger index;

};

//---------
#if __cplusplus >= 202002L
inline int Static_CellOneD::depth(indexInteger index) {
  return std::bit_width(index>>1);
}
#else

// static int Static_IndexOneD::depth_orignal(indexInteger index) {
//   int de = 0;
//   while(index!=1) {
//      de = de + 1;
//      index = index >> 1;
//   }
//   return de;
// }

inline int Static_CellOneD::depth(indexInteger index) {
    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_clz only works for ints or compile with cpp20 support");
    if(index==0) return 0;
    return __builtin_clz(index) ^ 31;
}
#endif

inline int CellOneD::depth() {
    return Static_CellOneD::depth(index);
}


double Static_CellOneD::coordinate(indexInteger index) {
/*
    double co = 0.0;
    while(index>0) {
        co = co * 0.5;
        if((index&1)==1) co = co + 1.0;
        else           co = co - 1.0;
        index = index>>1;
//cout << " co: " << co << endl;
    }
    return co;
*/

    if(index==0)return 0;
    int d = depth(index);
    //cout << d << endl;

    double h = 1.0/double(POW2(d));
    //cout << h << endl;
    indexInteger n= subtractFirstBit(index);
    //cout << n << endl;
    return  double(2*n+1)*h;


}
double CellOneD::coordinate() {
    return Static_CellOneD::coordinate(index);
}

#if __cplusplus >= 202002L
indexInteger Static_CellOneD::nextRight(indexInteger index) {

  if(index==~(1<<(sizeof(index)*8-1))) return 0;
  return index >> (std::countr_one(index)+1);
}
#else

indexInteger Static_CellOneD::nextRight(indexInteger index) {

    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_ctz only works for ints or compile with cpp20 support");
    unsigned inverted = ~index;
    unsigned firstRight = __builtin_ctz(inverted);
    return index >> (firstRight + 1 );
}
#endif

CellOneD CellOneD::nextRight() {
    return CellOneD(Static_CellOneD::nextRight(index));
}

#if __cplusplus >= 202002L
indexInteger Static_IndexOneD::nextLeft(indexInteger index) {
    if(index==(1<<(sizeof(index)*8-1)))
        return 0;
    return index >> (std::countr_zero(index)+1);
}
#else

indexInteger Static_CellOneD::nextLeft(indexInteger index) {
    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_ctz only works for ints or compile with cpp20 support");
    unsigned firstRight = __builtin_ctz(index);
    return index >> (firstRight + 1 );

}
#endif
CellOneD CellOneD::nextLeft() {
    return CellOneD(Static_CellOneD::nextLeft(index));
}

inline int Static_CellOneD::internDepth(indexInteger index) {
    if(index==0) return -1;
    return (sizeof(index)<<3) - __builtin_clz(index) - 1;
/*    int de = 0;
  while(index!=1) {
     de = de + 1;
     index = index >> 1;
  }
  return de;*/
}


indexInteger Static_CellOneD::nextRight(indexInteger index, int level) {
    int myDepth = internDepth(index);

    level = level - myDepth;
    myAssert(level>=0);

    if(level==0) return nextRight(index);
    index = (index << 1) + 1;
    for(int i=1;i<level;++i) {
        index = index << 1;
    }
    return index;
}

indexInteger Static_CellOneD::nextLeft(indexInteger index, int level) {
    int myDepth = depth(index);

    level = level - myDepth;
    myAssert(level>=0);

    if(level==0) return nextLeft(index);
    index = index << 1;
    for(int i=1;i<level;++i) {
        index = (index << 1) + 1;
    }
    return index;
}


CellOneD CellOneD::nextRightLevel(int level) {
    indexInteger indexRight = Static_CellOneD::nextRight(index,level);
    return CellOneD(Static_CellOneD::nextRight(indexRight,level));
}


class CellDimension : private Static_CellOneD {
public:
    inline CellDimension(); ///> constructs one center point


    inline indexInteger getIndex(int d) const;

    inline int getDepth(int d) const;

    inline bool operator==(const CellDimension I) const;

    inline double coordinate(int d) const;

    inline CellDimension nextLeftLevel(int d, int level);  ///> next left  index on same level in direction d

    inline CellDimension nextRightLevel(int d, int level);    ///> next right index on same level in direction d

    inline void replace(int d, indexInteger newIndex);  ///>  carefully use this function
    inline void replace(int d, CellOneD indRep);  ///>  replaces index in direction d

    inline void PrintCoordinates(){
        for(int d=0; d<DimensionSparseGrid; d++){
            cout <<  coordinate(d) << " " ;
        }
    };

    inline double getLeftBoundary(int d){
        return nextLeftLevel(d, getDepth(d)).coordinate(d);
    }

    inline double getRightBoundary(int d){
        return nextRightLevel(d, getDepth(d)).coordinate(d);
    }



private:

    inline CellDimension(const indexInteger index_org[DimensionSparseGrid]);

    indexInteger index[DimensionSparseGrid];
};


CellDimension::CellDimension(){

    for (int i = 0; i < DimensionSparseGrid; ++i) {
        index[i] = 2; // center point
    }


};

double CellDimension::coordinate(int d) const {
    return Static_CellOneD::coordinate(index[d]);
}


indexInteger CellDimension::getIndex(int d) const {
    assertDimension(d);
    return index[d];
}

int CellDimension::getDepth(int d) const {
    assertDimension(d);
    return Static_CellOneD::depth(index[d]);
}


inline bool CellDimension::operator==(const CellDimension I) const {
    for(int d=0;d<DimensionSparseGrid;++d) {
        if(index[d] != I.index[d]) return false;
    }
    return true;
}


CellDimension CellDimension::nextLeftLevel(int d, int level){
    assertDimension(d);

    CellDimension back(index);

    back.replace(d, Static_CellOneD::nextLeft(index[d],level));

    return back;
}

CellDimension CellDimension::nextRightLevel(int d, int level){
    assertDimension(d);

    CellDimension back(index);

    back.replace(d, Static_CellOneD::nextRight(index[d],level));

    return back;
}



inline void CellDimension::replace(int d, indexInteger newIndex) {
    index[d] = newIndex;
}


void CellDimension::replace(int d, CellOneD indRep) {
    index[d] = indRep.getInteger();
}


CellDimension::CellDimension(const indexInteger index_org[DimensionSparseGrid]) {
    for (int i = 0; i < DimensionSparseGrid; ++i) {
        index[i] = index_org[i];
    }
};


#endif