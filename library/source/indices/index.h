/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
#ifndef INDEXSPARSE_H
#define INDEXSPARSE_H

#include "../myAssert.h"
#include "../myMath/myMath.h"

#include <iostream>
#include <vector>

class Depth;

using namespace std;



enum Direction {
    Left, Right
};
enum Richtung {
    Mitte, Links, Rechts, LinksLinks, RechtsRechts
};



enum elementaryIntervalPoints {
    startUnitInterval, endUnitInterval, centerUnitInterval
};

/**
 * Elementary base class for one dimensional data
 * \verbatim
  Tree structure                          depth
  
         0 --------------- 1                0 
                          /
                         /
                        /
                       /
                      /  
                     /
                    /
                   /
                 10                         1
                /  \
               /    \
            100      101                    2
  
 
     this means:   left  adds 0, right adds 1
               /    \
              /      \
             0        1
             
             
 Indices (and coordinates) 
 
 
          0 --------------- 1                
                          /
                         /
                        /
                       /
                      /  
                     /
                    /
                   /
                  2
                (0.5)                          
                 / \
                /   \
               4     5
            (0.25) (0.75)                    
 
 \endverbatim
 **/


class Static_IndexOneD {
public:
    inline static double coordinate_test(indexInteger index);
    inline static double coordinate(indexInteger index);

protected:
    inline static int depth(indexInteger index);

    inline static indexInteger father(const indexInteger index);

    inline static indexInteger leftSon(const indexInteger index);

    inline static indexInteger rightSon(const indexInteger index);

    inline static indexInteger nextLeft(indexInteger index);   ///> next left index on same level
    inline static indexInteger nextRight(indexInteger index); ///> next right index on same level

    inline static indexInteger nextLeft(indexInteger index, int level);  ///> next left  index on level >= depth
    inline static indexInteger nextRight(indexInteger index, int level); ///> next right index on level >= depth



private:
    inline static int internDepth(indexInteger index);
};

/**
 * Elementary base class for one dimensional data
 * \verbatim
  Tree structure                          depth
  
         0 --------------- 1                0 
                          /
                         /
                        /
                       /
                      /  
                     /
                    /
                   /
                 10                         1
                /  \
               /    \
            100      101                    2
  
 
     this means:   left  adds 0, right adds 1
               /    \
              /      \
             0        1
             
             
 Indices (and coordinates) 
 
 
          0 --------------- 1                
                          /
                         /
                        /
                       /
                      /  
                     /
                    /
                   /
                  2
                (0.5)                          
                 / \
                /   \
               4     5
            (0.25) (0.75)                    
 
 \endverbatim
 **/

class IndexOneD : private Static_IndexOneD {
public:


    explicit inline IndexOneD(elementaryIntervalPoints x);

    inline IndexOneD(indexInteger ind) : index(ind) {};

    inline int depth();

    inline IndexOneD father();

    inline IndexOneD leftSon();

    inline IndexOneD rightSon();
    //! next left index on same level
    /*!
     * Example:
     * * IndexOneD::index = 2;
     * * IndexOneD::nextLeft()=0
     * 
     */
    inline IndexOneD nextLeft();  
    //! next right index on same level
    
    inline IndexOneD nextRight(); 
    //!  next left index on level \f$\geq\f$ depth 
    /*!
     * Example:
     * * IndexOneD::index=2;
     * * IndexOneD::nextLeft(2)=4;
     */
    inline IndexOneD nextLeft(int level);  ///> next left  index on level >= depth
    inline IndexOneD nextRight(int level); ///> next right index on level >= depth
    bool operator==(const IndexOneD I) const { return index==I.index; }
    bool operator!=(const IndexOneD I) const { return index!=I.index; }

    inline double coordinate();

    indexInteger getInteger() const { return index; }


private:
    indexInteger index;

};

/**
 * Berechnet 3^d zur Compile-Zeit
 * **/
template<int Di>
class MaxShift {
public:
    enum {
        value = 3 * MaxShift<Di - 1>::value
    };
};

template<>
class MaxShift<1> {
public:
    enum {
        value = 3
    };
};

/**
 * Berechnet 2^d zur Compile-Zeit
 * **/
template<int Di>
class IntCases {
public:
    enum {
        value = 2 * IntCases<Di - 1>::value
    };
};

template<>
class IntCases<1> {
public:
    enum {
        value = 2
    };
};


/**
 * Richtung fuer 2^d stencil
 *
 * Example: Print nextTwoStencil of indexNow;
 * i.e. if indexNow = (0.25,0.25),  then the code prints:
 *
 * (0.25,0.25),(0.25,0.75),(0.75,0.25),(0.75,0.75)
 * @code
 *
 * IndexDimension indexNow;
    for(MultiDimCompass mc;mc.goon();++mc) {
    IndexDimension J = indexNow.nextTwoStep(&mc, T);
    J.PrintCoord();
    cout << endl;
    }
 * @endcode
 **/
class MultiDimTwoCompass {
public:
    MultiDimTwoCompass() { shiftNumber = 0; }

    inline Richtung getRichtung(unsigned int d) const;

    void operator++() { ++shiftNumber; }

    unsigned int getMaxShift() { return maxShift; }

    unsigned int getShiftNumber() { return shiftNumber; };

    bool goon() { return shiftNumber < maxShift; }



private:
    unsigned int shiftNumber;
    static unsigned int maxShift;
};



/**
 * Richtung fuer 3^d stencil
 *
 * Example: Print nextTwoStencil of indexNow;
 * i.e. if indexNow = (0.25,0.25),  then the code prints:
 *
 * (0.25,0.25),(0.25,0.75),(0.75,0.25),(0.75,0.75)
 * @code
 *
 * IndexDimension indexNow;
    for(MultiDimCompass mc;mc.goon();++mc) {
    IndexDimension J = indexNow.nextTwoStep(&mc, T);
    J.PrintCoord();
    cout << endl;
    }
 * @endcode
 **//*

class MultiDimCompass {
public:
    MultiDimCompass() { shiftNumber = 0; }

    inline Richtung getRichtung(unsigned int d) const;

    void operator++() { ++shiftNumber; }

    unsigned int getMaxShift() { return maxShift; }

    unsigned int getShiftNumber() { return shiftNumber; };

    bool goon() { return shiftNumber < maxShift; }

    bool operator==(MultiDimTwoCompass& mc){
        for(int d=0; d<DimensionSparseGrid; d++){
            if(getRichtung(d)==mc.getRichtung(d))
                return true;
        }
        return false;

    }

    void PrintRichtung(int d){
        if(getRichtung(d)==Mitte)
            cout << "Mitte ";
        if(getRichtung(d)==Links)
            cout << "Links ";
        if(getRichtung(d)==Rechts)
            cout << "Rechts ";
    }

private:
    unsigned int shiftNumber;
    static unsigned int maxShift;
};
*/
/**
 * Richtung fuer 3^d stencil
 *
 * Example: Print nextTwoStencil of indexNow;
 * i.e. if indexNow = (0.25,0.25),  then the code prints:
 *
 * (0.25,0.25),(0.25,0.75),(0.75,0.25),(0.75,0.75)
 * @code
 *
 * IndexDimension indexNow;
    for(MultiDimCompass mc;mc.hasNext();++mc) {
    IndexDimension J = indexNow.nextTwoStep(&mc, T);
    J.PrintCoord();
    cout << endl;
    }
 * @endcode
 **/
class MultiDimCompass {
public:
    MultiDimCompass() {
        shiftNumber = 0;
        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            richtungen[i]=0;
        }
      /*  Richtungen = new unsigned int [DimensionSparseGrid*maxShift];

        for(int shiftN=0;shiftN<maxShift;shiftN++){

            for(size_t d=0; d<DimensionSparseGrid; d++){
                Richtungen[shiftN*DimensionSparseGrid+d]=richtungen[d];

            }


            for (size_t i = 0; i < DimensionSparseGrid; i++)
            {

                richtungen[i]=richtungen[i]+1;
                if(richtungen[i]==3){
                    richtungen[i] = 0;

                    goto goon_inner;
                }
                goto goon_outer;
                goon_inner:;
            }

            goon_outer:;
        }
*/
    }



    inline Richtung getRichtung (unsigned int d) const;

    void operator++() {
        ++shiftNumber;
        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            richtungen[i]=richtungen[i]+1;
            if(richtungen[i]==3){
                richtungen[i] = 0;
                continue;
            }
            break;
        }

    }

    static unsigned int getMaxShift() { return maxShift; }

    unsigned int getShiftNumber() { return shiftNumber; };

    bool hasNext() { return shiftNumber < maxShift; }

    bool goon() { return shiftNumber < maxShift; }

private:
    unsigned int richtungen[DimensionSparseGrid];
    unsigned int shiftNumber;
    static unsigned int maxShift;

};




/**
 * Anwendung:
 * @code
 *         for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
            IndexDimension J = indexNow.nextFive(&mc, T);
            }
  @endcode
 */
class MultiDimFiveCompass {
public:
    MultiDimFiveCompass() { shiftNumber = 0;
        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            directions[i]=0;
        }
    }

    inline Richtung getRichtung(unsigned int d){
        return (Richtung) directions[d];
    };

/*
    inline Richtung getRichtung(unsigned int d){
        unsigned int num = shiftNumber;
        for (int s = 0; s < d; ++s) num = num / 5;
        return (Richtung) (num % 5);
    };
*/

    void operator++() { ++shiftNumber;

        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            directions[i]=directions[i]+1;
            if(directions[i]==5){
                directions[i] = 0;
                continue;
            }
            break;
        }
    }

    unsigned int getMaxShift() { return maxShift; }
    bool hasNext() { return shiftNumber < maxShift; }

    bool goon() { return shiftNumber < maxShift; }

    void goToStart(){ shiftNumber = 0;}

    int getShiftNumber(){return shiftNumber;}

private:
    unsigned int shiftNumber;
    static unsigned int maxShift;
    unsigned int directions[DimensionSparseGrid];


};



class IndexDimension : private Static_IndexOneD {
public:
    inline IndexDimension(); ///> constructs one center point
    inline indexInteger getIndex(int d) const;

    inline int getDepth(int d) const;


    inline double coordinate(int d) const;

    inline IndexDimension nextLeft(int d);  ///> next left  index on same level in direction d
    inline IndexDimension nextLeft(int d) const;

    inline IndexDimension nextRight(int d) const; ///> next right index on same level in direction d
    inline IndexDimension nextRight(int d);

    inline IndexDimension nextLeft(int d, int level) const;  ///> next left  index on level in direction d
    inline IndexDimension nextRight(int d, int level) const;  ///> next right index on level in direction d
    inline IndexDimension father(int d);

    /**
     *
     * @param mc MultiCompass
     * @param T Tiefe Object
     * @param w Weite
     * @return  IndexDimension
     */
    IndexDimension
    nextTwoStep(MultiDimCompass *mc, Depth T) const;  ///> next index on level T

    IndexDimension
    nextFive(MultiDimFiveCompass *mc, Depth T) const;  ///> next index on level T

    IndexDimension
    nextFive(MultiDimFiveCompass *mc, Depth T,
             double *stencilvalue) const;  ///> next index on level T and corresponding prewavelet stencil value

    IndexDimension
    nextFiveP(MultiDimFiveCompass *mc, Depth& T,
             double *stencilvalue) const;

    IndexDimension
    nextFive(MultiDimFiveCompass *mc, Depth T,
             double *stencilvalue,
             bool *last) const;  ///> next index on level T and corresponding prewavelet stencil value

    inline IndexDimension
    nextFive(MultiDimFiveCompass *mc, int *T,
             double *stencilvalue) const;


    IndexDimension
    nextFive_Neumann(MultiDimFiveCompass *mc, Depth &T, double *stencilvalue) const;

    IndexDimension
    nextFive_Neumann(MultiDimFiveCompass *mc, Depth &T, double *stencilvalue, bool *last) const;

    inline IndexDimension
    nextFive_Neumann(MultiDimFiveCompass *mc, int *Tiefen, double *stencilvalue) const;


    IndexDimension
    nextThree_Stencil(MultiDimCompass *mc, Depth T, double *stencilvalue) const;

    inline IndexDimension
    nextThree(MultiDimCompass *mc, int *Tiefen) const;

    inline IndexDimension
    nextThree2(MultiDimCompass *mc, int *Tiefen) const;

    IndexDimension
    nextThree_boundary(MultiDimCompass *mc, Depth &T) const;

    inline IndexDimension
    nextThree_boundary(MultiDimCompass *mc, int *Tiefen) const;


    inline bool wasLeftChild(int d);

    inline bool wasRightChild(int d);

    inline bool isLeftGrandChildOf(IndexDimension Index);

    inline IndexDimension leftSon(int d) const;

    inline IndexDimension rightSon(int d) const;

    inline void replace(int d, indexInteger newIndex);  ///>  carefully use this function
    inline void replace(int d, IndexOneD indRep);  ///>  replaces index in direction d
    inline bool operator==(const IndexDimension I) const;

    inline bool operator!=(const IndexDimension I) const;

    inline void operator=(const IndexOneD indRep);

    inline void Print() const;

    inline void PrintCoord() const;

    inline bool isNotAtBoundary() const;

    inline bool isNotAtBoundary(int d) const;

    inline bool isOutOfBoundary(
            int d) const; ///> needed for prewavelet grids; otherwise grid points with index 3 or others with first two entrys 11 would be allowed.

    inline bool isAtRightBoundary(int d) const;

    inline bool isAtLeftBoundary(int d) const;

    inline bool isAtBoundary() const;

    inline bool isLinksRandNah(int d) const;

    inline bool isRechtsRandNah(int d) const;

    static IndexDimension Maximum(IndexDimension &A, IndexDimension &B);

    static IndexDimension Minimum(IndexDimension &A, IndexDimension &B);

    static IndexDimension Maximum(IndexDimension &max, IndexDimension &B, int *changes);

    static IndexDimension Minimum(IndexDimension &min, IndexDimension &B, int *changes);

    inline Richtung isRandnah(int d);

    inline bool isRandnah(Richtung* richtung);

private:
    inline IndexDimension(const indexInteger index_org[DimensionSparseGrid]);

    indexInteger index[DimensionSparseGrid];
};




////////////////////////////////////////////////

inline Richtung MultiDimTwoCompass::getRichtung(unsigned int d) const{
    unsigned int num = shiftNumber;
    for (int s = 0; s < d; ++s) num = num / 2;
    return (Richtung) ((num % 2)+1);
}

inline Richtung MultiDimCompass::getRichtung(unsigned int d) const{
    return (Richtung) richtungen[d];
    //return (Richtung) Richtungen[shiftNumber*DimensionSparseGrid+d];
/*    unsigned int num = shiftNumber;
    for (int s = 0; s < d; ++s) num = num / 3;
    return (Richtung) (num % 3)*/;
}

/*inline Richtung MultiDimFiveCompass::getRichtung(unsigned int d) {
    unsigned int num = shiftNumber;
    for (int s = 0; s < d; ++s) num = num / 5;
    return (Richtung) (num % 5);
}*/

IndexOneD::IndexOneD(elementaryIntervalPoints x) {

    if (x == startUnitInterval) index = 0;
    else if (x == endUnitInterval) index = 1;
    else index = 2;
}

//---------

indexInteger Static_IndexOneD::father(const indexInteger index) {
   return index>>1;
}

IndexOneD IndexOneD::father() {
   return  IndexOneD(Static_IndexOneD::father(index));
}

//---------

indexInteger Static_IndexOneD::leftSon(const indexInteger index) {
   return index<<1;
}

IndexOneD IndexOneD::leftSon() {
   return  IndexOneD(Static_IndexOneD::leftSon(index));
}

//---------

indexInteger Static_IndexOneD::rightSon(const indexInteger index) {
   return ((index<<1)+1);
}

IndexOneD IndexOneD::rightSon() {
   return IndexOneD(Static_IndexOneD::rightSon(index));
}

//---------

double Static_IndexOneD::coordinate_test(indexInteger index) {
   double co = 0.0;
   indexInteger newindex=index;
   while(index>0) {
     co = co * 0.5;
     if((index&1)==1) co = co + 1.0;
     else           co = co - 1.0;
     index = index>>1;   
//cout << " co: " << co << endl;
   }
/*   double newco = coordinate_test(newindex);
   if(abs(newco-co)>0){
       cout<< "error new coordinate" << endl;
       cout << newindex << "  " << co << "  " << newco << "  " << abs(newco-co) << endl;
       exit(0);
   }*/
   return co;
}

inline int subtractFirstBit(int num) {
    // Find the position of the most significant bit
    int msb_pos = 31 - __builtin_clz(num); // Assuming a 32-bit integer, use appropriate value for other sizes

    // Create a mask with all bits set except for the most significant bit
    int mask = ~(1 << msb_pos);

    // Perform bitwise AND operation to subtract the first bit
    return num & mask;
}


double Static_IndexOneD::coordinate(indexInteger index) {
    if(index==0)return 0;
    int d = depth(index);
    //cout << d << endl;

    double h = 1.0/double(POW2(d));
    //cout << h << endl;
    indexInteger n= subtractFirstBit(index);
    //cout << n << endl;
    return  double(2*n+1)*h;


}


double IndexOneD::coordinate() {
  return Static_IndexOneD::coordinate(index);
}

//---------

//---------

#if __cplusplus >= 202002L
indexInteger Static_IndexOneD::nextLeft(indexInteger index) {
    if(index==(1<<(sizeof(index)*8-1)))
        return 0;
    return index >> (std::countr_zero(index)+1);
}
#else

indexInteger Static_IndexOneD::nextLeft(indexInteger index) {
    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_ctz only works for ints or compile with cpp20 support");
    unsigned firstRight = __builtin_ctz(index);
    return index >> (firstRight + 1 );

}
#endif
IndexOneD IndexOneD::nextLeft() {
    return IndexOneD(Static_IndexOneD::nextLeft(index));
}

//--------

//--------


// static indexInteger nextRight_original(indexInteger index) {
//   if((index&1)==0) return (index>>1);
//   index = index>>1;
//   while((index&1)==1) {
//      index = index >> 1;
//   }
//   return (index>>1);
// }

#if __cplusplus >= 202002L
indexInteger Static_IndexOneD::nextRight(indexInteger index) {
  if(index==~(1<<(sizeof(index)*8-1))) return 0;
  return index >> (std::countr_one(index)+1);
}
#else

indexInteger Static_IndexOneD::nextRight(indexInteger index) {

    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_ctz only works for ints or compile with cpp20 support");
    unsigned inverted = ~index;
    unsigned firstRight = __builtin_ctz(inverted);
    return index >> (firstRight + 1 );
}
#endif

IndexOneD IndexOneD::nextRight() {
    return IndexOneD(Static_IndexOneD::nextRight(index));
}

//---------
#if __cplusplus >= 202002L
inline int Static_IndexOneD::depth(indexInteger index) {
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

inline int Static_IndexOneD::depth(indexInteger index) {
    static_assert(sizeof(indexInteger) == 4,
                  "__builtin_clz only works for ints or compile with cpp20 support");
    if(index==0) return 0;
    return __builtin_clz(index) ^ 31;
}
#endif

inline int IndexOneD::depth() {
    return Static_IndexOneD::depth(index);
}
//---------

inline int Static_IndexOneD::internDepth(indexInteger index) {
  if(index==0) return -1;
    return (sizeof(index)<<3) - __builtin_clz(index) - 1;
/*    int de = 0;
  while(index!=1) {
     de = de + 1;
     index = index >> 1;
  }
  return de;*/
}


//------  

indexInteger Static_IndexOneD::nextLeft(indexInteger index, int level) {
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

IndexOneD IndexOneD::nextLeft(int level) {
   return IndexOneD(Static_IndexOneD::nextLeft(index,level));
}

//------

indexInteger Static_IndexOneD::nextRight(indexInteger index, int level) {
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


IndexOneD IndexOneD::nextRight(int level) {
    return IndexOneD(Static_IndexOneD::nextRight(index,level));
}

//-------------------------------------------------------------------------
// IndexDimension
//-------------------------------------------------------------------------

IndexDimension::IndexDimension(const indexInteger index_org[DimensionSparseGrid]) {
    for (int i = 0; i < DimensionSparseGrid; ++i) {
        index[i] = index_org[i];
    }
};


void IndexDimension::Print() const {
   cout << " Print index for test: " << endl;
   for(int i=0;i<DimensionSparseGrid;++i) {
     cout << index[i] << endl;
  }
}

void IndexDimension::PrintCoord() const {
//   cout << " Print coordinate: " << endl;
   cout << "(" << Static_IndexOneD::coordinate(index[0]);
   for(int i=1;i<DimensionSparseGrid;++i) {
       cout << ", " << Static_IndexOneD::coordinate(index[i]);
  }
  cout << ")";
}


IndexDimension::IndexDimension() {
    for (int i = 0; i < DimensionSparseGrid; ++i) {
        index[i] = 2; // center point
    }


};

double IndexDimension::coordinate(int d) const {
  return Static_IndexOneD::coordinate(index[d]);
}

IndexDimension IndexDimension::leftSon(int d) const{
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::leftSon(index[d]));
   
   return back;
}
  
IndexDimension IndexDimension::rightSon(int d) const{
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::rightSon(index[d]));
   
   return back;
}
    
    
bool IndexDimension::wasLeftChild(int d){
    if((index[d]&1)==1) return false;
    return true;
}    


bool IndexDimension::wasRightChild(int d){
    if((index[d]&1)==1) return true;
    return false;
}




//---------------------   
    
IndexDimension IndexDimension::nextLeft(int d) const {
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextLeft(index[d]));
   
   return back;
}

IndexDimension IndexDimension::nextLeft(int d){
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextLeft(index[d]));
   
   return back;
}



IndexDimension IndexDimension::nextRight(int d) {
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextRight(index[d]));
   
   return back;  
}


IndexDimension IndexDimension::nextRight(int d) const {
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextRight(index[d]));
   
   return back;  
}
    
    
    
//---------------------   
    
    
    
IndexDimension IndexDimension::nextLeft(int d, int level) const{
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextLeft(index[d],level));
   
   return back;
}



IndexDimension IndexDimension::nextRight(int d, int level) const{
   assertDimension(d);

   IndexDimension back(index);
   back.replace(d, Static_IndexOneD::nextRight(index[d],level));
   
   return back;  
}

IndexDimension IndexDimension::father(int d){
    assertDimension(d);
    
    IndexDimension back(index);
    back.replace(d, Static_IndexOneD::father(index[d]));
    
    return back;
}

    
    
//---------------------   





inline void IndexDimension::replace(int d, indexInteger newIndex) {
    index[d] = newIndex;
}
 
 
void IndexDimension::replace(int d, IndexOneD indRep) {
   index[d] = indRep.getInteger();
}
      
//------------------      
      
bool IndexDimension::isNotAtBoundary() const {
   for(int d=0;d<DimensionSparseGrid;++d) {
       if(index[d]==0) return false;
       if(index[d]==1) return false;       
   }
   return true;
}

bool IndexDimension::isNotAtBoundary(int d) const {
       if(index[d]==0) return false;
       if(index[d]==1) return false;       
       return true;
}

  
bool IndexDimension::isOutOfBoundary(int d) const {
     int x =index[d];
     
     if (x == 3){return true;}
     if (x == 2) {return false;}
    
     int y = x;
    int bits;
    for (bits = 0; x != 0; ++bits) { x = x >> 1; }
    bits = bits - 2;
    y = y >> bits;

    if (y == 3) { return true; }
    if (y == 2) { return false; }

    return true;
}


bool IndexDimension::isAtRightBoundary(int d) const {
    if (index[d] == 1) return true;
    return false;
}

bool IndexDimension::isAtLeftBoundary(int d) const {
    if (index[d] == 0) return true;
    return false;
}

bool IndexDimension::isAtBoundary() const {
    for (int d = 0; d < DimensionSparseGrid; d++)
        if (index[d] == 0 || index[d] == 1) return true;
    return false;
}

bool IndexDimension::isLinksRandNah(int d) const {
    if (nextLeft(d).isAtLeftBoundary(d)) return true;
    return false;
}

bool IndexDimension::isRechtsRandNah(int d) const {
    if (nextRight(d).isAtRightBoundary(d)) return true;
    return false;
}



//-------------------------



void IndexDimension::operator=(const IndexOneD indRep) {
    indexInteger newIndex = indRep.getInteger();
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        index[d] = newIndex;
    }
}

indexInteger IndexDimension::getIndex(int d) const {
    assertDimension(d);
    return index[d];
}

int IndexDimension::getDepth(int d) const {
    assertDimension(d);
    return Static_IndexOneD::depth(index[d]);
}


inline bool IndexDimension::operator==(const IndexDimension I) const {
  for(int d=0;d<DimensionSparseGrid;++d) {
     if(index[d] != I.index[d]) return false;
  }
    return true;
}

inline bool IndexDimension::operator!=(const IndexDimension I) const {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (index[d] != I.index[d]) return true;
    }
    return false;
}

Richtung IndexDimension::isRandnah(int d) {
    if (nextLeft(d).isAtLeftBoundary(d))
        return Links;
    if (nextRight(d).isAtRightBoundary(d)) return Rechts;
    return Mitte;
}

bool IndexDimension::isRandnah(Richtung* richtung) {
    bool returnbool = false;
    for(int d=0; d<DimensionSparseGrid; d++) {
        if (nextLeft(d).isAtLeftBoundary(d)){
            richtung[d]=Links;
            returnbool=true;
            continue;
        }

        if (nextRight(d).isAtRightBoundary(d)){
            richtung[d]=Rechts;
            returnbool=true;
            continue;
        }

        richtung[d]=Mitte;

    }
    return returnbool;
}

//IndexDimension::IndexDimension(int d, indexInteger newIndex); 
IndexDimension IndexDimension::nextThree_boundary(MultiDimCompass *mc, int *Tiefen) const {

    IndexDimension J = (*this);

    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = Tiefen[d];
            Richtung r = mc->getRichtung(d);


            if (r == Links) {
                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
            }


            if (r == Rechts) {
                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        J = (*this);
        mc->operator++();
    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
    return J;
}


IndexDimension IndexDimension::nextFive_Neumann(MultiDimFiveCompass *mc, int *Tiefen, double *stencilvalue) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {
        for (int d = 0; d < DimensionSparseGrid; ++d) {
            int t = Tiefen[d];
            Richtung r = mc->getRichtung(d);


            if (r == Links) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
                if (J.isAtLeftBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtLeftBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;

            }

            if (r == LinksLinks) {
                if (J.getDepth(d) < 2) goto startagain;
                if (J.isAtLeftBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t).nextLeft(d, t);
                val = 0.1 * val;

            }

            if (r == Rechts) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
                if (J.isAtRightBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtRightBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;
            }
            if (r == RechtsRechts) {
                if (J.getDepth(d) < 2) goto startagain;

                if (J.isAtRightBoundary(d) || J.nextRight(d, t).isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t).nextRight(d, t);
                val = 0.1 * val;


            }

            if (r == Mitte) {


                if (J.getDepth(d) < 2) val = 1.0 * val;
                else {
                    if (J.nextLeft(d, t).isAtLeftBoundary(d) && !J.nextRight(d, t).isAtRightBoundary(d))
                        val = 1.1 * val;
                    if (J.nextRight(d, t).isAtRightBoundary(d) && !J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 1.1 * val;
                }
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.

    if (mc->goon() == false && J == *this) val = 0.0;

    *stencilvalue = val;

    return J;
}

inline IndexDimension IndexDimension::nextThree(MultiDimCompass *mc, int *Tiefen) const {

    IndexDimension J = (*this);

    while (mc->goon()) {

        for (size_t d = 0; d < DimensionSparseGrid; ++d) {

            int t = Tiefen[d];
            Richtung r = mc->getRichtung(d);
            //if(J.getDepth(d)==1) r = Mitte;



            if (r == Links) {
                if (J.getDepth(d) < 1) goto startagain;

                //if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);


            }


            if (r == Rechts) {
                if (J.getDepth(d) < 1) goto startagain;

                //if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);

            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:

        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.

    return J;
}


inline IndexDimension IndexDimension::nextThree2(MultiDimCompass *mc, int *Tiefen) const {

    IndexDimension J = (*this);
    for (size_t d = 0; d < DimensionSparseGrid && mc->goon(); ++d) {

            int t = Tiefen[d];
            Richtung r = mc->getRichtung(d);


            if (J.getDepth(d) < 1 && r !=Mitte){
                mc->operator++();
                d--;
                continue;
            }


            if (r == Links)
                J = J.nextLeft(d, t);



            if (r == Rechts)
                J = J.nextRight(d, t);




    }

    return J;
}

inline IndexDimension IndexDimension::nextFive(MultiDimFiveCompass *mc, int *T, double *stencilvalue) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {
        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T[d];
            Richtung r = mc->getRichtung(d);

            if (r == Links) {
                if (J.nextLeft(d, t).isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = J.nextLeft(d, t);
                } else // break outer for loop
                    goto startagain;

            }

            if (r == LinksLinks) {
                if (J.nextLeft(d, t).isNotAtBoundary(d)) {
                    J = J.nextLeft(d, t).nextLeft(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }
            if (r == Rechts) {
                if (J.nextRight(d, t).isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = J.nextRight(d, t);
                } else goto startagain;
            }
            if (r == RechtsRechts) {
                if (J.nextRight(d, t).isNotAtBoundary(d)) {
                    J = J.nextRight(d, t).nextRight(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }

            if (r == Mitte) {
                if (J.getDepth(d) > 1)
                    if (J.nextRight(d, t).isAtRightBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 0.9 * val;
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
    if (mc->goon() == false && J == *this) val = 0.0;
    *stencilvalue = val;

    return J;
}


#endif  // INDEXSPARSE_H

