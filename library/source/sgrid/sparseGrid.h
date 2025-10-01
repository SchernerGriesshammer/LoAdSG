/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef SPARSEGRID_H
#define SPARSEGRID_H

#include "../primes/prime.h"
#include "../indices/index.h"
#include "../iterator/RectangularIterator.h"
#include "omp.h"

#include <iostream>
#include <fstream>
#include "depth.h"
#include "multiDepthHashGrid.h"

class Coordinate;
class VectorSparseG;
class MultiLevelVector;
class ListOfDepthOrderedSubgrids;



class Process {

    friend VectorSparseG;
    friend MultiLevelVector;

private:
    int rank;
    bool all;
    int myrank;
public:

    Process(int *my_rank, int rank_) {
        myrank = *my_rank;
        rank = rank_;
        all = false;
    }

    Process() { all = true; };

    Process(int *my_rank) { myrank = *my_rank, all = true; };

    void runProcess(int rank_) {
        rank = rank_;
        all = false;
    }

    void runALL() { all = true; };

    void setMyRank(int *my_rank) { myrank = *my_rank; }

    int getMyRank() { return myrank; }

    int getRank() { return rank; }


};




template <class A> struct ExprSparseG;

enum type_sparseGrid { sparseGrid_withBoundary, sparseGrid_zeroBoundary };

  
/**
 *
 * adaptive sparse grid 
  \verbatim
                             |       |        |       |
                             |       |        |       | <--- minimalEmpty
                             | prime |        |second | 
  (index,...) ------------>  | table | -----> |table  | <--|
              hashfunction   |       |        |       |    |
                             |       |        |       | ---|
  \endverbatim
 **/

class AdaptiveSparseGrid_Base {




    friend VectorSparseG;
    friend MultiLevelVector;



    friend ListOfDepthOrderedSubgrids;

    //defines a non-member function and makes it a friend of this class
    template<class A, class B>
    friend double product(const ExprSparseG<A> &a, const ExprSparseG<B> &b);

    template<class A>
    friend double L_infty(const ExprSparseG<A> &a);

public:
    //!  Prints coordinates of adaptive sparse grid.     
    void printCoordinates();

    //!  Prints coordinates of adaptive sparse grid in gnu-file.
    void Print_gnu(string name);

    void PrintSlice_gnu(string name);

    void PrintSlice2_gnu(string name);

    void PrintSlice3_gnu(string name);

    void Print_vtk(std::ostream &Datei);

    //! prints active nodes (o) and hanging nodes (x)
    void PrintActiveHanging(int level);

    //! prints all nodes and highlights given Indices
    void PrintGridIndices(int level, IndexDimension* Indices, int numberofindices);


    int getMaxDepth(int d);

    int getMaxDepth(int d, IndexDimension Index);

    inline bool checkMultiDimFiveStencil(IndexDimension Index);

    virtual inline bool AddPoint(IndexDimension I);

    /**
     * schaut nach ob Punkt im Gitter ist
     * @param I index des Punktes
     * @param indexOfData return: index der daten, falls Punkt im Gitter ist
     * @return true, false je nachdem ob Punkt im Gitter ist
     **/
    inline bool occupied(unsigned long &indexOfData, IndexDimension I);

    inline bool occupied2(unsigned long &indexOfData, IndexDimension I);


    inline bool occupied(unsigned long &indexOfData, IndexDimension I, bool active);

    dataInteger *getSecondTable() { return secondTable; };

    dataInteger *getPrimeTable() { return primeTable; };

    bool *getActiveTable() { return isActiveNodeTable; };





    unsigned long getMaximalOccupiedSecondTable() { return maximalOccupiedSecondTable; };


    unsigned long getLengthSecondTable() { return secondTableLength; }


    /**
    * @param i kodierung des Index als langer unsigned long
    * @return IndexDimension "=(i_1,...,i_d)"
    **/
    inline IndexDimension getIndexOfTable(unsigned long i);

    inline IndexDimension getSupportIMin(unsigned long i);

    inline IndexDimension getSupportIMax(unsigned long i);



    /**
* @return unsigned key
**/
    inline unsigned long hash(IndexDimension index);

    inline unsigned long hashWithDepth(IndexDimension index, Depth T);


    Process *mpi;

    bool WorkOnHangingNodes;


    inline bool workonindex(unsigned long i) {
        if (WorkOnHangingNodes || ((!WorkOnHangingNodes) && isActiveNodeTable[i]))
            return true;
        return false;
    }

    unsigned long getPrimeTableLength() const { return primeTableLength; };;;


    MultiDepthHashGrid* getMultiDepthHashGrid() {return multiDepthHashGrid;}

    int getKey(){
        return grid_key;
    }

    int getLevel(){
        return max_LOne-DimensionSparseGrid+1;
    }

    int getMaxLOne(){
        return max_LOne;
    }






protected:
    friend MultiDepthHashGrid;
    MultiDepthHashGrid *multiDepthHashGrid;

    AdaptiveSparseGrid_Base();
    ~AdaptiveSparseGrid_Base();

    unsigned long primeTableLength;
    unsigned long secondTableLength;
    unsigned long minimalEmptySecondTable;
    unsigned long maximalOccupiedSecondTable;


    int *depthTable;
    dataInteger *primeTable;      ///> 0 means empty; v>0 means v-1 is array index
    dataInteger *secondTable;     ///> 0 means empty; 1 means occupied, but no next;  v>1 means v-2 is next array index
    bool *isActiveNodeTable;///> true means AdaptiveSparseGrid_Base::getIndexOfTable (i) is active node.


    indexInteger *indicesSecondTable; ///> Contains coded indices. Length: secondTableLength*DimensionSparseGrid
    indexInteger *indicesSupportMin;
    indexInteger *indicesSupportMax;


    /**
    * @return index k such that secondTable[k] is empty
    **/
    unsigned long getFreeSpaceNumberInSecondTable();


    void setIndexInTable(const IndexDimension I, unsigned long iSetz);

    int grid_key;

    int max_LOne=0;
    int max_level=0;



};

/***
 * Gitter mit Randpunkten
 ***/

/**
 *
 * adaptive sparse grid 
  \verbatim
                             |       |        |       |
                             |       |        |       | <--- minimalEmpty
                             | prime |        |second | 
  (index,...) ------------>  | table | -----> |table  | <--|
              hashfunction   |       |        |       |    |
                             |       |        |       | ---|
  \endverbatim
 **/



class AdaptiveSparseGrid : public AdaptiveSparseGrid_Base {
  public:
    /**
     * @param estimatedMaxNumberOfData geschaetzte maximale Anzahl von Speicher-Variablen
     **/
    AdaptiveSparseGrid() : AdaptiveSparseGrid_Base() {};


    void copy(AdaptiveSparseGrid& newgrid);

    void copy_inner(AdaptiveSparseGrid& newgrid);

    void add_outer(AdaptiveSparseGrid& newgrid);

    void add_RecursiveSonsBoundary(AdaptiveSparseGrid& newgrid);

    void clear();


    /**
     * \* 1.)key=Hashfunction(I);
     * \* 2.)ind = primeTable[key];
     * 
     * First Case: ind == 0 (that means primeTable[key] is free)
     * \verbatim
                      secondTable                           primeTable              Indices-
                                                                                    secondTable        
                          __                                   ___                      ______
        search free      |  |    store f+1                    |   |                    |      |
        space f  ------> |f |   (to make sure,                |   |                    |      |
                         |  |    that primeTable[key] ----->  |f+1| store indices ---> |f     |
                         |  |    isnt zero by accident)       |   | (i_1,..,i_d)       |f+1   |
                         |  |                                 |___|                    |...   |
                         |  |                                                          |f+d-1 |   
                         |  |                                                          |______| 
                         |__|                                                           
                         
      
      \endverbatim
     *
     *Second Case: ind=primeTable[key]!=0 (this key was already used)!\n
     *
     *\* Check if getIndexOfTable(ind-1)=I (remember that data was stored with shift)
     *\* else: Search in secondTable:
      \verbatim
                      secondTable                          
                                                                                           
                          ___   
                         |   |
                         |   |
                         |f  | <--- in secondTable[f] store 1   <----
                         |   |                                     | 
                         |   |                                     |
                         |   |                                     |
                         |   |
                         |ind| ---> 1 is stored because indicesSecondTable[ind] is already used;
                         |__ |     so now search for new freespace f, delete 1 and store f+2
                                   instead (and store indices (i_1,...,i_d) in indicesSecondTable)                                                   
                         
      
      \endverbatim
     **/
    bool AddPoint(const IndexDimension I);

    bool AddPoint(const IndexDimension I, bool hangingNode);

    void AddRecursiveSonsOfPoint(IndexDimension I, int level);





    int getDOFS();

    int getInnerDOFS();

    int getHangingNodes();

    virtual void completeNeumannGrid();

    virtual void completeDirichletGrid();

    virtual void   completeDirichletGrid_NEU();




    inline bool operator==(AdaptiveSparseGrid& second_grid){
        if(grid_key == second_grid.getKey())
            return true;
        else return false;
    }
    void AddPointsOfDepth(Depth T);
    void completeGridWithoutBoundary();
    std::vector<IndexDimension> pointsToAdd;
private:

    void completeGrid(); ///> adds all father points including boundary points


    /**
     * Einfachere Methode als CompleteToLocalTensorProductGrid(). Überprüft Definition. */
    void CompleteToLocalTensorProductGrid();

    void CompleteToLocalTensorProductGrid_NEU(ListOfDepthOrderedSubgrids* listOfDepthOrderedSubgrids);
    void addHangingNodes();




};


////////////////////////////////////////////////////////
// inline functions
////////////////////////////////////////////////////////


unsigned long AdaptiveSparseGrid_Base::hash(IndexDimension index) {
    Depth T(index);
    unsigned long value = index.getIndex(0);
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        value = value + index.getIndex(d) * PrimeNumbers::getPrimeForHash(d);
    }

    value = value + T.at(0) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid);
    for (int d = 1; d < DimensionSparseGrid; d++)
        value = value + T.at(d) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid + d);


    return value % primeTableLength;
}

unsigned long AdaptiveSparseGrid_Base::hashWithDepth(IndexDimension index, Depth T) {

    unsigned long value = index.getIndex(0);
    for (int d = 1; d < DimensionSparseGrid; ++d) {
        value = value + index.getIndex(d) * PrimeNumbers::getPrimeForHash(d);
    }

    value = value + T.at(0) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid);
    for (int d = 1; d < DimensionSparseGrid; d++)
        value = value + T.at(d) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid + d);


    return value % primeTableLength;
}



bool AdaptiveSparseGrid_Base::AddPoint(IndexDimension I) {
   unsigned long indArray =  hash(I);
   unsigned long i = primeTable[indArray];
   if(i == 0) { // in prime table is something free;
      unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
      primeTable[indArray] = freeSpace + 1;
      isActiveNodeTable[freeSpace]=true;
      setIndexInTable(I,freeSpace);
      return true;
   }
   else { // anaylse what is going on 
      i = i-1; // shift since data are stored with shift
      if(I == getIndexOfTable(i)){
          // data is already stored but we still want all these nodes to be active!?
          isActiveNodeTable[i]=true;
          return false;
      } else { // search in second table;
        unsigned long iNext = secondTable[i];
        while(iNext > 1) {
	       i = iNext - 2;
	       if(I == getIndexOfTable(i)){
               isActiveNodeTable[i]=true;
               return false;
           }
	       iNext = secondTable[i];
	    }
	    myAssert(iNext==1);
	    unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
        // shift, since we dont want to store 0 or 1 by accident, i.e. freespace = 1, but 1 means no further link
        secondTable[i] = freeSpace + 2;
        isActiveNodeTable[freeSpace]=true;
        setIndexInTable(I,freeSpace);
	    return true;
      }
   }
}


bool AdaptiveSparseGrid_Base::occupied(unsigned long& indexOfData, IndexDimension I) {

/*
    Depth T(I);
    SingleDepthHashGrid *depthGrid = this->getMultiDepthHashGrid()->tryGetGridForDepth(T);
    if(depthGrid){
        unsigned long j;
        if (depthGrid->occupied(j, I)) {
            indexOfData = j;
            return true;
        } else {
            return false;
        }
    }
    return false;

*/


    unsigned long indArray =  hash(I);
   unsigned long i = primeTable[indArray];
   if(i == 0) { // in prime table is something free;
      return false;
   }
   else { // anaylse what is going on 
      i = i-1; // shift since data are stored with shift
      if(I == getIndexOfTable(i)) {
          indexOfData = i;
          return true;
      }
      else { // search in second table;
          unsigned long iNext = secondTable[i];
          while (iNext > 1) {
              i = iNext - 2;
              if (I == getIndexOfTable(i)) {
                  indexOfData = i;
                  return true;
              }
              iNext = secondTable[i];
          }
          return false;
      }
   }
}




bool AdaptiveSparseGrid_Base::occupied(unsigned long& indexOfData, IndexDimension I, bool active) {
   unsigned long indArray =  hash(I);
   unsigned long i = primeTable[indArray];
   if(i == 0) { // in prime table is something free;
      return false;
   }
   else { // anaylse what is going on 
      i = i-1; // shift since data are stored with shift
      if(I == getIndexOfTable(i)) {
	indexOfData = i;
    active = isActiveNodeTable[i];
	return true;
      }
      else { // search in second table;
        unsigned long iNext = secondTable[i];
        while(iNext > 1) {
	  i = iNext - 2;
	  if(I == getIndexOfTable(i)) {
	     indexOfData = i;
          active = isActiveNodeTable[i];
          return true;
      }
            iNext = secondTable[i];
        }
          return false;
      }
   }
}


////////////////////////////////////////////////

bool AdaptiveSparseGrid_Base::checkMultiDimFiveStencil(IndexDimension Index) {
    unsigned long j;
    occupied(j,Index);
    return isActiveNodeTable[j];




/*    unsigned long k;
    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
        Depth T(Index);
        double value;
        IndexDimension StencilIndex = Index.nextFive_Neumann(&mc, T, &value);
        if (!(occupied(k, StencilIndex)))return false;
    }
    return true;*/
}


////////////////////////////////////////////// inline

inline IndexDimension AdaptiveSparseGrid_Base::getIndexOfTable(unsigned long i) {
    IndexDimension back;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        back.replace(d, indicesSecondTable[d + i * DimensionSparseGrid]);
    }
    return back;
}

inline IndexDimension AdaptiveSparseGrid_Base::getSupportIMin(unsigned long i) {
    IndexDimension back;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        back.replace(d, indicesSupportMin[d + i * DimensionSparseGrid]);
    }
    return back;
}

inline IndexDimension AdaptiveSparseGrid_Base::getSupportIMax(unsigned long i) {
    IndexDimension back;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        back.replace(d, indicesSupportMax[d + i * DimensionSparseGrid]);
    }
    return back;
}

bool AdaptiveSparseGrid_Base::occupied2(unsigned long &indexOfData, IndexDimension I) {
    Depth T(I);
    SingleDepthHashGrid &depthGrid = this->getMultiDepthHashGrid()->getGridForDepth(T);
    const auto &mapping = depthGrid._mapPosToGridPos;
    unsigned long j;
    if(depthGrid.occupied(j,I)){
        indexOfData=mapping[j];
        return true;
    }else{
        return false;
    }


}

#endif  // SPARSEGRID_H

