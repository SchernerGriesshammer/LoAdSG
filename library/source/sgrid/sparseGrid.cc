/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#include "../abbrevi.h"
#include "../myAssert.h"
#include "../indices/index.h"
#include "../primes/prime.h"

#include "sparseGrid.h"
#include "GridGeneration.h"
#include "../iterator/RectangularIterator.h"
#include "ListOfDepthOrderedGrids.h"
#include "../mympi.h"
#include "../iterator/depthIterator.h"

#include <iostream>
#include <sys/time.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////
// 1. Allgemeine Gitter
////////////////////////////////////////////////////////////////////////////////////////////
int key = 0;
 
AdaptiveSparseGrid_Base::AdaptiveSparseGrid_Base() {
    grid_key = key;
    key++;


    multiDepthHashGrid = new MultiDepthHashGrid(*this);


    WorkOnHangingNodes = false;


    static Process p;
    mpi = &p;


    primeTableLength = PrimeNumbers::getNextPrime(1500000);
    secondTableLength = 1500000;                                // Riccarda, Christoph:  Das müsst Ihr mal mit resize verbessern oder nicht??
    //numberOfData = estimatedMaxNumberOfData;
    //newNumberForData = -1;
   
   primeTable = new dataInteger[primeTableLength]; 
   for(unsigned long i=0;i<primeTableLength;++i) primeTable[i] = 0;

   secondTable = new dataInteger[secondTableLength]; 
   for(unsigned long i=0;i<secondTableLength;++i) secondTable[i] = 0;
   
   isActiveNodeTable= new bool[secondTableLength];
    for (unsigned long i = 0; i < secondTableLength; ++i) isActiveNodeTable[i] = false;



    indicesSecondTable = new indexInteger[secondTableLength *
                                          DimensionSparseGrid];///>for every point we need to store three indices i.e. i_x and i_y
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSecondTable[i] = 0;

    indicesSupportMin = new indexInteger[secondTableLength *
                                          DimensionSparseGrid];///>for every point we need to store three indices i.e. i_x and i_y
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMin[i] = 0;

    indicesSupportMax = new indexInteger[secondTableLength *
                                          DimensionSparseGrid];///>for every point we need to store three indices i.e. i_x and i_y
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMax[i] = 0;


    minimalEmptySecondTable = 0;
    maximalOccupiedSecondTable = 0;
}









AdaptiveSparseGrid_Base::~AdaptiveSparseGrid_Base() {
    delete[] primeTable;
    delete[] secondTable;
    delete[] isActiveNodeTable;

    delete[] indicesSecondTable;

    delete[] indicesSupportMin;
    delete[] indicesSupportMax;

    delete multiDepthHashGrid;



}


/*
int AdaptiveSparseGrid_Base::getNewDataNumber() {
    newNumberForData++;
    if(newNumberForData<numberOfData) return newNumberForData;
    
    // zweiter hashtable muss vergroesert werden
    numberOfData++;
    double* dataTable_new = new double[numberOfData * secondTableLength];
    for(unsigned long i=0;i<secondTableLength;++i) {
        for(int j=0;j<numberOfData-1;++j) {
	    // kopieren vom alten hash table
	    dataTable_new[i * numberOfData + j] = dataTable[i * (numberOfData-1) + j];
	}
	// neuer default Wert 0.0
	dataTable_new[i * numberOfData + (numberOfData-1)] = 0.0;
    }
    delete[] dataTable;
    dataTable = dataTable_new;
    return newNumberForData;
}
*/
    
void AdaptiveSparseGrid_Base::printCoordinates() {
  cout << " Print coordinates of adaptive sparse grid: " << endl;
  cout << " ------------------------------------------ " << endl;
  for(unsigned long  i=0;i<maximalOccupiedSecondTable;++i) {
      if(secondTable[i]!=0) {
	 for(int d=0;d<DimensionSparseGrid;++d) {
	     cout <<  IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << ",  ";
	 }
	 cout << endl;
      }
  }
}


void AdaptiveSparseGrid_Base::PrintActiveHanging(int level) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << "     " << "\n";
    cout << "active (o) and hanging (x) nodes on grid " << endl;
    cout << "------------------------- " << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);

    for (int d = 2; d < DimensionSparseGrid; d++)
        I.replace(d, centerUnitInterval);


    unsigned long indexOfData;
    //unsigned int numberOfData = sparseGrid->numberOfData;
    //double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (occupied(indexOfData, I)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";	      

                    if (isActiveNodeTable[indexOfData]) {
                        cout << " o ";
                    } else {
                        cout << " x ";

                    //if(typeOfHangingNode[indexOfData]==0){cout << " x ";}
                    //else cout << " + ";
                }


            } else {
                if(I.isAtBoundary())
                    cout << " + ";
                else
                    cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        std::cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        if (DimensionSparseGrid > 1) I = I.nextRight(1, level);
        else break;

    }
}

void AdaptiveSparseGrid_Base::PrintGridIndices(int level, IndexDimension* Indices, int numberofindices){
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << "     " << "\n";
    cout << "active (o), hanging (x) and neumann-hanging (+) nodes on grid " << endl;
    cout << "------------------------- " << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);

    for (int d = 2; d < DimensionSparseGrid; d++)
        I.replace(d, centerUnitInterval);


    unsigned long indexOfData;
    //unsigned int numberOfData = sparseGrid->numberOfData;
    //double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (occupied(indexOfData, I)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";
                bool goon = true;
                for(int k=0; k < numberofindices; k++){
                    if( I == Indices[k]) {
                        cout << " O ";
                        goon = false;
                    }
                }
                if(goon){
                    if (isActiveNodeTable[indexOfData]) {
                        cout << " o ";
                    } else {
                        cout << " x ";

                        //if(typeOfHangingNode[indexOfData]==0){cout << " x ";}
                        //else cout << " + ";
                    }
                }


            } else {
                bool goon = true;
                for(int k=0; k < numberofindices; k++){
                    if( I == Indices[k]) {
                        cout << " O ";
                        goon = false;
                    }
                }
                if(goon) {
                    if (I.isAtBoundary())
                        cout << " + ";
                    else
                        cout << "   ";
                }
            }
            I = I.nextRight(0, level);
        }
        std::cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        if (DimensionSparseGrid > 1) I = I.nextRight(1, level);
        else break;

    }
}



void AdaptiveSparseGrid_Base::Print_gnu(string name) {

    std::ofstream Datei;

    Datei.open(name, std::ios::out);


    Datei << "# Print coordinates of adaptive sparse grid: " << endl;
    for (unsigned long i = 0; i < maximalOccupiedSecondTable; ++i) {
        if (secondTable[i] != 0) {
            for (int d = 0; d < DimensionSparseGrid; ++d) {
                Datei << IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
            }
            if (getActiveTable()[i]) Datei << 1 << endl;
            else Datei << 0 << endl;
        }

    }

    Datei.close();
}


void AdaptiveSparseGrid_Base::PrintSlice_gnu(string name) {

    std::ofstream Datei;

    Datei.open(name, std::ios::out);


    Datei << "# Print coordinates of adaptive sparse grid: " << endl;
    for (unsigned long i = 0; i < maximalOccupiedSecondTable; ++i) {
        if (secondTable[i] != 0) {
            IndexDimension I = getIndexOfTable(i);

            if(I.getIndex(0)== I.getIndex(2) && I.getIndex(1)==I.getIndex(3)) {
                for (int d = 0; d < 2; ++d) {
                    Datei << IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
                }
                if (getActiveTable()[i]) Datei << 1 << endl;
                else Datei << 0 << endl;
            }
        }

    }

    Datei.close();
}
void AdaptiveSparseGrid_Base::PrintSlice2_gnu(string name) {

    std::ofstream Datei;

    Datei.open(name, std::ios::out);


    Datei << "# Print coordinates of adaptive sparse grid: " << endl;
    for (unsigned long i = 0; i < maximalOccupiedSecondTable; ++i) {
        if (secondTable[i] != 0) {
            IndexDimension I = getIndexOfTable(i);

            if(abs(I.coordinate(0)-(-I.coordinate(2)+1.0))>1e-10 && abs(I.coordinate(1)-(-I.coordinate(3)+1.0))>1e-10) {
                for (int d = 0; d < 2; ++d) {
                    Datei << IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
                }
                if (getActiveTable()[i]) Datei << 1 << endl;
                else Datei << 0 << endl;
            }
        }

    }

    Datei.close();
}
void AdaptiveSparseGrid_Base::PrintSlice3_gnu(string name) {

    std::ofstream Datei;

    Datei.open(name, std::ios::out);


    Datei << "# Print coordinates of adaptive sparse grid: " << endl;
    for (unsigned long i = 0; i < maximalOccupiedSecondTable; ++i) {
        if (secondTable[i] != 0) {
            IndexDimension I = getIndexOfTable(i);

            if(I.getIndex(0)==2 &&I.getIndex(1)==2) {
                for (int d = 0; d < DimensionSparseGrid; ++d) {
                    if(d==2 || d==3 )
                    Datei << IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
                }
                if (getActiveTable()[i]) Datei << 1 << endl;
                else Datei << 0 << endl;
            }
        }

    }

    Datei.close();
}



void AdaptiveSparseGrid_Base::Print_vtk(std::ostream& Datei) {
  // Teil 0: Write information
  Datei << "# vtk DataFile Version 2.0\n"
        << "SparseGrid" << endl
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";

  // Teil 1: Kopfzeile schreiben
  int num_total = 0;
    for(unsigned long  i=0;i<maximalOccupiedSecondTable;++i) {
      if (secondTable[i]!=0) num_total++;
  }
  Datei << "POINTS " << num_total << " float\n";

  
  // Teil 2: Koordinaten der Punkte ausgeben
  
  for(unsigned long  i=0;i<maximalOccupiedSecondTable;++i) {
      if(secondTable[i]!=0) {
	 for(int d=0;d<DimensionSparseGrid;++d) {
	     Datei <<  IndexOneD(indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
	 }
	 Datei << endl;
      }
  }
  // Teil 3: Zellen ausgeben
  Datei << "\nCELLS " << num_total << " " <<  num_total*2  << "\n";
  for(int i = 0; i < num_total; i++)
      Datei << 1 << " " << i << endl;
  Datei << "\nCELL_TYPES " << num_total << "\n";
 for(int i = 0; i < num_total; i++)
      Datei << "1" << endl;
 // Teil 4:
 Datei << "\nPOINT_DATA " << num_total << "\n";
 Datei << "SCALARS " << "SparseGrid" << " float 1\n";
 Datei << "LOOKUP_TABLE default\r\n";
 for(int i = 0; i < num_total; i++)
      Datei << 0.5 << endl;
  
}
      
      
unsigned long AdaptiveSparseGrid_Base::getFreeSpaceNumberInSecondTable() {
   while(secondTable[minimalEmptySecondTable]!=0) {
     ++minimalEmptySecondTable;
     if(minimalEmptySecondTable>=secondTableLength) {
        cout << " resize of second table needed! not implemented! " << endl;
         cout << "minimalEmptySecondTable" << minimalEmptySecondTable << endl;

	    assert(false);
     }
   }
   return minimalEmptySecondTable;
}
    
void AdaptiveSparseGrid_Base::setIndexInTable(const IndexDimension I, unsigned long iSetz) {
    Depth T(I);
    int t=0;
for(int d=0; d<DimensionSparseGrid;d++){
    t+=T.at(d);
}
if(max_LOne<t)max_LOne=t;


    IndexDimension Imin;
    IndexDimension Imax;

    support(I,Imin,Imax);

   for(int d=0;d<DimensionSparseGrid;++d) {
       indicesSecondTable[d + iSetz * DimensionSparseGrid] = I.getIndex(d);
       indicesSupportMin[d+iSetz*DimensionSparseGrid]=Imin.getIndex(d);
       indicesSupportMax[d+iSetz*DimensionSparseGrid]=Imax.getIndex(d);

   }


    myAssert(secondTable[iSetz] == 0);
   secondTable[iSetz] = 1;  // this means it is set, but no further link
   if(maximalOccupiedSecondTable <= iSetz) maximalOccupiedSecondTable = iSetz + 1;

    multiDepthHashGrid->addPoint(I,iSetz); // TODO is this safe to assume???
}
        
        


int AdaptiveSparseGrid_Base::getMaxDepth(int d) {
    int depth = 0;
    for (unsigned long i = 0; i < maximalOccupiedSecondTable; ++i) {
        if (secondTable[i] != 0) {
            IndexDimension I = getIndexOfTable(i);
            if (I.getDepth(d) > depth) depth = I.getDepth(d);
        }
    }
    return depth;
}




int AdaptiveSparseGrid_Base::getMaxDepth(int d, IndexDimension Index) {
    int depth = Index.getDepth(d);

    unsigned long k;
    if (depth > 0)
        while (occupied(k, Index.leftSon(d))) {
            Index = Index.leftSon(d);

            depth++;
        }


    return depth;
}


////////////////////////////////////////////////////////////////////////////////////////////
// 2. Gitter mit Rand
////////////////////////////////////////////////////////////////////////////////////////////

	
bool AdaptiveSparseGrid:: AddPoint(const IndexDimension I) {
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
        secondTable[i] = freeSpace + 2;
        isActiveNodeTable[freeSpace]=true;
        setIndexInTable(I,freeSpace);
	    return true;
      }
   }
}


bool AdaptiveSparseGrid::AddPoint(const IndexDimension I, bool hangingNode) {

   unsigned long indArray =  hash(I);
   unsigned long i = primeTable[indArray];
   if(i == 0) { // in prime table is something free;
      unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
      primeTable[indArray] = freeSpace + 1;
      setIndexInTable(I,freeSpace);
      if (hangingNode) isActiveNodeTable[freeSpace]=true;

      return true;
   }
   else { // anaylse what is going on 
      i = i-1; // shift since data are stored with shift
      if(I == getIndexOfTable(i)) return false;
      else { // search in second table;
        unsigned long iNext = secondTable[i];
        while(iNext > 1) {
	       i = iNext - 2;
	       if(I == getIndexOfTable(i)) return false;
	       iNext = secondTable[i];
	    }
	    myAssert(iNext==1);
	    unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
        secondTable[i] = freeSpace + 2;
        if (hangingNode) isActiveNodeTable[freeSpace]=true;
        setIndexInTable(I,freeSpace);

	    return true;
      }
   }
}






void AdaptiveSparseGrid::AddRecursiveSonsOfPoint(IndexDimension I, int level) {
    AddPoint(I);

    if (level > 0) {
        for (int d = 0; d < DimensionSparseGrid; ++d) {
            AddRecursiveSonsOfPoint(I.leftSon(d), level - 1);
            AddRecursiveSonsOfPoint(I.rightSon(d), level - 1);
        }
    }
}


void AdaptiveSparseGrid::AddPointsOfDepth(Depth T) {
    IndexDimension minI;
    IndexDimension maxI;

    for (int dd = 0; dd < DimensionSparseGrid; dd++) {
        minI.replace(dd, 0);
        maxI.replace(dd, 1);
    }


    for (RectangularIterator iter(minI, maxI, T); iter.goon();++iter) {
        IndexDimension I =iter.getIndex();
        AddPoint(I);
    }

}


void AdaptiveSparseGrid::addHangingNodes() {
    unsigned long maxocc = maximalOccupiedSecondTable;
//#pragma omp parallel for
    for (unsigned long i = 0; i < maxocc; ++i) {
        if (secondTable[i] != 0 && isActiveNodeTable[i]) {
            IndexDimension I = getIndexOfTable(i);
            Depth Tlocal(I);
            for (int d = 0; d < DimensionSparseGrid; ++d) {

                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {

                    IndexDimension J = I.nextFive(&mc, Tlocal);
                    unsigned long j;

                    if (!occupied(j, J))
//#pragma omp critical
                        if (AddPoint(J)) {


                            unsigned long k;
                            if (occupied(k, J)) {
                                isActiveNodeTable[k] = false;
                            }

                        }

                }


            }
        }
    }
}



int AdaptiveSparseGrid::getDOFS() {
    int dofs = 0;
    for (unsigned long k = 0; k < maximalOccupiedSecondTable; k++) {
        if (secondTable[k] != 0)
            if (isActiveNodeTable[k]) dofs++;
    }
    return dofs;
}

int AdaptiveSparseGrid::getInnerDOFS() {
    int dofs = 0;
    for (unsigned long k = 0; k < maximalOccupiedSecondTable; k++) {
            if (isActiveNodeTable[k]){
                Depth T(this->getIndexOfTable(k));
                if(T>>0)
                dofs++;
            }
    }
    return dofs;
}

int AdaptiveSparseGrid::getHangingNodes() {
    int dofs = 0;
    for (unsigned long k = 0; k < maximalOccupiedSecondTable; k++) {
        if (secondTable[k] != 0)
            if (!isActiveNodeTable[k]) dofs++;
    }
    return dofs;
}


void AdaptiveSparseGrid::completeNeumannGrid() {
/*    if(!isNeumannGrid)
        typeOfHangingNode = new int[secondTableLength];

    isNeumannGrid = true;



    for (unsigned long i = 0; i < secondTableLength; ++i) typeOfHangingNode[i]=0;*/

    completeGrid();

/*
    CompleteToLocalTensorProductGrid();
    completeGrid();


    // completeRandNah();


    addHangingNodes();*/


}

void AdaptiveSparseGrid::completeGrid() {
    unsigned long k;
    for (bool changed = true; changed;) {
        changed = false;
        unsigned long end = maximalOccupiedSecondTable;
        for (unsigned long i = 0; i < end; ++i) {
            if (secondTable[i] != 0 && isActiveNodeTable[i]) {
                IndexDimension I = getIndexOfTable(i);
                for (int d = 0; d < DimensionSparseGrid; ++d) {
                    if (I.getDepth(d) > 0) {
                        if (AddPoint(I.nextLeft(d),true)) {
                            changed = true;

                        } else {
                            if (occupied(k, I.nextLeft(d)))
                                if (!isActiveNodeTable[k]) {
                                    isActiveNodeTable[k] = true;
                                }
                        }
                        if (AddPoint(I.nextRight(d),true)) changed = true;
                        else {
                            if (occupied(k, I.nextRight(d)))
                                if (!isActiveNodeTable[k]) {
                                    isActiveNodeTable[k] = true;
                                }
                        }
                    }
                    if (I.getDepth(d) == 0) {
                        IndexDimension J = I;
                        if (I.getIndex(d) == 0)
                            J.replace(d, 1);

                        if (AddPoint(J,true)) {
                            changed = true;
                        } else {
                            if (occupied(k, J))
                                if (!isActiveNodeTable[k]) {
                                    isActiveNodeTable[k] = true;
                                }
                        }

                        J = I;
                        if (I.getIndex(d) == 1)
                            J.replace(d, 0);

                        if (AddPoint(J,true)) {
                            changed = true;
                        } else {
                            if (occupied(k, J))
                                if (!isActiveNodeTable[k]) {
                                    isActiveNodeTable[k] = true;
                                }
                        }
                    }
                }
            }
       }
   }
   //updateMaxDepth();
}

void AdaptiveSparseGrid::completeGridWithoutBoundary() {
unsigned long k;
    for (bool changed = true; changed;) {
        changed = false;
        unsigned long end=maximalOccupiedSecondTable;
        for (unsigned long i = 0; i < end; ++i) {
            if (secondTable[i] != 0) {
                IndexDimension I = getIndexOfTable(i);
                for (int d = 0; d < DimensionSparseGrid; ++d) {
                    IndexDimension J = I.nextLeft(d);
                    if (J.isNotAtBoundary()&& !occupied(k,J)) {

                        if (AddPoint(I.nextLeft(d))) {

                            changed = true;
                        }
                    }
                    J = I.nextRight(d);

                    if (J.isNotAtBoundary()&&!occupied(k,J))
                        if (AddPoint(I.nextRight(d))){

                            changed = true;
                        }
                }
            }
        }
    }
    //updateMaxDepth();
}

void AdaptiveSparseGrid::completeDirichletGrid() {


    completeGridWithoutBoundary();

   CompleteToLocalTensorProductGrid();

       completeGridWithoutBoundary();


       addHangingNodes();

       addPoints2(*this);



}

void AdaptiveSparseGrid::completeDirichletGrid_NEU() {


    ListOfDepthOrderedSubgrids listOfDepthOrderedSubgrids(*this);


    int rank;


    completeGridWithoutBoundary();




    CompleteToLocalTensorProductGrid();





    completeGridWithoutBoundary();




    addHangingNodes();




    DepthList list2(*this);



    addPoints2_NEU(*this,list2);





}

void AdaptiveSparseGrid::copy(AdaptiveSparseGrid& newgrid) {




    for(unsigned long i=0;i<primeTableLength;++i) primeTable[i] = 0;

    for(unsigned long i=0;i<secondTableLength;++i) secondTable[i] = 0;


    for (unsigned long i = 0; i < secondTableLength; ++i) isActiveNodeTable[i] = false;





    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSecondTable[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMin[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMax[i] = 0;


    minimalEmptySecondTable = 0;
    maximalOccupiedSecondTable = 0;

    delete multiDepthHashGrid;

    multiDepthHashGrid = new MultiDepthHashGrid(*this);


    for(unsigned long k=0; k<newgrid.getMaximalOccupiedSecondTable(); k++){
        IndexDimension I = newgrid.getIndexOfTable(k);
        if(newgrid.getActiveTable()[k]){
                        AddPoint(I,true);        }else{
            AddPoint(I,false);
        }
    }

}


void AdaptiveSparseGrid::copy_inner(AdaptiveSparseGrid& newgrid) {




    for(unsigned long i=0;i<primeTableLength;++i) primeTable[i] = 0;

    for(unsigned long i=0;i<secondTableLength;++i) secondTable[i] = 0;


    for (unsigned long i = 0; i < secondTableLength; ++i) isActiveNodeTable[i] = false;





    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSecondTable[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMin[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMax[i] = 0;

    minimalEmptySecondTable = 0;
    maximalOccupiedSecondTable = 0;

    delete multiDepthHashGrid;

    multiDepthHashGrid = new MultiDepthHashGrid(*this);


    for(unsigned long k=0; k<newgrid.getMaximalOccupiedSecondTable(); k++){
        if(newgrid.getActiveTable()[k]){
            IndexDimension I = newgrid.getIndexOfTable(k);
            Depth T(I);
            if(T>>0) AddPoint(I);

        }
    }

}

void AdaptiveSparseGrid::add_outer(AdaptiveSparseGrid& newgrid) {

    for(unsigned long k=0; k<newgrid.getMaximalOccupiedSecondTable(); k++){
        if(newgrid.getActiveTable()[k]){
            IndexDimension I = newgrid.getIndexOfTable(k);
            Depth T(I);
            if(!(T>>0)) AddPoint(I);

        }
    }

}


void AdaptiveSparseGrid::add_RecursiveSonsBoundary(AdaptiveSparseGrid &newgrid) {


}

void AdaptiveSparseGrid::clear(){
    for(unsigned long i=0;i<primeTableLength;++i) primeTable[i] = 0;

    for(unsigned long i=0;i<secondTableLength;++i) secondTable[i] = 0;


    for (unsigned long i = 0; i < secondTableLength; ++i) isActiveNodeTable[i] = false;





    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSecondTable[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMin[i] = 0;
    for (unsigned long i = 0; i < secondTableLength * DimensionSparseGrid; ++i) indicesSupportMax[i] = 0;


    minimalEmptySecondTable = 0;
    maximalOccupiedSecondTable = 0;

    delete multiDepthHashGrid;

    multiDepthHashGrid = new MultiDepthHashGrid(*this);


}

