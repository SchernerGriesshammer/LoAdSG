  
/**********************************************************************************
 * Copyright 2015 Christoph Pflaum 
 *              Department Informatik Lehrstuhl 10 - Systemsimulation
 *              Friedrich-Alexander Universität Erlangen-Nürnberg
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **********************************************************************************/

#include "../abbrevi.h"
#include "../myAssert.h"
#include "../indices/index.h"
#include "../primes/prime.h"
#include "sparseGrid.h"
#include "depth.h"
#include "ListOfDepthOrderedGrids.h"
#include <iostream>
#include <omp.h>

#include "localTensorProduct.h"

using namespace std;


/////////////////////////////////////////////////////////////////////////////////
// Implementierung der Methoden
/////////////////////////////////////////////////////////////////////////////////


LocalForTensor::LocalForTensor(AdaptiveSparseGrid *grid_) : T(0) {
    grid = grid_;
    secondTableLength = grid->getLengthSecondTable();
    marker = new Marker[secondTableLength];
   for(unsigned long i=0;i<secondTableLength;++i)
       marker[i] = nicht;     
}

LocalForTensor::~LocalForTensor() {
   delete[] marker; 
}


unsigned long LocalForTensor::findAnicht(unsigned long iNow) {
    unsigned long iMax = grid->getMaximalOccupiedSecondTable();
    for (unsigned long i = iNow + 1; i < iMax; ++i) {
        if (marker[i] == nicht) return i;
    }
    for (unsigned long i = 0; i < iMax; ++i) {
        if (marker[i] == nicht) return i;
    }
    return iMax; // also nichts gefunden
}

unsigned long LocalForTensor::findAnichtWithDepth(unsigned long iNow, Depth &T_) {
    unsigned long iMax = grid->getMaximalOccupiedSecondTable();
    for (unsigned long i = iNow + 1; i < iMax; ++i) {
        IndexDimension Index = grid->getIndexOfTable(i);
        if (marker[i] == nicht && Depth(Index) == T_ && grid->getActiveTable()[i]) {
            return i;
        }
    }
/*        for(unsigned long i=0;i<iMax;++i) {
            if(marker[i]==nicht) return i;
        }*/
    return iMax; // also nichts gefunden
}


IndexDimension LocalForTensor::StartSearchComponent(unsigned long iNow) {
    maxIndex = grid->getIndexOfTable(iNow);
    minIndex = maxIndex;
    T = Depth(maxIndex);
    return minIndex;
}

void LocalForTensor::recursiveMarkArbeite(IndexDimension indexNow, unsigned long iNow) {
    marker[iNow] = arbeite;

    //Rekursion
    for (MultiDimCompass mc; mc.goon(); ++mc) {
        IndexDimension J = indexNow.nextTwoStep(&mc, T);
        unsigned long indexOfData;
        if (grid->occupied(indexOfData, J)) {
            Marker m = marker[indexOfData];
            if (m != arbeite)
                recursiveMarkArbeite(J, indexOfData);
        }
    }
    //max min erneuern -> wichtig: Tiefe ändert sich dadurch nicht
    maxIndex = IndexDimension::Maximum(maxIndex, indexNow);
    minIndex = IndexDimension::Minimum(minIndex, indexNow);
}


bool LocalForTensor::markArbeiteAndAddRechteck() {

    bool touch = false;

    for (RectangularIterator iter(minIndex, maxIndex); iter.goon(); ++iter) {
        unsigned long indexOfData;
        IndexDimension test = iter.getIndex();
        bool exists = grid->occupied(indexOfData, test);
        if (exists) {
            if (marker[indexOfData] == fertig) touch = true;
        } else {
            grid->AddPoint(iter.getIndex());
            bool noError = grid->occupied(indexOfData, iter.getIndex());
            if (noError == false) {
                std::cout << " Error in void LocalForTensor::markRechteckArbeite()" << std::endl;
            }
        }
        marker[indexOfData] = arbeite;
    }

    return touch;
}


void LocalForTensor::markRechteck(Marker mm) {
    for (RectangularIterator iter(minIndex, maxIndex); iter.goon(); ++iter) {
        unsigned long indexOfData;
        bool noError = grid->occupied(indexOfData, iter.getIndex());
        if (noError == false) {
            std::cout << " Error in void LocalForTensor::markRechteckFertig()" << std::endl;
        }
        marker[indexOfData] = mm;
    }
}
///////////////////////////////////////////////////////////////////////////////////

void LocalForTensor::recursiveMarkArbeite2(IndexDimension indexNow, unsigned long iNow) {
    marker[iNow] = arbeite;

    //Rekursion
    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
        IndexDimension J = indexNow.nextFive(&mc, T);


        unsigned long indexOfData;
        if (grid->occupied(indexOfData, J)) {
            if (grid->getActiveTable()[indexOfData]) {

                maxIndex = IndexDimension::Maximum(maxIndex, J);
                minIndex = IndexDimension::Minimum(minIndex, J);
            }

        }


    }


}


bool LocalForTensor::checkRechteck(IndexDimension indexNow) {


    bool changed = false;


    for (RectangularIterator iter(minIndex, maxIndex, T); iter.goon(); ++iter) {
        unsigned long indexOfData;
        IndexDimension test = iter.getIndex();
        bool exists = grid->occupied(indexOfData, test);
        if (!exists) {
            grid->AddPoint(iter.getIndex(), true);
            changed = true;


            bool noError = grid->occupied(indexOfData, iter.getIndex());
            if (noError == false) {
                std::cout << " Error in void LocalForTensor::markRechteckArbeite()" << std::endl;
            }

        } else {
            grid->getActiveTable()[indexOfData] = true;
        }

    }
    return changed;

}

//////////////////////////////////////////////////////////////////////////////////////////





void AdaptiveSparseGrid::CompleteToLocalTensorProductGrid() {

    /* // Grid hat immer zuerst hierarchische Struktur und dann noch als TP erweitert
     completeGridWithoutBoundary();*/


    //1. Start
    LocalForTensor constructorTensor(this);


    ListOfDepthOrderedSubgrids orderedSubgrids(*this);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    outerIterationsSubgrids.gotoBegin();



    //Gehe alle Tiefen von Coarse to Fine durch
    do {

        Depth T = outerIterationsSubgrids.getSubgrid()->getT();


        unsigned long iNow = constructorTensor.findAnichtWithDepth(-1, T);
        while (iNow < maximalOccupiedSecondTable) {


            IndexDimension indexNow = constructorTensor.StartSearchComponent(iNow);


            // Zusammenhangskomponente wird untersucht und bereits vorhandene Punkte auf arbeite gesetzt
            constructorTensor.recursiveMarkArbeite2(indexNow, iNow);

            constructorTensor.checkRechteck(indexNow);


            constructorTensor.checkStorage();
            iNow = constructorTensor.findAnichtWithDepth(iNow, T);

        }


    } while (outerIterationsSubgrids.next());


}


void AdaptiveSparseGrid::CompleteToLocalTensorProductGrid_NEU(ListOfDepthOrderedSubgrids* listOfDepthOrderedSubgrids) {

    /* // Grid hat immer zuerst hierarchische Struktur und dann noch als TP erweitert
     completeGridWithoutBoundary();*/


    //1. Start
    LocalForTensor constructorTensor(this);
    //ListOfDepthOrderedSubgrids orderedSubgrids(*this);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(*listOfDepthOrderedSubgrids);
    outerIterationsSubgrids.gotoBegin();





    //Gehe alle Tiefen von Coarse to Fine durch
    do {

        Depth T = outerIterationsSubgrids.getSubgrid()->getT();
        vector<pair<IndexDimension,IndexDimension>> minmaxvector;

        unsigned long iNow = constructorTensor.findAnichtWithDepth(-1, T);
        while (iNow < maximalOccupiedSecondTable) {


            IndexDimension indexNow = constructorTensor.StartSearchComponent(iNow);


            // Zusammenhangskomponente wird untersucht und bereits vorhandene Punkte auf arbeite gesetzt
            constructorTensor.recursiveMarkArbeite2(indexNow, iNow);

            //constructorTensor.checkRechteck(indexNow);
            minmaxvector.emplace_back(constructorTensor.getMin(),constructorTensor.getMax());



            constructorTensor.checkStorage();
            iNow = constructorTensor.findAnichtWithDepth(iNow, T);

        }
        if (!minmaxvector.empty()) {
            int end_v = minmaxvector.size();

//#pragma omp parallel
            {
//#pragma omp for
                for (int i = 0; i < end_v; i++) {

                    pair<IndexDimension, IndexDimension> minmax;

//#pragma omp critical
                    {
                        minmax = minmaxvector.at(i);  // Accessing shared data within a critical section
                    }
                    checkRechteck(minmax.first, minmax.second, T, this);

                }
            }
        }


    } while (outerIterationsSubgrids.next());





}


bool checkRechteck(IndexDimension minIndex,IndexDimension maxIndex, Depth T, AdaptiveSparseGrid* grid){
    bool changed = false;


    for (RectangularIterator iter(minIndex, maxIndex, T); iter.goon(); ++iter) {

        unsigned long indexOfData;
        IndexDimension test = iter.getIndex();
        bool exists;
#pragma omp critical
        {
            exists = grid->occupied(indexOfData, test);

            if (!exists) {

                grid->AddPoint(iter.getIndex(), true);
                changed = true;


                bool noError = grid->occupied(indexOfData, iter.getIndex());
                if (noError == false) {
                    std::cout << " Error in void LocalForTensor::markRechteckArbeite()" << std::endl;
                }

            } else {
                grid->getActiveTable()[indexOfData] = true;
            }


        };


    }

    return changed;
};