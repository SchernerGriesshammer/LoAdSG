/****************************
 * Test for ordered grids
 ***************************/

#include "../abbrevi.h"
#include "../myAssert.h"
#include "../mympi.h"
#include "../indices/index.h"
#include "sparseGrid.h"
#include "depth.h"
#include "ListOfDepthOrderedGrids.h"
#include "../tests/testing.h"


SubgridFixedDepth::SubgridFixedDepth(Depth &T_, IndexDimension p, unsigned long i) : T(T_) {
    // warum dynamische liste? liste ist automatisch dynamisch
    //listPoints =  new std::list<SubgridFixedDepth::IndexPlusi>;


    listPoints.clear();


    listPoints.emplace_back(SubgridFixedDepth::IndexPlusi(p, i));

}

SubgridFixedDepth::SubgridFixedDepth(Depth &T_) : T(T_) {
    //listPoints = std::list<SubgridFixedDepth::IndexPlusi>;
    listPoints.clear();
}

bool SubgridFixedDepth::addPoint(IndexDimension &p, unsigned long i) {
    //for (list<SubgridFixedDepth::IndexPlusi>::iterator it = listPoints.begin(); it != listPoints.end(); it++)
    for (auto it:listPoints) {
        SubgridFixedDepth::IndexPlusi vergleich = it;
        SubgridFixedDepth::IndexPlusi vergleich2(p, i);
        if (vergleich == vergleich2) return false;
    }
    SubgridFixedDepth::IndexPlusi index(p, i);

    listPoints.push_back(index);
    return true;
}

bool SubgridFixedDepth::isOccupied(IndexDimension &p, unsigned long i) {
    // for (list<SubgridFixedDepth::IndexPlusi>::iterator it = listPoints.begin(); it != listPoints.end(); it++)
    for (auto it:listPoints) {
        SubgridFixedDepth::IndexPlusi vergleich = it;
        SubgridFixedDepth::IndexPlusi vergleich2(p, i);
        if (vergleich == vergleich2) return true;
    }
    return false;
}

SubgridFixedDepth::iterator::iterator(SubgridFixedDepth &dataStructure_) :
        dataStructure(&dataStructure_) {

    it = dataStructure->listPoints.begin();
}

void SubgridFixedDepth::iterator::gotobegin() {
    it = dataStructure->listPoints.begin();
}

bool SubgridFixedDepth::iterator::next() {
    ++it;
    return it != dataStructure->listPoints.end();
}

bool SubgridFixedDepth::iterator::empty() {
    return it == dataStructure->listPoints.end();
}


///////////////////////////////////////////////////////////////////////////////////////////////////

ListOfDepthOrderedSubgrids::ListOfDepthOrderedSubgrids(AdaptiveSparseGrid_Base &grid) {

    TiefeVector.resize(0);

    for (int d = 0; d < DimensionSparseGrid; d++) {
        MaximalDepth[d] = grid.getMaxDepth(d);
    }


    unsigned long endIndex = grid.maximalOccupiedSecondTable;
    //double* dataTable         = sparseGrid->dataTable;
    dataInteger *secondTable = grid.secondTable;
    //unsigned int numberOfData = sparseGrid->numberOfData;

    for (unsigned long i = 0; i < endIndex; ++i) {
        IndexDimension p = grid.getIndexOfTable(i);
        if (secondTable[i] != 0)
            addPoint(p, i);
    }

    //Entferne leere Subgrids


    for (auto it = TiefeVector.begin(); it != TiefeVector.end(); ++it) {
        std::list<SubgridFixedDepth> &listOfSubgrids = *it;

        if (listOfSubgrids.empty()) {
            TiefeVector.erase(it);
            it--;
        }
    }


};

ListOfDepthOrderedSubgrids::ListOfDepthOrderedSubgrids() {
  TiefeVector.resize(0);

};


void ListOfDepthOrderedSubgrids::addPoint(IndexDimension& p, unsigned long i) {
    Depth T(p);

    unsigned int t = T.LoneNorm();

    if (TiefeVector.size() <= t) TiefeVector.resize(t + 1);


    std::list<SubgridFixedDepth> &listOfSubgrids = TiefeVector.at(t);

    for (auto it = listOfSubgrids.begin();
         it != listOfSubgrids.end();
         ++it) {
        if (it->getT() == T) {
            it->addPoint(p, i);
            return;
        }
    }


    listOfSubgrids.emplace_back(SubgridFixedDepth(T, p, i));
}

bool ListOfDepthOrderedSubgrids::isIncluded(Depth &T) {
    iterator iter(*this);
    iter.gotoBegin();
    do {
        Depth Tlocal = iter.getSubgrid()->getT();
        if (Tlocal == T) return true;
    } while (iter.next());
    return false;
}


void ListOfDepthOrderedSubgrids::recursiveDepthCoarse(int d, Depth *T) {

    if (d >= 0) {
        for (int t = MaximalDepth[d]; t >= 0; t--) {
            T->set(t, d);
            recursiveDepthCoarse(d - 1, T);
        }
    }

    if (d < 0) {
        if (isIncluded(*T))
            SortierteTiefen.push_back(*T);
    }
}


ListOfDepthOrderedSubgrids::ListOfDepthOrderedSubgrids(AdaptiveSparseGrid_Base &grid, bool *restrictions) {

    TiefeVector.resize(0);

    unsigned long endIndex = grid.maximalOccupiedSecondTable;
    //double* dataTable         = sparseGrid->dataTable;
    dataInteger *secondTable = grid.secondTable;
    //unsigned int numberOfData = sparseGrid->numberOfData;

    for (unsigned long i = 0; i < endIndex; ++i) {
        if (secondTable[i] != 0) {
            IndexDimension p = grid.getIndexOfTable(i);
            addPoint(p, i);
        }
    }
    for (int d = 0; d < DimensionSparseGrid; d++) {
        MaximalDepth[d] = grid.getMaxDepth(d);
    }

    sortDepths(restrictions);
}


ListOfDepthOrderedSubgrids::iterator::iterator(ListOfDepthOrderedSubgrids &dataStructure_) :
        dataStructure(&dataStructure_) {
    t = int(dataStructure->TiefeVector.size() - 1);
    it = std::list<SubgridFixedDepth>::iterator(dataStructure->TiefeVector.at(t).begin());
}


void ListOfDepthOrderedSubgrids::iterator::gotoBegin() {
    t = 0;
   it = std::list<SubgridFixedDepth>::iterator(dataStructure->TiefeVector.at(t).begin());

}

void ListOfDepthOrderedSubgrids::iterator::gotoEnd() {
    t = int(dataStructure->TiefeVector.size() - 1);
    it = std::list<SubgridFixedDepth>::iterator(dataStructure->TiefeVector.at(t).end());
    --it;

    /*for(int i=t; i >=0 ;i--) {
        std::list<SubgridFixedDepth>::iterator testbegin(dataStructure->TiefeVector.at(i).begin());
        std::list<SubgridFixedDepth>::iterator testend(dataStructure->TiefeVector.at(i).end());
        std::list<SubgridFixedDepth>::iterator it2 = testend;
        while (it2 !=testbegin)
        {
            --it2;
            *//*do stuff *//*
          cout << it2->getT().at(0) << it2->getT().at(1) << "  ";


      }*/

    //}


}

bool ListOfDepthOrderedSubgrids::iterator::empty() {
   return (t<0 || t>=dataStructure->TiefeVector.size());
}

bool ListOfDepthOrderedSubgrids::iterator::previous() {
    if (empty()) return false;

    --it;
    if (it == --dataStructure->TiefeVector.at(t).begin()) {
        --t;
        if (t >= 0) {
            it = std::list<SubgridFixedDepth>::iterator(dataStructure->TiefeVector.at(t).end());
            --it;
        } else {
            return false;
        }
    }
    return true;
}



bool ListOfDepthOrderedSubgrids::iterator::next() {
   if(empty()) return false;

   ++it;
   if(it==dataStructure->TiefeVector.at(t).end()) {
      ++t;
      if(t<dataStructure->TiefeVector.size()) {
	 it = std::list<SubgridFixedDepth>::iterator(dataStructure->TiefeVector.at(t).begin());
      }
      else {
	 return false;
      }
   }
   return true;
}




