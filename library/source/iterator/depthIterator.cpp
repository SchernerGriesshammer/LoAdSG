//
// Created by to35jepo on 5/8/23.
//

#include "depthIterator.h"

 void DepthList::recursiveDepth(bool *Restrictions, bool checkProlongation, int d, Depth *T) {
    if (checkProlongation && d > 0) {
        if (!Restrictions[d]) {
            for (int t = 1; t <= MaxDepth.at(d); t++) {
                T->set(t, d);
                recursiveDepth(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (checkProlongation && d == 0) {
        if (!Restrictions[d]) {
            for (int t = 1; t <= MaxDepth.at(d); t++) {
                T->set(t, d);
                recursiveDepth(Restrictions, false, DimensionSparseGrid - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, false, DimensionSparseGrid - 1, T);
        }
    }

    if ((!checkProlongation) && d >= 0) {
        if (Restrictions[d]) {
            for (int t = MaxDepth.at(d); t > 1; t--) {
                T->set(t, d);
                recursiveDepth(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (d < 0 && isIncluded(*T)) {
        sortedlist.push_back(*T);
      }

}

void DepthList::recursiveDepthNeumann(bool *Restrictions, bool checkProlongation, int d, Depth *T) {
    if (checkProlongation && d > 0) {
        if (!Restrictions[d]) {
            for (int t = 0; t <= MaxDepth.at(d); t++) {
                T->set(t, d);
                recursiveDepthNeumann(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepthNeumann(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (checkProlongation && d == 0) {
        if (!Restrictions[d]) {
            for (int t = 0; t <= MaxDepth.at(d); t++) {
                T->set(t, d);
                recursiveDepthNeumann(Restrictions, false, DimensionSparseGrid - 1, T);
            }
        } else {
            recursiveDepthNeumann(Restrictions, false, DimensionSparseGrid - 1, T);
        }
    }

    if ((!checkProlongation) && d >= 0) {
        if (Restrictions[d]) {
            for (int t = MaxDepth.at(d); t > 0; t--) {
                T->set(t, d);
                recursiveDepthNeumann(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepthNeumann(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (d < 0 && isIncluded(*T)) {
        sortedlist.push_back(*T);
    }

}

void  DepthList::sortDepths(bool *Restrictions) {

    sortedlist.clear();
    Depth T(0);
    recursiveDepth(Restrictions, true, DimensionSparseGrid - 1, &T);

}

void  DepthList::sortDepthsNeumann(bool *Restrictions) {

    sortedlist.clear();
    Depth T(0);
    recursiveDepthNeumann(Restrictions, true, DimensionSparseGrid - 1, &T);

}


void  DepthList::sortFineToCoarse() {

    sortedlist.clear();
    sortedlist = mylist;

   sortedlist.sort([](Depth a, Depth b){
        int la=0;
        int lb=0;
        for(int d=0; d<DimensionSparseGrid; d++){
            la+=a.at(d);
            lb+=b.at(d);
        }
        return la > lb; // Custom comparison for descending order
    });



}