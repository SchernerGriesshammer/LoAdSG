//
// Created by to35jepo on 5/8/23.
//

#ifndef RUN_DEPTHITERATOR_H
#define RUN_DEPTHITERATOR_H


#include "../sgrid/sparseGrid.h"

class DepthList {

public:
    DepthList(AdaptiveSparseGrid_Base& grid){

        for(unsigned long k=0; k<grid.getMaximalOccupiedSecondTable(); k++){
            IndexDimension I = grid.getIndexOfTable(k);
            Depth T(I);


            bool addT = true;
            for (auto it = mylist.begin(); it != mylist.end(); ++it) {
                if(*it==T){
                    addT = false;
                    break;
                }
            }
            //auto it = std::find(mylist.begin(), mylist.end(), T);

            if (addT) {
                mylist.push_back(T);
            }
        }

        int ret=0;


        for (auto it = mylist.begin(); it != mylist.end(); ++it) {
            int sum = 0;
            for(int d=0;d<DimensionSparseGrid;d++){
                sum+=it->at(d);
                if(MaxDepth.at(d)<it->at(d))
                    MaxDepth.set(it->at(d),d);
            }
            if(ret<sum)ret=sum;

        }

        maxLOne=ret;
    }

    bool isIncluded(Depth T){
        auto it = std::find(mylist.begin(), mylist.end(), T);
        if (it != mylist.end()) {
            return true;
        }
        return false;
    }


    int getMaxLOne() {
        return maxLOne;
    }

    Depth getMaxDepth(){
        return MaxDepth;
    }

    std::list<Depth>::iterator begin_all(){
        return mylist.begin();
    }

    std::list<Depth>::iterator end_all(){
        return mylist.end();
    }

    std::list<Depth>::iterator begin_sorted(){
        return sortedlist.begin();
    }

    std::list<Depth>::iterator end_sorted(){
        return sortedlist.end();
    }


    void recursiveDepth(bool *Restrictions, bool checkProlongation, int d, Depth *T);

    void recursiveDepthNeumann(bool *Restrictions, bool checkProlongation, int d, Depth *T);

    void sortDepths(bool *Restrictions);

    void sortDepthsNeumann(bool *Restrictions);

    void sortFineToCoarse();


    std::list<Depth> *getAllDepths() { return &mylist; }
    std::list<Depth> *getSortierteTiefen() { return &sortedlist; }

private:
    std::list<Depth> mylist;
    std::list<Depth> sortedlist;
    int  maxLOne;
    Depth MaxDepth;

};


#endif //RUN_DEPTHITERATOR_H
