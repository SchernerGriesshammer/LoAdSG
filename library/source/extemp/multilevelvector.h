//
// Created by scherner on 31.05.21.
//

#ifndef GRUN_MULTILEVELVECTOR_H
#define GRUN_MULTILEVELVECTOR_H

#include "vector.h"
#include "../sgrid/multilevelSparseGrid.h"


class GetNextSmallerDepthsIterator{
    private:
    Depth startDepth;
    Depth lastDepth;
    public:
    GetNextSmallerDepthsIterator(const Depth& T):startDepth(T),lastDepth(T){
        // next();
    }
    Depth operator*() const{
        return lastDepth;
    };
    bool next(){
        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            if(lastDepth.at(i)==0)lastDepth.set(startDepth.at(i),i);
            else {
                lastDepth.set(lastDepth.at(i)-1,i);
                return true;
            }
        }
        return false;
    }
    void Print(){
        do{
            lastDepth.Print();
        }while (next());
    }
};

class GetNextSmallerDepthsIteratorInner{
private:
    Depth startDepth;
    Depth lastDepth;
public:
    GetNextSmallerDepthsIteratorInner(const Depth& T):startDepth(T),lastDepth(T){
        // next();
    }
    Depth operator*() const{
        return lastDepth;
    };
    bool next(){
        for (size_t i = 0; i < DimensionSparseGrid; i++)
        {
            if(lastDepth.at(i)==1)lastDepth.set(startDepth.at(i),i);
            else {
                lastDepth.set(lastDepth.at(i)-1,i);
                return true;
            }
        }
        return false;
    }
    void Print(){
        do{
            lastDepth.Print();
        }while (next());
    }
};



class MultiLevelVector {


public:
    MultiLevelVector(MultiLevelAdaptiveSparseGrid &grid);

    ~MultiLevelVector();


    double getValue(unsigned long i) {
        return dataTableVector[i];
    }

    void setValue(unsigned long i, double val) {
        dataTableVector[i]=val;
    }


    inline void setMultiLevelValues(VectorSparseG &vector, Depth &T);

    inline void setMultiLevelValues2(VectorSparseG &vector, Depth &T);

    inline void setMultiLevelValuesOMP(VectorSparseG &vector, Depth &T, MultiDepthHashGrid& multiDepthHashGrid);

    inline void setMultiLevelValues(VectorSparseG &vector, Depth &T, ListOfDepthOrderedSubgrids& list);

    inline void addMultiLevelValues(VectorSparseG &vector, Depth &T);

    inline void setMultiLevelValuesInVector(VectorSparseG &vector, Depth &T, ListOfDepthOrderedSubgrids &list);

    MultiLevelAdaptiveSparseGrid *getSparseGrid() { return sparseGrid; };


    ///// Operatoren

    inline void operator=(double x) {
        for (unsigned long i = 0; i < sparseGrid->getMaximalOccupiedSecondTable(); i++)
            dataTableVector[i] = x;
    };

    void operator+=(MultiLevelVector &vector) {
        for (unsigned long i = 0; i < sparseGrid->getMaximalOccupiedSecondTable(); i++) {
            dataTableVector[i] = vector.getValue(i);
        }
    }



    ///// MPI

    void Broadcast(int rank);

    bool mpi_doit();

    void sendTo(int torank);

    void ReduceSum(int rank);


    ///// PrintFunktionen
    void PrintDoubleTwoD(int level);

    void PrintDoubleTwoD(int level, Depth T);

    ///// weiteres
    bool workonindex(unsigned long i);

private:
    MultiLevelAdaptiveSparseGrid *sparseGrid;
    double *dataTableVector;


};

inline void MultiLevelVector::setMultiLevelValues(VectorSparseG &vector, Depth &T) {

    double *data2 = vector.getDatatableVector();
    AdaptiveSparseGrid_Base *grid = vector.getSparseGrid();
    // cout <<endl;
    // cout <<  grid->getMaximalOccupiedSecondTable() << endl;
    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;
        // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
        for (size_t i = 0; i < mapping.size(); i++)
        {
            unsigned long k;
            IndexDimension I = depthGrid._map.getIndexOfTable(i);
            if (sparseGrid->occupied(k, I, T))
                dataTableVector[k] = data2[mapping[i]];
        }
        // cout <<endl;

    }while(iter.next());


};

inline void MultiLevelVector::setMultiLevelValues2(VectorSparseG &vector, Depth &T) {

double *data2 = vector.getDatatableVector();
    AdaptiveSparseGrid_Base *grid = vector.getSparseGrid();
    // cout <<endl;
    // cout <<  grid->getMaximalOccupiedSecondTable() << endl;
    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;

        if(depthGrid.getNumberOfEntries()>0) {
            int end = mapping.size();
            for (size_t i = 0; i < end; i++) {
                unsigned long k;
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (sparseGrid->occupied(k, I, T))
                    dataTableVector[k] = data2[mapping[i]];
            }
        }
            // cout <<endl;

    }while(iter.next());

};


inline void MultiLevelVector::setMultiLevelValuesOMP(VectorSparseG &vector, Depth &T, MultiDepthHashGrid& multiDepthHashGrid) {

    double *data2 = vector.getDatatableVector();

    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = multiDepthHashGrid.getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;

        if(depthGrid.getNumberOfEntries()>0) {
            auto end = mapping.size();
            for (size_t i = 0; i < end; i++) {
                unsigned long k;
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (sparseGrid->occupied(k, I, T))
                    dataTableVector[k] = data2[mapping[i]];
            }
        }


    }while(iter.next());

};


inline void MultiLevelVector::setMultiLevelValues(VectorSparseG &vector, Depth &T, ListOfDepthOrderedSubgrids &list) {
    unsigned long k;
    double *data2 = vector.getDatatableVector();
    ListOfDepthOrderedSubgrids::iterator iter(list);
    do {
        Depth T_new = iter.getDepth();
        if(T_new <= T){
            SubgridFixedDepth::iterator iter_inner(*list.getSubgrid(T_new));
            do {
                unsigned long i = iter_inner.geti();
                IndexDimension I = iter_inner.getPoint();
                if (sparseGrid->occupied(k, I, T))
                    dataTableVector[k] = data2[i];

            } while (iter_inner.next());

        }


    } while (iter.next());

   /* double *data2 = vector.getDatatableVector();
    unsigned long k;
    AdaptiveSparseGrid_Base *grid = vector.getSparseGrid();
    unsigned long maxocc = grid->getMaximalOccupiedSecondTable();
    unsigned long *secondTable = grid->getSecondTable();

    for (unsigned long i = 0; i < maxocc; i++) {
        if (secondTable[i] != 0) {
            IndexDimension I = grid->getIndexOfTable(i);
            Depth Tlocal(I);
            if (Tlocal <= T) {
                if (sparseGrid->occupied(k, I, T))
                    dataTableVector[k] = data2[i];
            }
        }
    }*/
};


inline void MultiLevelVector::addMultiLevelValues(VectorSparseG &vector, Depth &T) {


    double *data2 = vector.getDatatableVector();
    AdaptiveSparseGrid_Base *grid = vector.getSparseGrid();
    // cout <<endl;
    // cout <<  grid->getMaximalOccupiedSecondTable() << endl;

    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;

        SingleDepthHashGrid& depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;
        // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
        if(depthGrid.getNumberOfEntries()>0) {
            for (size_t i = 0; i < mapping.size(); i++) {
                unsigned long k;
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (sparseGrid->occupied(k, I, T))
                    dataTableVector[k] = dataTableVector[k] + data2[mapping[i]];
            }
        }
        // cout <<endl;

    }while(iter.next());






};

void MultiLevelVector::setMultiLevelValuesInVector(VectorSparseG &vector, Depth &T, ListOfDepthOrderedSubgrids &list) {
    ListOfDepthOrderedSubgrids::iterator iter(list);
    iter.gotoBegin();
    do {
        Depth Tlocal = iter.getDepth();
        if (Tlocal <= T) {
            SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
            inneriter.gotobegin();
            do {
                IndexDimension Index = inneriter.getPoint();
                unsigned long k;
                if (sparseGrid->occupied(k, Index, T))
                    vector.setValue(Index, dataTableVector[k]);

            } while (inneriter.next());
        }
    } while (iter.next());

/*    double *data2 = vector.getDatatableVector();
    unsigned long k;
    AdaptiveSparseGrid_Base* grid = vector.getSparseGrid();
    unsigned long maxocc =  grid->getMaximalOccupiedSecondTable();
    unsigned long* secondTable = grid->getSecondTable();

    for (unsigned long i = 0; i < maxocc; i++) {
        if(secondTable[i]!=0) {
            IndexDimension I = grid->getIndexOfTable(i);
            Depth Tlocal(I);
            if (Tlocal <= T) {
                if (sparseGrid->occupied(k, I, T))
                    data2[i]=dataTableVector[k];
            }
        }
    }*/
}

#endif //GRUN_MULTILEVELVECTOR_H
