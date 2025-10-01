//
// Created by scherner on 21.09.21.
//

#include "multilevelSparseGrid.h"
#include "sparseGrid.h"


void MultiLevelAdaptiveSparseGrid::createbysparsegrid(AdaptiveSparseGrid_Base *grid) {

    for (unsigned long k = 0; k < grid->getMaximalOccupiedSecondTable(); k++) {
        IndexDimension I = grid->getIndexOfTable(k);
        Depth T(I);

        AddPointDepth(I, T);

        bool found = (std::find(liste.begin(), liste.end(), T) != liste.end());
        if (!found)
            liste.push_back(T);

    }


    for (Depth T: liste) {

        for (unsigned long k = 0; k < grid->getMaximalOccupiedSecondTable(); k++) {
            IndexDimension Index = getIndexOfTable(k);
            Depth Tlocal(Index);
            if (Tlocal <= T) {
                AddPointDepth(Index, T);
            }
        }
    }
};

void MultiLevelAdaptiveSparseGrid::PrintActiveHangingDepth(int level, Depth T) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << "     " << "\n";
    cout << "active (o) and hanging (x) nodes on grid with ";
    T.Print();
    cout << endl;
    cout << "------------------------- " << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);
    unsigned long indexOfData;
    //unsigned int numberOfData = sparseGrid->numberOfData;
    //double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (occupied(indexOfData, I, T)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";

                if (isActiveNodeTable[indexOfData]) {
                    cout << " o ";
                } else { cout << " x "; }


            } else {
//	      cout << " - ";
                cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        std::cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        I = I.nextRight(1, level);

    }
};


bool MultiLevelAdaptiveSparseGrid::AddPointDepth(IndexDimension &Index, Depth &Tfine) {

/*        Depth Tlocal(Index);
        if(Tlocal>Tfine) return  false;*/

    unsigned long indArray = hash(Index, Tfine);

    unsigned long i = primeTable[indArray];

    if (i == 0) { // in prime table is something free;
        unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
        primeTable[indArray] = freeSpace + 1;
        isActiveNodeTable[freeSpace] = true;
        setIndexAndDepthInTable(Index, Tfine, freeSpace);
        return true;
    } else { // anaylse what is going on
        i = i - 1; // shift since data are stored with shift
        if (Index == getIndexOfTable(i) && Tfine == getDepthOfTable(i)) {
            // data is already stored but we still want all these nodes to be active!?
            isActiveNodeTable[i] = true;
            return false;
        } else { // search in second table;
            unsigned long iNext = secondTable[i];
            while (iNext > 1) {
                i = iNext - 2;
                if (Index == getIndexOfTable(i) && Tfine == getDepthOfTable(i)) {
                    isActiveNodeTable[i] = true;
                    return false;
                }
                iNext = secondTable[i];
            }
            myAssert(iNext == 1);
            unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
            secondTable[i] = freeSpace + 2;
            isActiveNodeTable[freeSpace] = true;
            setIndexAndDepthInTable(Index, Tfine, freeSpace);
            return true;
        }
    }
    if (i == 0) { // in prime table is something free;
        unsigned long freeSpace = getFreeSpaceNumberInSecondTable();
        primeTable[indArray] = freeSpace + 1;
        isActiveNodeTable[freeSpace] = true;

        setIndexAndDepthInTable(Index, Tfine, freeSpace);

        return true;
    }
    return false;
}

void
MultiLevelAdaptiveSparseGrid::setIndexAndDepthInTable(const IndexDimension Index, Depth Tfine, unsigned long iSetz) {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        indicesSecondTable[d + iSetz * DimensionSparseGrid] = Index.getIndex(d);
    }
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        depthTable[d + iSetz * DimensionSparseGrid] = Tfine.at(d);
    }
    myAssert(secondTable[iSetz] == 0);
    secondTable[iSetz] = 1;  // this means it is set, but no further link
    if (maximalOccupiedSecondTable <= iSetz) maximalOccupiedSecondTable = iSetz + 1;
}



unsigned long MultiLevelAdaptiveSparseGrid::getFreeSpaceNumberInSecondTable() {
    while(secondTable[minimalEmptySecondTable]!=0) {
        ++minimalEmptySecondTable;
        if(minimalEmptySecondTable>=secondTableLength) {
            cout << " resize of second table needed! not implemented! " << endl;
            cout << "minimalEmptySecondTable " << minimalEmptySecondTable << endl;

            assert(false);
        }
    }
    return minimalEmptySecondTable;
}


int  MultiLevelAdaptiveSparseGrid::getDOFS() {
    int dofs = 0;
    for (unsigned long k = 0; k < maximalOccupiedSecondTable; k++) {
        if (secondTable[k] != 0)
            if (isActiveNodeTable[k]) dofs++;
    }
    return dofs;
}