//
// Created by scherner on 31.05.21.
//

#include "multilevelvector.h"



MultiLevelVector::MultiLevelVector(MultiLevelAdaptiveSparseGrid &grid) {
    sparseGrid = &grid;
    //number = sparseGrid->getNewDataNumber();
    unsigned long length = sparseGrid->getLengthSecondTable();

    dataTableVector = new double[length];
};


MultiLevelVector::~MultiLevelVector() {
    delete[] dataTableVector;
}


void MultiLevelVector::Broadcast(int rank) {
    int length = int(sparseGrid->getLengthSecondTable());
    MPI_Bcast(dataTableVector, length, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}





void MultiLevelVector::ReduceSum(int rank) {
    int length = int(sparseGrid->getLengthSecondTable());

    auto *data = new double[length];
    MPI_Reduce(dataTableVector, data, length, MPI_DOUBLE, MPI_SUM, rank, MPI_COMM_WORLD);
    dataTableVector = data;

}


void MultiLevelVector::PrintDoubleTwoD(int level, Depth T) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << " integer values on grid " << endl;

    cout << "------------------------- " << endl;
    cout << endl;
    IndexDimension I;

    I = IndexOneD(startUnitInterval);
    unsigned long int indexOfData = 0;
    ////unsigned int numberOfData = sparseGrid->numberOfData;
    ////double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (sparseGrid->occupied(indexOfData, I, T)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
                //     cout << " " << indexOfData << " ";
                double c = dataTableVector[indexOfData];
                std::cout << c << "  ";
            } else {
//	      cout << " - ";
                cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        I = I.nextRight(1, level);

    }
    cout << endl;
    cout << "------------------------- " << endl;
    cout << endl;

}

void MultiLevelVector::PrintDoubleTwoD(int level) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << " integer values on grid " << endl;

    cout << "------------------------- " << endl;
    cout << endl;
    IndexDimension I;

    I = IndexOneD(startUnitInterval);
    unsigned long int indexOfData = 0;
    ////unsigned int numberOfData = sparseGrid->numberOfData;
    ////double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Depth T(I);
            if (sparseGrid->occupied(indexOfData, I, T)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
                //     cout << " " << indexOfData << " ";
                double c = dataTableVector[indexOfData];
                std::cout << c << "  ";
            } else {
//	      cout << " - ";
                cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        I = I.nextRight(1, level);

    }
    cout << endl;
    cout << "------------------------- " << endl;
    cout << endl;

}






