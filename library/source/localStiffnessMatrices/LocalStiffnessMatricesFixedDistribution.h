//
// Created by to35jepo on 7/12/24.
//

#ifndef RUN_LOCALSTIFFNESSMATRICESFIXEDDISTRIBUTION_H
#define RUN_LOCALSTIFFNESSMATRICESFIXEDDISTRIBUTION_H


#include <utility>

#include "LocalStiffnessMatrices.h"

class LocalStiffnessMatricesFixedDistribution : public LocalStiffnessMatrices{


public:
    LocalStiffnessMatricesFixedDistribution(AdaptiveSparseGrid &sg, StencilTemplate& stencilClass, int number_processes_, DistributedDepthsHashtable distributedDepthsHashtable_ )
    :LocalStiffnessMatrices(sg, stencilClass, number_processes_){
        distributedDepthsHashtable = distributedDepthsHashtable_;
        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
            Depth T = *it;

            int node = distributedDepthsHashtable.getNodeForDepth(T);
            int node_vergleich = distributedDepthsHashtable_.getNodeForDepth(T);
            if (node != node_vergleich) {
                cout << " err " << endl;
                exit(1);
            }
        }

        int num_tasks = 1;
        int rank = 0;

#ifdef MY_MPI_ON
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#else
        cout << "USE LOCALSTIFFNESSMATRICES ONLY WITH MPI" << endl;
        //exit(1);
#endif

        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
            Depth T = *it;

            int node = distributedDepthsHashtable.getNodeForDepth(T);
            int node_vergleich = distributedDepthsHashtable_.getNodeForDepth(T);
            if (node != node_vergleich) {
                cout << " err " << endl;
                exit(1);
            }

            if (rank == node) {

                stencilClass.initialize(T);

                Depth Tcell = T;
                ++Tcell;




                const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tcell);
                const auto &map = depthGrid._map;
                const auto end = depthGrid.getNumberOfEntries();



                std::vector<LocalStiffnessMatrixFixedDepthSymmetric> vector;

#pragma omp parallel for schedule(dynamic)
                for (size_t i = 0; i < end; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    LocalStiffnessMatrixFixedDepthSymmetric localStiffnessMatrixFixedDepthSymmetric(cellDimension, sg, stencilClass);
                    numbercells++;


#pragma omp critical
                    vector.push_back(localStiffnessMatrixFixedDepthSymmetric);
                }
                pairedDepthsLocalStiffnessMatrices.push_back(std::make_pair(T, vector));


            }

        }

    }

    void operator+=(LocalStiffnessMatrices& localStiffnessMatrices){
        for(auto& item : pairedDepthsLocalStiffnessMatrices){
            Depth T = item.first;

            auto it = std::find_if(
                    localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->begin(),
                    localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->end(),
                    [&T](const std::pair<Depth,std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>>& pair) {
                        return pair.first == T;
                    }
            );

            if (it !=  localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->end()) {
                // Found the matching entry
                std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric> &matrices_add = it->second;

                for(auto& matrix : item.second){
                    for (auto &matrix_add: matrices_add){
                        CellDimension cellCompare=*matrix.getCell();
                     if(cellCompare==*matrix_add.getCell()){
                         matrix+=matrix_add;
                     }
                    }
                }
            }
        }
    }

    LocalStiffnessMatrices& operator+(LocalStiffnessMatrices& localStiffnessMatrices){
        for(auto& item : pairedDepthsLocalStiffnessMatrices){
            Depth T = item.first;

            auto it = std::find_if(
                    localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->begin(),
                    localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->end(),
                    [&T](const std::pair<Depth,std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>>& pair) {
                        return pair.first == T;
                    }
            );

            if (it !=  localStiffnessMatrices.getPairedDepthsLocalStiffnessMatrices()->end()) {
                // Found the matching entry
                std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric> &matrices_add = it->second;

                for(auto& matrix : item.second) {
                    for (auto &matrix_add: matrices_add){
                        CellDimension cellCompare=*matrix.getCell();
                        if(cellCompare==*matrix_add.getCell()){
                            matrix+=matrix_add;
                        }
                    }
                }
            }
        }
        return *this;
    }
};


#endif //RUN_LOCALSTIFFNESSMATRICESFIXEDDISTRIBUTION_H
