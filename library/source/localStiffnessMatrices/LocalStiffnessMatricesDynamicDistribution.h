//
// Created by to35jepo on 6/26/24.
//

#ifndef RUN_LOCALSTIFFNESSMATRICESDYNAMICDISTRIBUTION_H
#define RUN_LOCALSTIFFNESSMATRICESDYNAMICDISTRIBUTION_H


#include "LocalStiffnessMatrices.h"
#include "../mympi.h"
#include <sys/time.h>

#include <iostream>
#include <thread>
#include <chrono>
#include <random>

// Function to simulate random time passing
void randomSleep() {
    // Seed with a real random value, if available
    std::random_device rd;
    // Initialize random number generator
    std::mt19937 gen(rd());
    // Define the range for random sleep time (e.g., 1 to 5 seconds)
    std::uniform_int_distribution<> dis(0, 2);

    // Generate a random sleep duration
    int sleepTime = dis(gen);


    // Sleep for the randomly generated duration
    std::this_thread::sleep_for(std::chrono::seconds(sleepTime));


}

enum LS_TAGS{TAG_READY_WORKER=11, TAG_START_WORK=12, TAG_TERMINATE_WORK=13};

class LocalStiffnessMatricesDynamicDistribution: public LocalStiffnessMatrices{
public:
    LocalStiffnessMatricesDynamicDistribution(AdaptiveSparseGrid &sg, StencilTemplate& stencilClass, int number_processes_)
    :LocalStiffnessMatrices(sg, stencilClass, number_processes_){
        int rank;
        int size;
#ifdef MY_MPI_ON
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        MPI_Barrier(MPI_COMM_WORLD);



        int count=0;
        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
                  count++;
        }

        if(rank==0){

            master();


        }else if(rank>=size-number_processes_){

            worker();


        }


        MPI_Barrier(MPI_COMM_WORLD);



        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {

            Depth T = *it;
            int node;
            if(rank==0){
                node = distributedDepthsHashtable.getNodeForDepth(T);
            }

            MPI_Bcast(&node, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(rank>0){
                distributedDepthsHashtable.setNodeDepth(node,T);
            }

        }




        MPI_Barrier(MPI_COMM_WORLD);

#else
        cout << "LocalStiffness only with MPI!" <<endl;
        MPI_Finalize();
        exit(1);
#endif
      }

    void master(){
#ifdef MY_MPI_ON

        MPI_Status status;
        int d[DimensionSparseGrid];
        auto it = depthList.begin_all();
        int message=0;
        for (; it != depthList.end_all(); ++it) {
            Depth T = *it;

            MPI_Probe(MPI_ANY_SOURCE,TAG_READY_WORKER,MPI_COMM_WORLD,&status);
            MPI_Recv(&message,1, MPI_INT,status.MPI_SOURCE, TAG_READY_WORKER, MPI_COMM_WORLD, &status);
            for(int j=0; j<DimensionSparseGrid; j++)d[j]=T.at(j);

            MPI_Send(d,DimensionSparseGrid, MPI_INT,status.MPI_SOURCE,TAG_START_WORK,MPI_COMM_WORLD);
            distributedDepthsHashtable.setNodeDepth(status.MPI_SOURCE,T);
        }

        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        for(int i=size-number_processes; i < size; i++){

            MPI_Probe(i, TAG_READY_WORKER, MPI_COMM_WORLD, &status);
            MPI_Recv(&message,1, MPI_INT,status.MPI_SOURCE, TAG_READY_WORKER, MPI_COMM_WORLD, &status);

            MPI_Send(&message,1, MPI_INT,status.MPI_SOURCE,TAG_TERMINATE_WORK,MPI_COMM_WORLD);

        }

#endif
    }

    void worker(){

        #ifdef MY_MPI_ON

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        while(1){
            int message = 1; // Worker ready message
            int master_rank = 0; // Assuming master is rank 0


            // Send a message (with a tag) to the master indicating the worker is ready

            MPI_Send(&message, 1, MPI_INT, master_rank, TAG_READY_WORKER, MPI_COMM_WORLD);

            int d[DimensionSparseGrid];
            Depth D;
            MPI_Status status;
            int flag;
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD,&status);



            if(status.MPI_TAG == TAG_START_WORK){
                    MPI_Recv(d, DimensionSparseGrid, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    for(int j=0; j<DimensionSparseGrid; j++)D.set(d[j],j);
                    calculateEntries(D);
            }
            if(status.MPI_TAG == TAG_TERMINATE_WORK){

                MPI_Recv(&message, 1, MPI_INT, 0, TAG_TERMINATE_WORK, MPI_COMM_WORLD, &status);

                break;
            }

        }
#endif
    }


    void calculateEntries(Depth& depth);


    DistributedDepthsHashtable getDistributedDepthsHashtable(){return distributedDepthsHashtable;};


    void addStencil(StencilTemplate& newStencil){
        for(auto& item : pairedDepthsLocalStiffnessMatrices) {
            Depth T = item.first;
            newStencil.initialize(T);
            auto& matrices=item.second;
            int end=matrices.size();

#pragma omp parallel for
            for (int i=0; i<end;i++) {
                auto& matrix = matrices.at(i);
                CellDimension cell = *matrix.getCell();
                LocalStiffnessMatrixFixedDepthSymmetric matrix_add(cell, grid, newStencil);
                matrix+=matrix_add;

            }
        }
    }

    bool operator==(LocalStiffnessMatrices &localStiffnessMatrices){

        bool retbool = true;
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
                            if(!(matrix_add==matrix)){

                            retbool= false;
                            return false;
                            break;}
                        }
                    }
                }
            }
        }
        return retbool;
    }







};

void LocalStiffnessMatricesDynamicDistribution::calculateEntries(Depth &depth) {


    struct timeval begin_time, end_time;
    long seconds;
    long microseconds;
    gettimeofday(&begin_time, 0);


    stencilClass.initialize(depth);

    Depth Tcell = depth;
    ++Tcell;


    std::vector<LocalStiffnessMatrixFixedDepthSymmetric> vector;

    const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tcell);
    const auto &map = depthGrid._map;
    const auto end = depthGrid.getNumberOfEntries();


#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < end; i++) {
        CellDimension cellDimension;
        cellDimension = map.getIndexOfTable(i);
        LocalStiffnessMatrixFixedDepthSymmetric localStiffnessMatrixFixedDepthSymmetric(
                cellDimension, grid, stencilClass);


#pragma omp critical
        vector.push_back(localStiffnessMatrixFixedDepthSymmetric);
        numbercells++;
    }
    pairedDepthsLocalStiffnessMatrices.emplace_back(depth, vector);

    gettimeofday(&end_time, 0);
    seconds = end_time.tv_sec - begin_time.tv_sec;
    microseconds = end_time.tv_usec - begin_time.tv_usec;
    double t = seconds + microseconds * 1e-6;
    if(t>time) {
        time = t;
    }
}


#endif //RUN_LOCALSTIFFNESSMATRICESDYNAMICDISTRIBUTION_H
