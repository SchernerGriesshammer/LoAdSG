//
// Created by to35jepo on 6/26/24.
//

#ifndef RUN_LOCALSTIFFNESSMATRICESMEMORYDISTRIBUTION_H
#define RUN_LOCALSTIFFNESSMATRICESMEMORYDISTRIBUTION_H


#include <utility>

#include "LocalStiffnessMatrices.h"
#include <sys/time.h>

class LocalStiffnessMatricesMemoryDistribution: public LocalStiffnessMatrices{

public:
    LocalStiffnessMatricesMemoryDistribution(AdaptiveSparseGrid &sg, StencilTemplate& stencilClass, int number_processes_)
            :LocalStiffnessMatrices(sg, stencilClass, number_processes_){



        int num_tasks = 1;
        int rank = 0;

#ifdef MY_MPI_ON
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#else
        cout << "USE LOCALSTIFFNESSMATRICES ONLY WITH MPI" << endl;
        //exit(1);
#endif
        struct timeval begin_time, end_time;
        long seconds;
        long microseconds;
        gettimeofday(&begin_time, 0);

        Container container[number_processes];
        distribute_LocalStiffnessMatrices(number_processes, sg, container);

        distributedDepthsHashtable = DistributedDepthsHashtable(container, number_processes);


        gettimeofday(&end_time, 0);
        seconds = end_time.tv_sec - begin_time.tv_sec;
        microseconds = end_time.tv_usec - begin_time.tv_usec;
        double t2 = seconds + microseconds * 1e-6;
        if(rank==0)cout <<"time generate hashtable " << t2 << endl;

 /*       // Step 1: Collect keys
        std::vector<Depth> keys;
        keys.reserve(distributedDepthsHashtable.getMap()->size());
        for (const auto& entry : *distributedDepthsHashtable.getMap()) {
            keys.push_back(entry.first);
        }
        //printEstimatedStorage();
        int desiredValue = rank;


#pragma omp parallel for schedule(dynamic)
        for(std::size_t i = 0; i < keys.size(); ++i){
            Depth T = keys[i];
            int value = distributedDepthsHashtable.getNodeForDepth(T);
            if(value == desiredValue)
            {


                stencilClass.initialize(T);

                Depth Tcell = T;
                ++Tcell;


                std::vector<LocalStiffnessMatrixFixedDepthSymmetric> vector;

                const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tcell);
                const auto &map = depthGrid._map;
                const auto end = depthGrid.getNumberOfEntries();



                for (size_t i = 0; i < end; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    LocalStiffnessMatrixFixedDepthSymmetric localStiffnessMatrixFixedDepthSymmetric(cellDimension, sg, stencilClass);
                    numbercells++;



                    vector.push_back(localStiffnessMatrixFixedDepthSymmetric);
                }
                #pragma omp critical
                pairedDepthsLocalStiffnessMatrices.push_back(std::make_pair(T, vector));
            }
        }*/

 int j=0;
        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
            Depth T = *it;

            int node = distributedDepthsHashtable.getNodeForDepth(T);

            if (rank == node) {

                struct timeval begin_time, end_time, end_time2;
                long seconds;
                long microseconds;
                gettimeofday(&begin_time, 0);

                stencilClass.initialize(T);

                Depth Tcell = T;
                ++Tcell;




                const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tcell);
                const auto &map = depthGrid._map;
                const auto end = depthGrid.getNumberOfEntries();

                //LocalStiffnessMatrixFixedDepthSymmetric ls_dummy(grid,stencilClass);

                //std::vector<LocalStiffnessMatrixFixedDepthSymmetric> vector(end,ls_dummy);
                std::vector<LocalStiffnessMatrixFixedDepthSymmetric> vector;
                gettimeofday(&begin_time, 0);
#pragma omp parallel for schedule(dynamic)
                for (size_t i = 0; i < end; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    LocalStiffnessMatrixFixedDepthSymmetric localStiffnessMatrixFixedDepthSymmetric(cellDimension, sg, stencilClass);
                    numbercells++;

//vector[i]=localStiffnessMatrixFixedDepthSymmetric;
#pragma omp critical
                    vector.push_back(localStiffnessMatrixFixedDepthSymmetric);
                }
                gettimeofday(&end_time2, 0);
                pairedDepthsLocalStiffnessMatrices.push_back(std::make_pair(T, vector));

                seconds = end_time2.tv_sec - begin_time.tv_sec;
                microseconds = end_time2.tv_usec - begin_time.tv_usec;
                double t = seconds + microseconds * 1e-6;
                if(t>time ) {
                    time = t;
                    numbercells_max=end;
                    Tmax=T;
                }
                if(t<time_min ){
                    time_min=t;
                    numbercells_min=end;

                }

            }
        j++;
        }

        if(rank==0)cout << " number of depths " << j << endl;

    };
};


#endif //RUN_LOCALSTIFFNESSMATRICESMEMORYDISTRIBUTION_H
