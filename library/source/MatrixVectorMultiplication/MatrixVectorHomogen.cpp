//
// Created by to35jepo on 12/7/22.
//

#include "MatrixVectorHomogen.h"




void MatrixVectorHomogen::calcNodal(VectorSparseG &prew, Depth &Tiefe) {


    MultiDimFiveCompass mc;
    SingleDepthHashGrid &depthGrid = prew.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(Tiefe);
    const auto &mapping = depthGrid._mapPosToGridPos;


    for (size_t i = 0; i < mapping.size(); i++) {
        if (prew.getSparseGrid()->getActiveTable()[mapping[i]]) {

            IndexDimension Index = depthGrid._map.getIndexOfTable(i);
            double coeff = prew.getValue(mapping[i]);

            IndexDimension J;
            mc.goToStart();
            for (; mc.hasNext(); ++mc) {
                double basis_coeff = 1.0;
                J = Index.nextFiveP(&mc, Tiefe, &basis_coeff);
                if (mc.goon()) {
                    unsigned long k;
                    if (prew.getSparseGrid()->occupied(k, J)) {
                        double val = u.getValue(k) + coeff * basis_coeff;
                        u.setValue(k, val);
                    }
                }
            }


        }
    }
}

void MatrixVectorHomogen::calcNodal(VectorSparseG &nodal_u, VectorSparseG &prew, Depth &Tiefe) {


    MultiDimFiveCompass mc;
    SingleDepthHashGrid &depthGrid = prew.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(Tiefe);
    const auto &mapping = depthGrid._mapPosToGridPos;


    for (size_t i = 0; i < mapping.size(); i++) {
        if (prew.getSparseGrid()->getActiveTable()[mapping[i]]) {

            IndexDimension Index = depthGrid._map.getIndexOfTable(i);
            double coeff = prew.getValue(mapping[i]);

            IndexDimension J;
            mc.goToStart();
            for (; mc.hasNext(); ++mc) {
                double basis_coeff = 1.0;
                J = Index.nextFiveP(&mc, Tiefe, &basis_coeff);
                if (mc.goon()) {
                    unsigned long k;
                    if (prew.getSparseGrid()->occupied(k, J)) {
                        double val = nodal_u.getValue(k) + coeff * basis_coeff;
                        nodal_u.setValue(k, val);
                    }
                }
            }


        }
    }
}

void MatrixVectorHomogen::master(int num_workers, int total_tasks) {
#ifdef MY_MPI_ON
    int task_index = 0;
    int worker_rank;
    MPI_Status status;

    int old_index;
    // Distribute initial tasks to workers
    for (worker_rank = 1; worker_rank <= num_workers; worker_rank++) {

        MPI_Send(&task_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD);

        task_index++;
    }



    // Receive results from workers and assign new tasks


    while (task_index < total_tasks) {
        MPI_Recv(&old_index, 1, MPI_INT, MPI_ANY_SOURCE, TAG_TASK, MPI_COMM_WORLD, &status);
        worker_rank = status.MPI_SOURCE;
        MPI_Send(&task_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD);
        //MPI_Recv(&confirmation, 1, MPI_INT, worker_rank, TAG_CONFIRMATION, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        task_index ++;
    }


    // Send termination signal to workers
    for (worker_rank = 1; worker_rank <= num_workers; worker_rank++) {
        //receive message, that worker has completed all of its tasks
        MPI_Recv(&old_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //send termination signal to worker
        MPI_Send(&task_index, 1, MPI_INT, worker_rank, TAG_TERMINATE, MPI_COMM_WORLD);


    }
#endif

}




void MatrixVectorHomogen::master_start(int num_workers){
#ifdef MY_MPI_ON




    MPI_Status status;
    MPI_num_workers=min(num_workers,MPI_total_tasks);
    MPI_task_index=0;
    //Distribute initial tasks to workers



    for (int worker_rank = mpi_per_worker; (worker_rank < num_workers) && (MPI_task_index < MPI_total_tasks); worker_rank+=mpi_per_worker) {

        for(int j=0; j < mpi_per_worker; j++){

           MPI_Send(&MPI_task_index, 1, MPI_INT, worker_rank+j, TAG_TASK, MPI_COMM_WORLD);
        }

        MPI_task_index++;
    }



#endif

}

void MatrixVectorHomogen::master_start_onlyCases(int num_workers) {
#ifdef MY_MPI_ON




    MPI_Status status;
    MPI_num_workers=min(num_workers,MPI_total_tasks);
    MPI_task_index=0;

    int num_tasks;




    for (int worker_rank = 1; (worker_rank < MPI_num_workers) && (MPI_task_index < MPI_total_tasks); worker_rank++) {

        MPI_Send(&MPI_task_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD);


        MPI_task_index++;
    }


#endif

}



void MatrixVectorHomogen::master_distribute() {

#ifdef MY_MPI_ON


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Status status;
    int flag=0;
    int old_index;





    for (int worker = mpi_per_worker; worker < MPI_num_workers && MPI_task_index < MPI_total_tasks; worker+=mpi_per_worker) {

         bool received_message=false;
        for(int j=0; j< mpi_per_worker; j++){

            MPI_Iprobe(worker+j, TAG_TASK, MPI_COMM_WORLD, &flag, &status);

            if (flag){
                //received flag by worker+j, now has to wait for all corresponding workers, to send also their message

                received_message=true;

                break;
            }
        }

        if(received_message){
               for(int k=0; k< mpi_per_worker; k++){

                    MPI_Recv(&old_index, 1, MPI_INT, worker+k, TAG_TASK, MPI_COMM_WORLD, &status);

                    MPI_Send(&MPI_task_index, 1, MPI_INT, worker+k, TAG_TASK, MPI_COMM_WORLD);
                }
               MPI_task_index++;
        }

    }







#endif

}


void MatrixVectorHomogen::master_end_distribute() {

#ifdef MY_MPI_ON




        int old_index;


        // Send termination signal to workers
        for (int worker_rank = mpi_per_worker; worker_rank < MPI_num_workers; worker_rank++) {



            MPI_Recv(&old_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            //send termination signal to worker#

            MPI_Send(&MPI_task_index, 1, MPI_INT, worker_rank, TAG_TERMINATE, MPI_COMM_WORLD);


        }



#endif

}




void MatrixVectorHomogen::master_distribute_onlyCases() {
#ifdef MY_MPI_ON



    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Status status;
    int flag=0;
    int old_index;





    for (int worker =1; worker < MPI_num_workers && MPI_task_index < MPI_total_tasks; worker++) {

            MPI_Iprobe(worker, TAG_TASK, MPI_COMM_WORLD, &flag, &status);

            if (flag){

                    MPI_Recv(&old_index, 1, MPI_INT, worker, TAG_TASK, MPI_COMM_WORLD, &status);

                    MPI_Send(&MPI_task_index, 1, MPI_INT, worker, TAG_TASK, MPI_COMM_WORLD);
                    MPI_task_index++;
            }


    }







#endif

};




void MatrixVectorHomogen::distributeAndApplyLocalStiffnessMatrix(LocalStiffnessMatrices &localStiffnessMatrices) {

#ifdef MY_MPI_ON
    int num_tasks = 1;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);



    int flag=0;
    MPI_Status status;
    std::vector<int> vec;
    for(int n=0; n< num_tasks; n++){
        if(n!=rank){
            vec.push_back(n);
            int sendI;
            MPI_Send(&sendI,1, MPI_INT, n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD);
        }
    }



    while(vec.size()>0){
        for (auto it = vec.begin(); it != vec.end();++it ){
            int n=*it;
            MPI_Iprobe(n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                localStiffnessMatrices.receiveApplySend(n);
            }
            MPI_Iprobe(n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &flag, &status);
            if (flag) vec.erase(it);

        }
    }
#endif
}
