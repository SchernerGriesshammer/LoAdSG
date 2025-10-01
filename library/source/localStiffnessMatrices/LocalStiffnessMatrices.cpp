//
// Created by to35jepo on 6/6/24.
//

#include "LocalStiffnessMatrices.h"












double bytesToGigabytes(unsigned long long int bytes) {
    const unsigned long long bytesPerGigabyte = 1073741824; // 1024^3
    return static_cast<double>(bytes) / bytesPerGigabyte;
}

void distribute_LocalStiffnessMatrices(int numberofprocesses, AdaptiveSparseGrid &sg, Container *containers_sorted) {
    DepthList depthList(sg);

    CellData cellData(sg);


    // Depth + number of cells per Depth
    std::vector<std::pair<Depth,int>> boxes;
    for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
        Depth T = *it;
        Depth Tcell= T;
        ++Tcell;
        const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tcell);
        int size = depthGrid.getNumberOfEntries();

        boxes.push_back(std::make_pair(T,size));
    }

    // Sort the boxes in descending order based on size
    std::sort(boxes.begin(), boxes.end(), [](const std::pair<Depth, int>& a, const std::pair<Depth, int>& b) {
        return a.second > b.second;
    });


    // Create a priority queue to hold the containers by their total size (smallest first)
    std::priority_queue<Container, std::vector<Container>, std::greater<Container>> containers;

    // Initialize the containers
    int num_tasks = 1;
#ifdef MY_MPI_ON
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif

    for (int i = 0; i < numberofprocesses; ++i) {
        containers.push({i, 0});
    }

    // Distribute the boxes
    for (const auto& box : boxes) {
        // Get the container with the smallest total size
        Container smallest = containers.top();
        containers.pop();

        if(bytesToGigabytes(smallest.totalSize*TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value*25) >= 200){
            cout << "Container (smallest) on node " << smallest.id << " is full! Use more nodes! " << endl;
            exit(1);
        }

        // Add the box to this container
        smallest.totalSize += box.second;

        smallest.boxes.push_back(box);


        // Put the container back into the priority queue
        containers.push(smallest);
    }




    while (!containers.empty()) {
        Container c = containers.top();
        containers.pop();
        containers_sorted[c.id]=c;
    }

}

int DistributedDepthsHashtable::getNodeForDepth(const Depth &D) {
    auto it = _map.find(D);
    if (it != _map.end()) {
        return it->second;
    }
    // Return a default value or throw an exception if the depth is not found
    return 0; // Default value
}

DistributedDepthsHashtable::DistributedDepthsHashtable(Container *containers_sorted, int n) : _map{}{

    int num_tasks = 1;
#ifdef MY_MPI_ON
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif
    for(int i=1;i<=n;i++){
        Container c = containers_sorted[i-1];

        for (const auto& box : c.boxes) {
            Depth T= box.first;
            int node = num_tasks-i;
            _map.insert({T, node});
        }
    }
}


void LocalStiffnessMatrices::applyLocalStiffnessMatricesDepth(VectorSparseG &input, VectorSparseG &output,
                                                              Depth &depth){

    if(distribute){

#ifdef MY_MPI_ON
        int num_tasks = 1;
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

        int flag=0;


        int node = distributedDepthsHashtable.getNodeForDepth(depth);
        MPI_Status status;

        int length = int(input.getSparseGrid()->getMaximalOccupiedSecondTable());

        MPI_Send(input.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD);
        MPI_Send(output.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD);


        int depth_node_int[DimensionSparseGrid];
        for(int d=0; d<DimensionSparseGrid; d++) depth_node_int[d]=depth.at(d);
        int dim = int(DimensionSparseGrid);

        MPI_Send(depth_node_int, dim , MPI_INT, node, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD);

        MPI_Recv(input.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD,&status);
        MPI_Recv(output.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD,&status);




#else
        applyLocalStiffnessMatricesFixedDepth_onNode(input,output,depth);
#endif

    }else{
        applyLocalStiffnessMatricesFixedDepth_onNode(input,output,depth);
    }



}


void LocalStiffnessMatrices::applyLocalStiffnessMatricesFixedDepth(VectorSparseG &input, VectorSparseG &output,
                                                                   Depth &depth){



#ifdef MY_MPI_ON
    int num_tasks = 1;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Status status;
    int flag=0;
    MPI_Request send_requests[3], recv_requests[3];

    Depth depth_node;
    int depth_node_int[DimensionSparseGrid];




    int node = distributedDepthsHashtable.getNodeForDepth(depth);

    if(rank == node){

        applyLocalStiffnessMatricesFixedDepth_onNode(input,output,depth);
    }else{

        MPI_Status status;

        int length = int(input.getSparseGrid()->getMaximalOccupiedSecondTable());

        MPI_Isend(input.getDatatableVector(), length, MPI_DOUBLE, 1, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &send_requests[0]);
        MPI_Isend(output.getDatatableVector(), length, MPI_DOUBLE, 1, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &send_requests[1]);


        int depth_node_int2[DimensionSparseGrid];
        for(int d=0; d<DimensionSparseGrid; d++) depth_node_int2[d]=depth.at(d);
        int dim = int(DimensionSparseGrid);

        MPI_Isend(depth_node_int2, dim , MPI_INT, node, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD,&send_requests[2]);



        MPI_Wait(&send_requests[0], &status);
        MPI_Wait(&send_requests[1], &status);
        MPI_Wait(&send_requests[2], &status);

        MPI_Iprobe(node, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            receiveApplySend(node);
        }


        MPI_Irecv(input.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD, &recv_requests[0]);
        MPI_Irecv(output.getDatatableVector(), length, MPI_DOUBLE, node, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD, &recv_requests[1]);



        // Wait for all send and receive operations to complete
        MPI_Wait(&recv_requests[0], &status);
        MPI_Wait(&recv_requests[1], &status);





    }





#endif





}


void LocalStiffnessMatrices::applyLocalStiffnessMatricesFixedDepth_onNode(VectorSparseG &input, VectorSparseG &output,
                                                                          Depth &depth) {
    auto it = std::find_if(
            pairedDepthsLocalStiffnessMatrices.begin(),
            pairedDepthsLocalStiffnessMatrices.end(),
            [&depth](const std::pair<Depth,std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>>& pair) {
                return pair.first == depth;
            }
    );

    if (it != pairedDepthsLocalStiffnessMatrices.end()) {
        // Found the matching entry
        std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>& matrices = it->second;

        //all local stiffness matrices
        for (auto& matrix : matrices) {
            matrix.applyLocalMatrix(input,output);
        }
    } else {
        std::cout << "Error LocalStiffnessMatrix : No matching entry found." << std::endl;
        exit(1);
    }
}
double LocalStiffnessMatrices::applyLocalStiffnessMatricesFixedDepthIndex_onNode(VectorSparseG &input, VectorSparseG &output,
                                                                                 Depth &depth, IndexDimension &Index) {


    auto it = std::find_if(
            pairedDepthsLocalStiffnessMatrices.begin(),
            pairedDepthsLocalStiffnessMatrices.end(),
            [&depth](const std::pair<Depth,std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>>& pair) {
                return pair.first == depth;
            }
    );


    if (it != pairedDepthsLocalStiffnessMatrices.end()){
        // Found the matching entry
        std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>& matrices = it->second;

        int i=0;
        int maxNumberCellsTouchingIndex = PowerOfTwo<DimensionSparseGrid>::value;

        // Range-based for loop
        for (auto& matrix : matrices) {
            if(matrix.indexInMatrix(Index)){
                matrix.applyLocalMatrixIndex(input,output,Index);
                i++;
            }
            if(i>maxNumberCellsTouchingIndex) break;
        }
    } else {
        std::cout << "Error LocalStiffnessMatrix : No matching entry found." << std::endl;
        exit(1);
    }

    double return_value=output.getValue(Index);
    output=0.0;
    return return_value;



}
void LocalStiffnessMatrices::applyLocalStiffnessMatricesOnIndices_onNode(VectorSparseG &input,
                                                                         VectorSparseG &output,
                                                                         Depth &depth,
                                                                         vector<IndexDimension> &Indices) {

    auto it = std::find_if(
            pairedDepthsLocalStiffnessMatrices.begin(),
            pairedDepthsLocalStiffnessMatrices.end(),
            [&depth](const std::pair<Depth,std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>>& pair) {
                return pair.first == depth;
            }
    );


    if (it != pairedDepthsLocalStiffnessMatrices.end()){
        // Found the matching entry
        std::vector<LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric>& matrices = it->second;

        int i=0;
        int maxNumberCellsTouchingIndex = PowerOfTwo<DimensionSparseGrid>::value;
        maxNumberCellsTouchingIndex*=int(Indices.size());

        bool apply=false;
        // Range-based for loop
        for (auto& matrix : matrices) {
            for(auto& Index : Indices) {
                if (matrix.indexInMatrix(Index)){
                    apply=true;
                    break;
                }
            }
            if(apply){
                i++;
                matrix.applyLocalMatrix(input, output);
            }
            apply=false;
            if(i>maxNumberCellsTouchingIndex)break;
        }
    } else {
        std::cout << "Error LocalStiffnessMatrix : No matching entry found." << std::endl;
        exit(1);
    }

}




void LocalStiffnessMatrices::receiveApplySend(int n) {
#ifdef MY_MPI_ON
    int depth_node_int[DimensionSparseGrid];
    Depth depth_node;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);



    // receive input, output, depth
    int length = int(input_node.getSparseGrid()->getMaximalOccupiedSecondTable());
    MPI_Recv(input_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD,&status);
    MPI_Recv(output_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD,&status);

    int dim = int(DimensionSparseGrid);
    MPI_Recv(depth_node_int, dim, MPI_INT, n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD,&status);

    for (int d = 0; d < DimensionSparseGrid; d++) depth_node.set(depth_node_int[d], d);
    applyLocalStiffnessMatricesFixedDepth_onNode(input_node, output_node, depth_node);



    MPI_Send(input_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD);
    MPI_Send(output_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD);

#endif
}


void LocalStiffnessMatrices::receiveApplySendOnActiveWorkers() {
#ifdef MY_MPI_ON
    MPI_Status status;
    int flag;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(pairedDepthsLocalStiffnessMatrices.size()>0) {
        for (auto it = active_worker.begin(); it != active_worker.end();) {
            int n = *it;
            flag = 0;
            MPI_Iprobe(n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                receiveApplySend(n);
            }
            flag = 0;
            MPI_Iprobe(n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int index;
                MPI_Recv(&index, 1, MPI_INT, status.MPI_SOURCE, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &status);

                it = active_worker.erase(it);


            } else {
                ++it;
            }
        }
    }
    else {
        for (auto it = active_worker.begin(); it != active_worker.end();) {
            int n = *it;
            flag = 0;
            MPI_Iprobe(n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
                int index;
                MPI_Recv(&index, 1, MPI_INT, status.MPI_SOURCE, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &status);

                it = active_worker.erase(it);


            } else {
                ++it;
            }
        }
    }
#endif
}



void LocalStiffnessMatrices::resetActiveWorkers() {
    active_worker.clear();
#ifdef MY_MPI_ON
    int num_tasks = 1;
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
        for(int i=0; i< num_tasks-number_processes; i++) active_worker.push_back(i);
#endif
}



void LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric::setValues(CellDimension &cell) {

    int j=0;



    for (CellIndexIterator outerIter2(&cell); outerIter2.goon(); ++outerIter2){

        IndexDimension p = outerIter2.getIndex();
        unsigned long kp;
        bool occP =grid.occupied(kp, p);
        CellIndexDirection dirP = outerIter2.getCellIndexDirection();


        for (CellIndexIterator innerIter(outerIter2); innerIter.goon(); ++innerIter) {
            IndexDimension q = innerIter.getIndex();

            if (q.isNotAtBoundary()) {

                CellIndexDirection dirQ = innerIter.getCellIndexDirection();

                unsigned long kq;

                bool occQ = grid.occupied(kq, q);
                if (occQ && occP) {
                    occupied[j] = true;

                    storage_p[j]=kp;
                    storage_q[j]=kq;
                    double val = stencilTemplate.integration(cell, dirP, dirQ);
                    entries[j] = val;

                } else {
                    occupied[j] = false;
                }

            } else {
                occupied[j] = false;
            }
            j++;
        }
    }






}


void LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric::setValues(CellDimension &cell, StencilTemplate& stencil) {

    int j=0;



    for (CellIndexIterator outerIter2(&cell); outerIter2.goon(); ++outerIter2){

        IndexDimension p = outerIter2.getIndex();
        unsigned long kp;
        bool occP =grid.occupied(kp, p);
        CellIndexDirection dirP = outerIter2.getCellIndexDirection();


        for (CellIndexIterator innerIter(outerIter2); innerIter.goon(); ++innerIter) {
            IndexDimension q = innerIter.getIndex();

            if (q.isNotAtBoundary()) {

                CellIndexDirection dirQ = innerIter.getCellIndexDirection();

                unsigned long kq;

                bool occQ = grid.occupied(kq, q);
                if (occQ && occP) {
                    occupied[j] = true;

                    storage_p[j]=kp;
                    storage_q[j]=kq;
                    double val = stencil.integration(cell, dirP, dirQ);
                    entries[j] = val;

                } else {
                    occupied[j] = false;
                }

            } else {
                occupied[j] = false;
            }
            j++;
        }
    }






}


void LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric::applyLocalMatrix(VectorSparseG &input,
                                                                                       VectorSparseG &output) {

    for (int i = 0; i < array_size; i++) {
        if (occupied[i]) {
            double val =entries[i];


            unsigned long kp = storage_p[i];
            unsigned long kq = storage_q[i];

            if (kp==kq) {
                val = val * input.getValue(kp);

                output.addToValue(kq, val);
            } else {
                double val_p = val * input.getValue(kp);

                output.addToValue(kq, val_p);

                double val_q = val * input.getValue(kq);

                output.addToValue(kp, val_q);
            }

        }
    }
}



void LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric::applyLocalMatrixIndex(VectorSparseG &input,
                                                                                            VectorSparseG &output,
                                                                                            IndexDimension &Index) {

    unsigned long storageIndex;
    if(input.getSparseGrid()->occupied(storageIndex,Index)){
        for (int i = 0; i < array_size; i++) {
            if (occupied[i]) {
                double val = entries[i];

                unsigned long kp = storage_p[i];
                unsigned long kq = storage_q[i];


                if (kp==kq) {
                    val = val * input.getValue(kp);
                    output.addToValue(kq, val);
                } else {
                    double val_p = val * input.getValue(kp);
                    output.addToValue(kq, val_p);

                    double val_q = val * input.getValue(kq);
                    output.addToValue(kp, val_q);
                }


            }
        }
    }
}

bool LocalStiffnessMatrices::LocalStiffnessMatrixFixedDepthSymmetric::indexInMatrix(IndexDimension &Index) {
    unsigned long storageIndex;
    if(grid.occupied(storageIndex,Index)){
        for (int i = 0; i < array_size; i++){
            if (occupied[i]) {
                unsigned long kp = storage_p[i];
                unsigned long kq = storage_q[i];
                if(kp==storageIndex || kq==storageIndex)return true;
            }
        }
    }
    return false;
}




void LocalStiffnessMatrices::printEstimatedStorage() {
    int a = int(numbercells);
    int b = int(TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value);


    cout <<a<< " x " <<b<<" x 25"<< " = " <<a*b*25 <<" " << bytesToGigabytes(a*b*25)<<" GB " <<  endl;

    for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
        Depth T = *it;

        int node = distributedDepthsHashtable.getNodeForDepth(T);
        cout << " Depth " << T.at(0) << T.at(1) << " stored in node " << node << endl;
    }

}
