//
// Created by to35jepo on 6/24/24.
//

#ifndef RUN_ADDLOCALSTIFFNESSMATRICES_H
#define RUN_ADDLOCALSTIFFNESSMATRICES_H

#include "LocalStiffnessMatrices.h"


class ADDLocalStiffnessMatrices{
public:
    //ghost functions
    virtual inline void initialize(Depth &T_){};
    virtual inline void applyStencilOnCell(CellDimension& cell,VectorSparseG& input, VectorSparseG& output){};
    virtual inline double returnValue(const IndexDimension &Index, const MultiDimCompass &mc){return 0.0;};
    void applyLocalStiffnessMatricesOnIndices_onNode(VectorSparseG& input, VectorSparseG& output, Depth& D, std::vector<IndexDimension>& Indices){
        matricesA->applyLocalStiffnessMatricesOnIndices_onNode(input,output,D,Indices);
        matricesB->applyLocalStiffnessMatricesOnIndices_onNode(input,output,D,Indices);
    };

    ADDLocalStiffnessMatrices(LocalStiffnessMatrices* matricesA_, LocalStiffnessMatrices* matricesB_, AdaptiveSparseGrid& grid): matricesA(matricesA_),
        matricesB(matricesB_), input_node(grid), output_node(grid){
        number_processes=matricesB_->getNumberProcesses();
        resetActiveWorkers();
    };

    void applyLocalStiffnessMatricesDepth(VectorSparseG& input, VectorSparseG& output, Depth& depth){

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




    int nodeA = matricesA->getNodeForDepth(depth);
    int nodeB = matricesB->getNodeForDepth(depth);

    if(nodeA!=nodeB){cout << "ADDLocalStiffnessMatrices: ADD Local Stiffness Matrices only implemented for MemoryDistribution" << endl; exit(1);}

    if(rank == nodeA){


        matricesA->applyLocalStiffnessMatricesFixedDepth_onNode(input,output,depth);
        matricesB->applyLocalStiffnessMatricesFixedDepth_onNode(input,output,depth);

    }else{



        MPI_Status status;

        int length = int(input.getSparseGrid()->getMaximalOccupiedSecondTable());

        MPI_Isend(input.getDatatableVector(), length, MPI_DOUBLE, nodeA, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &send_requests[0]);
        MPI_Isend(output.getDatatableVector(), length, MPI_DOUBLE, nodeA, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &send_requests[1]);


        int depth_node_int2[DimensionSparseGrid];
        for(int d=0; d<DimensionSparseGrid; d++) depth_node_int2[d]=depth.at(d);
        int dim = int(DimensionSparseGrid);

        MPI_Isend(depth_node_int2, dim , MPI_INT, nodeA, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD,&send_requests[2]);



        MPI_Wait(&send_requests[0], &status);
        MPI_Wait(&send_requests[1], &status);
        MPI_Wait(&send_requests[2], &status);

        //braucht man das? man wird ja nie von nodeA einen Apply Auftrage bekommen
        /*MPI_Iprobe(nodeA, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            receiveApplySend(nodeA);
        }*/


        MPI_Irecv(input.getDatatableVector(), length, MPI_DOUBLE, nodeA, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD, &recv_requests[0]);
        MPI_Irecv(output.getDatatableVector(), length, MPI_DOUBLE, nodeA, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD, &recv_requests[1]);



        // Wait for all send and receive operations to complete
        MPI_Wait(&recv_requests[0], &status);
        MPI_Wait(&recv_requests[1], &status);





    }





#endif


    };

    void receiveApplySend(int n) {
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


    matricesA->applyLocalStiffnessMatricesFixedDepth_onNode(input_node, output_node, depth_node);
    matricesB->applyLocalStiffnessMatricesFixedDepth_onNode(input_node, output_node, depth_node);



    MPI_Send(input_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD);
    MPI_Send(output_node.getDatatableVector(), length, MPI_DOUBLE, n, TAG_LOCALSTIFFNESS_SENDBACK, MPI_COMM_WORLD);

#endif
    }
    TypeMatrixVectorMultiplication getTypeMatrixVectorMultiplication(){return  typeMatrixVectorMultiplication;};

    int getNodeForDepth(Depth& T){return matricesA->getNodeForDepth(T);};


    double applyLocalStiffnessMatricesFixedDepthIndex_onNode(VectorSparseG &input, VectorSparseG &output,
                                                                                                   Depth &depth, IndexDimension &Index) {

        double value = 0.0;
        value += matricesA->applyLocalStiffnessMatricesFixedDepthIndex_onNode(input,output,depth,Index);
        value += matricesB->applyLocalStiffnessMatricesFixedDepthIndex_onNode(input,output,depth,Index);


        return value;



    }


    void resetActiveWorkers(){
#ifdef MY_MPI_ON
        int num_tasks = 1;
            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
            active_worker.clear();
            for(int i=0; i< num_tasks-number_processes; i++) active_worker.push_back(i);
#endif
    };

    int getNumberProcesses(){return  number_processes;};




    void receiveApplySendOnActiveWorkers() {
#ifdef MY_MPI_ON
        MPI_Status status;
    int flag;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (auto it = active_worker.begin(); it != active_worker.end();){
        int n=*it;
        flag=0;
        MPI_Iprobe(n, TAG_LOCALSTIFFNESS_SEND, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            receiveApplySend(n);
        }
        flag=0;
        MPI_Iprobe(n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            int index;
            MPI_Recv(&index, 1, MPI_INT, status.MPI_SOURCE, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD, &status);

            it = active_worker.erase(it);


        }else{
            ++it;
        }
    }

#endif
    }
    std::vector<int> active_worker;
private:
    LocalStiffnessMatrices* matricesA;
    LocalStiffnessMatrices* matricesB;
    TypeMatrixVectorMultiplication typeMatrixVectorMultiplication = StoreLocalStiffnessMatrix;

    int number_processes=1;
    VectorSparseG input_node,output_node;
};

#endif //RUN_ADDLOCALSTIFFNESSMATRICES_H
