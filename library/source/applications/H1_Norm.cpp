//
// Created by to35jepo on 7/30/24.
//
#include "H1_Norm.h"
#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"


double H1::getValue(VectorSparseG &u){

    Poisson stencilPoisson(sparseGrid);
    int world_rank;
    int world_size;
#ifdef MY_MPI_ON

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



    VectorSparseG v(sparseGrid);
    VectorSparseG Av(sparseGrid);
    calcPrewByNodal(v,u);

    MatrixVectorHomogen matrixVectorHomogen(sparseGrid, mgrid, 1, 0);
    matrixVectorHomogen.multiplication(v,Av,stencilPoisson);

    double sum=0.0;
    for (size_t i = 0; i < sparseGrid.getMaximalOccupiedSecondTable(); i++) {
        if (sparseGrid.getActiveTable()[i]) {
            sum+=v.getValue(i)*Av.getValue(i);

        }


    }
    return sqrt(sum);



}


double L2::getValue(VectorSparseG &u){

    HelmHoltz stencilPoisson(sparseGrid);
    int world_rank;
    int world_size;
#ifdef MY_MPI_ON

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



    VectorSparseG v(sparseGrid);
    VectorSparseG Av(sparseGrid);
    calcPrewByNodal(v,u);

    MatrixVectorHomogen matrixVectorHomogen(sparseGrid, mgrid, 1, 0);
    matrixVectorHomogen.multiplication(v,Av,stencilPoisson);

    double sum=0.0;
    for (size_t i = 0; i < sparseGrid.getMaximalOccupiedSecondTable(); i++) {
        if (sparseGrid.getActiveTable()[i]) {
            sum+=v.getValue(i)*Av.getValue(i);

        }


    }
    return sqrt(sum);



}