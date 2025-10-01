//
// Created by to35jepo on 1/15/24.
//
#include "mycomm.h"
#include "abbrevi.h"

MPI_Comm case_comm;
MPI_Comm case_comm_base;
void create_my_communicator() {
//#ifdef MY_MPI_ON

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int returnvalue = 1;
    int base = 2;
    for (int j = 0; j < DimensionSparseGrid; j++) {
        returnvalue = returnvalue * base;
    }
    int num_cases=returnvalue;
    int color = world_rank % num_cases; // Determine color based on row
    if(world_size%num_cases!=0){
        std::cout << "ERROR: number of MPI_PROCESSES has to be a multiple of 2^DIMENSION" << std::endl;
        MPI_Finalize();
        exit(0);
    }

// Split the communicator based on the color and use the
// original rank for ordering
    /*MPI_Comm case_comm;*/
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &case_comm);

    int row_rank, row_size;
    MPI_Comm_rank(case_comm, &row_rank);
    MPI_Comm_size(case_comm, &row_size);



    //create group
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int *base_ranks = new int[num_cases];
    for (int i = 0; i < num_cases; i++) {
         base_ranks[i] = i;
    }
    MPI_Group case_group;
    MPI_Group_incl(world_group, num_cases,base_ranks, &case_group);

    MPI_Comm_create_group(MPI_COMM_WORLD, case_group, 0, &case_comm_base);




    delete[] base_ranks;

    int new_rank = -1;
    if (case_comm_base != MPI_COMM_NULL) {
        MPI_Comm_rank(case_comm_base, &new_rank);
    }

    // Print out information about the new communicator
/*
   if (new_rank != -1) {
            std::cout << "Process " << world_rank << " (original rank) is in the new communicator. "
                      << "New rank: " << new_rank << std::endl;
        } else {
            std::cout << "Process " << world_rank << " (original rank) is not in the new communicator." << std::endl;
   }
*/


int data=1,result=0;
//MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, case_comm);
/*    if (case_comm_base != MPI_COMM_NULL) {
        MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, case_comm_base);
    }
    std::cout << result << std::endl;*/
//MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);



/*
    if(row_rank==0)MPI_Comm_split(MPI_COMM_WORLD, 0,world_rank, &case_comm_base);
    else MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED,world_rank, &case_comm_base);
int data;
int result;*/
//MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, case_comm);
//MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, case_comm_base);
//MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //teste MPI_Allreduce-Aufruf mit dem neuen Kommunikator
    /*if(case_comm_base != MPI_UNDEFINED) {

            MPI_Allreduce(&data, &result, 1, MPI_INT, MPI_SUM, case_comm_base);




 *//*       int row_rank_base, row_size_base;
        MPI_Comm_rank(case_comm_base, &row_rank_base);
        MPI_Comm_size(case_comm_base, &row_size_base);

        std::cout << "WORLD RANK/SIZE: " << world_rank << " " << world_size << "  case_comm_rank/size " << row_rank << " " << row_size<< "  case_comm_base_rank/size " << row_rank_base << " " << row_size_base << std::endl;
*//*
    }
    //MPI_Comm_free(&case_comm);*/


//#endif

}