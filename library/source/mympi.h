/**********************************************************************************
*
 **********************************************************************************/

//////////////////////////////////////////////////////////
//   dummympi.h  (include file for library)
/////////////////////////////////////////////////////////

#ifndef RUN_MPI
#define RUN_MPI
enum Type_of_task { TAG_TASK=0, TAG_TERMINATE=1,TAG_LOCALSTIFFNESS_SEND=2,  TAG_LOCALSTIFFNESS_SENDBACK=3,TAG_LOCALSTIFFNESS_END=4 };


#ifdef _MPI_PARALLEL
#include <mpi.h>


#endif
#include <iostream>





#ifndef _MPI_PARALLEL


#ifndef DUMMYMPI_H_
#define DUMMYMPI_H_


#define MPI_MAX_PROCESSOR_NAME 6

inline int MPI_Get_processor_name(char s[MPI_MAX_PROCESSOR_NAME], int *i) {
    return 0;
}

#define MPI_Comm         int
#define MPI_Group         int
#define MPI_COMM_WORLD   0
#define MPI_Request      int
#define MPI_Datatype     int
#define MPI_IN_PLACE     double
#define MPI_UNDEFINED    int
#define MPI_COMM_NULL    0
#define MPI_Status  int



#define MPI_INT        0
#define MPI_DOUBLE     0
#define MPI_LOR        0
#define MPI_SUM        0
#define MPI_MIN        0
#define MPI_BOOL       0
#define MPI_OR         0
#define MPI_DOUBLE_COMPLEX        0
#define MPI_BYTE        0
#define MPI_MAX        0
#define MPI_ANY_SOURCE 0

inline void Error_Dummy_mpi() {
    std::cout << " Error_Dummy_mpi! " << std::endl;
}


inline void MPI_Comm_rank(int, int *my_rank) { *my_rank = 0; }

inline void MPI_Comm_size(int, int *p) { *p = 1; }

template<class A>
inline void MPI_Irecv(A *, int, int, int, int, int, int *) {}

template<class A>
inline void MPI_Isend(A *, int, int, int, int, int, int *) {}

template<class A>
inline void MPI_Recv(A *, int, int, int, int, int) {}

template<class A>
inline void MPI_Send(A *, int, int, int, int, int) {}

template<class A>
inline void MPI_Recv(A *, int, int, int, int, int, int *) {}

template<class A>
inline void MPI_Send(A *, int, int, int, int, int, int *) {}

inline void MPI_Waitall(int, int *, int *) {}

inline void MPI_Barrier(int) {}

inline void MPI_Gather(int *, int, int, int *, int, int, int, int) {}

inline void MPI_Type_commit(int *) {}

inline void MPI_Type_free(int *) {}

inline void MPI_Type_vector(int, int, int, int, int *) {}

template<class A>
inline void MPI_Bcast(A *, int, int, int, int *) {}

template<class A>
inline void MPI_Reduce(A *, A *, int, int, int, int, int *) {}

inline int MPI_Init(int *, char ***) {
    return 0;
}

inline double MPI_Wtime() {
    return clock() / CLOCKS_PER_SEC;
}

inline void MPI_Allreduce(int *input, int *output, int, int, int, int) {
    *output = *input;
}

inline void MPI_Allreduce(double *input, double *output, int, int, int, int) {
    *output = *input;
}

inline void MPI_Allreduce(bool *input, bool *output, int, int, int, int) {
    *output = *input;
}

inline int MPI_Finalize() {
    std::cout << " This is the serial version of EXPDE on sparse grids! " << std::endl;
    return 0;
}
inline int MPI_Comm_split(MPI_Comm comm,int color,int key,MPI_Comm *My_Comm);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm * newcomm);

inline int MPI_Comm_free(MPI_Comm *My_Comm);
#endif
#endif
#endif