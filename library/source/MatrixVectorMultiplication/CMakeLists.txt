add_library(
        MatrixVectorMultiplication
        MatrixVectorInhomogen.h
        MatrixVectorInhomogen.cpp
        MatrixVectorHomogen.h
        MatrixVectorHomogen.cpp
        RHS.h
        RHS.cpp
)

if (MPI_FOUND)
    target_link_libraries(MatrixVectorMultiplication MPI::MPI_CXX OpenMP::OpenMP_CXX)
else (MPI_FOUND)
    target_link_libraries(MatrixVectorMultiplication OpenMP::OpenMP_CXX)
endif (MPI_FOUND)
