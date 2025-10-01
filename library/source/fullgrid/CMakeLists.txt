

add_library(
       Fullgrid
       FullGrid.h
       FullGrid.cpp
       MatrixVectorMultiplication_Fullgrid.h
       MatrixVectorMultiplication_Fullgrid.cpp)


if (MPI_FOUND)
target_link_libraries(Fullgrid MPI::MPI_CXX)
else(MPI_FOUND)
    target_link_libraries(Fullgrid)
endif(MPI_FOUND)

