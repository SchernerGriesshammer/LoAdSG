
#find_package(MPI REQUIRED)

add_library(
        SGoperation



        matrix_operations.cc
        matrix_operations.h


        matrix.cpp
        matrix.h

)


if (MPI_FOUND)
target_link_libraries(SGoperation MPI::MPI_CXX )
else(MPI_FOUND)
    target_link_libraries(SGoperation )
endif(MPI_FOUND)

