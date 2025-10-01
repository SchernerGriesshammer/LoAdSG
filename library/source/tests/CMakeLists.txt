
add_library(
       Test


        old_versions/MatrixVectorMultiplicationPrewavelets.cpp
        old_versions/MatrixVectorMultiplicationPrewavelets.h


        testing.cpp
        testing.h )

if (MPI_FOUND)
    target_link_libraries(Test MPI::MPI_CXX)
else (MPI_FOUND)
    target_link_libraries(Test)
endif (MPI_FOUND)

