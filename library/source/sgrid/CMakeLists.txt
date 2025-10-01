

add_library(
        SGrid


        multiDepthHashGrid.cpp
        multiDepthHashGrid.h
        sparseGrid.cc
        sparseGrid.h

        depth.h
        ListOfDepthOrderedGrids.h

        ListOfDepthOrderedGrids.cc
        ListOfDepthOrderedGrids.h

        localTensorProduct.cc

        komponente.cpp
        komponente.h

        multilevelSparseGrid.h
        multilevelSparseGrid.cpp

        localTensorProduct.h


        simpleMultiHash.h



        GridGeneration.h
        GridGeneration.cpp

        refinement.h
        refinement.cpp
)




if (MPI_FOUND)
    target_link_libraries(SGrid MPI::MPI_CXX OpenMP::OpenMP_CXX Indices Iterator)
else (MPI_FOUND)
    target_link_libraries(SGrid OpenMP::OpenMP_CXX Indices Iterator)
endif (MPI_FOUND)