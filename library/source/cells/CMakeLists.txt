
#find_package(MPI REQUIRED)

add_library(
        cells

        celldimension.h
        celldimension.cpp


        MultiDepthHashCellStructure.cpp
        MultiDepthHashCellStructure.h
        iterator_neu.cpp
        iterator_neu.h
        SingleDepthCellStructure.cpp
        SingleDepthCellStructure.h
        CellStructure.cpp
        CellStructure.h
)


if (MPI_FOUND)
target_link_libraries(cells MPI::MPI_CXX )
else(MPI_FOUND)
    target_link_libraries(cells )
endif(MPI_FOUND)

