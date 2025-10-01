add_library(
        Iterator
        RectangularIterator.h
        RectangularIterator.cpp

        depthIterator.h
        depthIterator.cpp
       )

if (MPI_FOUND)
    target_link_libraries(Iterator MPI::MPI_CXX OpenMP::OpenMP_CXX)
else (MPI_FOUND)
    target_link_libraries(Iterator OpenMP::OpenMP_CXX)
endif (MPI_FOUND)
