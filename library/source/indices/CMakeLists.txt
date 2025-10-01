


add_library(
        Indices
        index.h
        index.cc
)

if (MPI_FOUND)
    target_link_libraries(Indices MPI::MPI_CXX )
else(MPI_FOUND)
    target_link_libraries(Indices )
endif(MPI_FOUND)

