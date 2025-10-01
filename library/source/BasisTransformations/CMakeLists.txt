
add_library(
       BasisTransformations

       BasisTransformations.cpp
       BasisTransformations.h

       BasisTransformations_inhomog.cpp
       BasisTransformations_inhomog.h


)

if (MPI_FOUND)
    target_link_libraries(BasisTransformations MPI::MPI_CXX   SGoperation SGrid Iterator)
else (MPI_FOUND)
    target_link_libraries(BasisTransformations  SGoperation SGrid Iterator)
endif (MPI_FOUND)

