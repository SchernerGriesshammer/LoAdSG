add_library(
       LocalStiffnessMatrices


      LocalStiffnessMatrices.h
      LocalStiffnessMatrices.cpp
        ADDLocalStiffnessMatrices.cpp
        ADDLocalStiffnessMatrices.h
        LocalStiffnessMatricesDynamic.cpp
        LocalStiffnessMatricesDynamicDistribution.h
        LocalStiffnessClassMemoryDistribution.cpp
        LocalStiffnessMatricesMemoryDistribution.h
        LocalStiffnessMatricesFixedDistribution.cpp
        LocalStiffnessMatricesFixedDistribution.h


)

if (MPI_FOUND)
    target_link_libraries(LocalStiffnessMatrices MPI::MPI_CXX OpenMP::OpenMP_CXX Extemp)
else (MPI_FOUND)
    target_link_libraries(LocalStiffnessMatrices OpenMP::OpenMP_CXX Extemp)
endif (MPI_FOUND)
