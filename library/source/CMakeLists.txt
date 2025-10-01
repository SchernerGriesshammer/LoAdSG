


add_library(

        Source

        primes/prime.h
        primes/prime.cc

        myAssert.h
        mympi.h
        abbrevi.h


        myMath/myMath.h



        mycomm.h
        mycomm.cpp
        functions/functions.cpp
        functions/functions.h


)


add_subdirectory(extemp)
add_subdirectory(sgoperation)
add_subdirectory(sgrid)
add_subdirectory(indices)
add_subdirectory(applications)
add_subdirectory(MatrixVectorMultiplication)
add_subdirectory(fullgrid)

add_subdirectory(iterator)
add_subdirectory(BasisTransformations)
add_subdirectory(tests)
add_subdirectory(cells)
add_subdirectory(stencils)
add_subdirectory(localStiffnessMatrices)




if (MPI_FOUND)
    target_link_libraries(Source MPI::MPI_CXX OpenMP::OpenMP_CXX  Extemp SGrid SGoperation  Indices  Applications MatrixVectorMultiplication LocalStiffnessMatrices Stencils Fullgrid Iterator BasisTransformations Test cells )
else (MPI_FOUND)
    target_link_libraries(Source OpenMP::OpenMP_CXX Extemp SGrid SGoperation Indices Applications MatrixVectorMultiplication LocalStiffnessMatrices Stencils Fullgrid Iterator BasisTransformations Test cells )
endif (MPI_FOUND)


 

