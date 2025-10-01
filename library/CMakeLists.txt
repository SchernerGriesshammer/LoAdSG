cmake_minimum_required(VERSION 3.5)
project(run)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")


# Include directories
#include_directories("/local/to35jepo/SGpp/base/src")

# Link directories
#link_directories("/local/to35jepo/SGpp/lib")

if(MPI_ON)
    find_package(MPI REQUIRED)
else(MPI_ON)
    set(MPI_FOUND false)
    set(MY_MPI_ON false)
endif(MPI_ON)

if(OPENMP_ON)
    find_package(OpenMP REQUIRED)
else(OPENMP_ON)
    set(OPENMP_FOUND false)
endif(OPENMP_ON)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (MPI_FOUND)
    message("This application will be compiled with MPI")
    add_compile_definitions(_MPI_PARALLEL)
    add_compile_definitions(MY_MPI_ON)
else (MPI_FOUND)
    message("This application will be compiled without MPI")
endif (MPI_FOUND)


if(MAIN_FILE)
    add_executable(${PROJECT_NAME} ${MAIN_FILE})
else(MAIN_FILE)
    add_executable(${PROJECT_NAME} "helmholtz_variable_coeff.cpp")
endif(MAIN_FILE)




if (MPI_FOUND)
     target_link_libraries(${PROJECT_NAME} MPI::MPI_CXX)
endif (MPI_FOUND)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("This application will be compiled with OMP")
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()



add_subdirectory(source)


#target_link_libraries(${PROJECT_NAME} Source sgppbase)
target_link_libraries(${PROJECT_NAME} Source)


# Set the location of the dynamically linked libraries
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the path to LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)




