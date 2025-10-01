


add_library(
        Extemp

        coordinate.cc
        coordinate.h
        extempAlg.h
        funktor.h
        funktorExample.cc
        funktorExample.h

        shift.cc
        shift.h
        system.h
        vector.cc
        vector.h

        multilevelvector.cpp multilevelvector.h)

if (MPI_FOUND)
target_link_libraries(Extemp MPI::MPI_CXX )
else(MPI_FOUND)
    target_link_libraries(Extemp )
endif(MPI_FOUND)

