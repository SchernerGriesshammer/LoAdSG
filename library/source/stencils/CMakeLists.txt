add_library(
        Stencils

        Stencil.h

        PoissonStencil.h
        PoissonStencil.cpp

        MassStencil.h
        MassStencil.cpp




        InterfaceIntegration/interfaceMatrices.h
        InterfaceIntegration/constantIntegrators.h
        InterfaceIntegration/interatorBasisFunction.h
        InterfaceIntegration/HelmholtzIntegrator/helmIntegrator.h
        StencilMC.cpp
        StencilMC.h


)

if (MPI_FOUND)
    target_link_libraries(Stencils MPI::MPI_CXX OpenMP::OpenMP_CXX)
else (MPI_FOUND)
    target_link_libraries(Stencils OpenMP::OpenMP_CXX)
endif (MPI_FOUND)
