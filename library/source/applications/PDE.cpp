//
// Created by to35jepo on 9/19/23.
//

#include "PDE.h"
#include "../stencils/MassStencil.h"
#include "../MatrixVectorMultiplication/RHS.h"

void PDE_Poisson_InHomRHS::SetUp(VectorSparseG& f_neumann, VectorSparseG& g_neumann){
// - lap(u)=f, u=g on boundary and f also defined on boundary --> set up RHS for -lap(u-g)=-lap(u)+lap(g)=f+lap(g)=f_inner+f_outer+lap(g)
// split f into inner and boundary basis
    CalcHierarchicalBasis(f_neumann);
    for (unsigned long i = 0; i < neumann_grid.getMaximalOccupiedSecondTable(); i++) {
        if (neumann_grid.getActiveTable()[i]) {
            IndexDimension I = neumann_grid.getIndexOfTable(i);
            Depth T(I);
            if((T>>0)){
                f_dirichlet.setValue(I,f_neumann.getValue(i));
                f_neumann.setValue(i,0.0);
            }
        }
    }




    // f_dirichlet only defined in inner grid points -> use standard matrix vector multiplication
    // calculates f_inner*phi_p
    CalcNodalBasis(f_dirichlet);
    calcPrewByNodal(prew,f_dirichlet);
    //MassStencil stencil1;
    HelmHoltz stencil1(*prew.getSparseGrid());


    f_dirichlet=0.0;
    matrix.multiplication<HelmHoltz>(prew,f_dirichlet,stencil1);




    // calculates f_outer*phi_p
    InHomoBoundaryRHSHierarchical test_hier(dirichlet_grid);
    test_hier.multiply_mass(f_neumann,dirichlet);

    f_dirichlet=f_dirichlet+dirichlet;

    // transforms g_neumann (given only on boundary points) into hierarchical basis
    // calculates then grad(g_neumann)*grad(phi)
    CalcHierarchicalBasisForRHS(g_neumann);
    dirichlet = 0.0;
    test_hier.multiply(g_neumann,dirichlet);
    f_dirichlet=f_dirichlet-dirichlet;
};

