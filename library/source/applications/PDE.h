//
// Created by to35jepo on 9/19/23.
//

#ifndef RUN_PDE_H
#define RUN_PDE_H


#include "../sgrid/sparseGrid.h"
#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "../MatrixVectorMultiplication/RHS.h"
#include "../stencils/Stencil.h"



class PDE_Poisson_InHomRHS{
public:
    PDE_Poisson_InHomRHS(AdaptiveSparseGrid& n_grid, AdaptiveSparseGrid& d_grid, MatrixVectorHomogen& matrix_):neumann_grid(n_grid),dirichlet_grid(d_grid), f_dirichlet(d_grid), prew(d_grid),matrix(matrix_),dirichlet(d_grid){};
    void SetUp(VectorSparseG& f_neumann, VectorSparseG& g_neumann);

    VectorSparseG& getRHS(){
        return f_dirichlet;
    }
private:
    AdaptiveSparseGrid& neumann_grid;
    AdaptiveSparseGrid& dirichlet_grid;
    VectorSparseG f_dirichlet;
    VectorSparseG prew;
    MatrixVectorHomogen& matrix;
    VectorSparseG dirichlet;
};

/*
double var_coeff2( const std::array<double, DimensionSparseGrid>& x )
{
    double val = 1.0;
    for(size_t d= 0;d<DimensionSparseGrid; d++ ){
        double t = x.at(d);
        val *= (1.0-(t*t));
    }

    return val;

}
*/


class PDE_Helmholtz_VariableCoeff{
public:
    PDE_Helmholtz_VariableCoeff(AdaptiveSparseGrid& n_grid, AdaptiveSparseGrid& d_grid, MatrixVectorHomogen& matrix_):neumann_grid(n_grid),dirichlet_grid(d_grid), f_dirichlet(d_grid), prew(d_grid),matrix(matrix_),dirichlet(d_grid){};
    /**
 * @f[
 * -\Delta u + c(x) u = f \quad \text{in}\; \Omega = (0,1)^d \\
 * u = g \quad \text{on}\; \partial \Omega
 * @f]
 * transformed to
 * @f[
 * -\Delta (u-g) + c(x) (u-g) = f+\Delta g - c(x)g \quad \text{in}\; \Omega = (0,1)^d \\
 * u - g= 0 \quad \text{on}\; \partial \Omega
 * @f]
 *
 * @param f_neumann
 * @param g_neumann
 */
    void SetUp(VectorSparseG& f_neumann, VectorSparseG& g_neumann) {

        CalcHierarchicalBasis(f_neumann);


        for (unsigned long i = 0; i < neumann_grid.getMaximalOccupiedSecondTable(); i++) {
            if (neumann_grid.getActiveTable()[i]) {
                IndexDimension I = neumann_grid.getIndexOfTable(i);
                Depth T(I);
                if((T>>0)){
                    f_dirichlet.setValue(I,f_neumann.getValue(i));
                    f_neumann.setValue(i, 0.0);
                }
            }
        }




        CalcNodalBasis(f_dirichlet);
        calcPrewByNodal(prew,f_dirichlet);
        //MassStencil stencil1;
        HelmHoltz stencil1(*f_dirichlet.getSparseGrid());

        f_dirichlet=0.0;
        matrix.multiplication(prew,f_dirichlet,stencil1);


        InHomoBoundaryRHSHierarchical test_hier(dirichlet_grid);
        test_hier.multiply_mass(f_neumann, dirichlet);

        f_dirichlet=f_dirichlet+dirichlet;

        CalcHierarchicalBasisForRHS(g_neumann);
        dirichlet = 0.0;
        test_hier.multiply(g_neumann,dirichlet);

        f_dirichlet=f_dirichlet-dirichlet;

        dirichlet=0.0;

        test_hier.multiply_mass_coeff(g_neumann,dirichlet);




        f_dirichlet = f_dirichlet-dirichlet;



    }


    VectorSparseG& getRHS(){
        return f_dirichlet;
    }
private:
    AdaptiveSparseGrid& neumann_grid;
    AdaptiveSparseGrid& dirichlet_grid;
    VectorSparseG f_dirichlet;
    VectorSparseG prew;
    MatrixVectorHomogen& matrix;
    VectorSparseG dirichlet;
};



class PDE_Helmholtz_VariableCoeff_interface{
public:
    PDE_Helmholtz_VariableCoeff_interface(AdaptiveSparseGrid& n_grid, AdaptiveSparseGrid& d_grid, MatrixVectorHomogen& matrix_):neumann_grid(n_grid),dirichlet_grid(d_grid), f_dirichlet(d_grid), prew(d_grid),matrix(matrix_),dirichlet(d_grid){};
    /**
 * @f[
 * -\Delta u + c(x) u = f \quad \text{in}\; \Omega = (0,1)^d \\
 * u = g \quad \text{on}\; \partial \Omega
 * @f]
 * transformed to
 * @f[
 * -\Delta (u-g) + c(x) (u-g) = f+\Delta g - c(x)g \quad \text{in}\; \Omega = (0,1)^d \\
 * u - g= 0 \quad \text{on}\; \partial \Omega
 * @f]
 *
 * @param f_neumann
 * @param g_neumann
 */
    template<class F>
    void SetUp(VectorSparseG& f_neumann, VectorSparseG& g_neumann,const F& f) {

        CalcHierarchicalBasis(f_neumann);


        for (unsigned long i = 0; i < neumann_grid.getMaximalOccupiedSecondTable(); i++) {
            if (neumann_grid.getActiveTable()[i]) {
                IndexDimension I = neumann_grid.getIndexOfTable(i);
                Depth T(I);
                if((T>>0)){
                    f_dirichlet.setValue(I,f_neumann.getValue(i));
                    f_neumann.setValue(i, 0.0);
                }
            }
        }




        CalcNodalBasis(f_dirichlet);
        calcPrewByNodal(prew,f_dirichlet);
        //MassStencil stencil1;
        HelmHoltz stencil1(f_dirichlet.getSparseGrid());


        f_dirichlet=0.0;
        matrix.multiplication(prew,f_dirichlet,stencil1);


        InHomoBoundaryRHSHierarchical test_hier(dirichlet_grid);
        test_hier.multiply_mass(f_neumann, dirichlet);

        f_dirichlet=f_dirichlet+dirichlet;

        CalcHierarchicalBasisForRHS(g_neumann);
        dirichlet = 0.0;
        test_hier.multiply(g_neumann,dirichlet);

        f_dirichlet=f_dirichlet-dirichlet;



        dirichlet=0.0;
        test_hier.multiply_mass_coeff_interface(g_neumann,dirichlet,f);



        f_dirichlet = f_dirichlet-dirichlet;


    }


    VectorSparseG& getRHS(){
        return f_dirichlet;
    }
private:
    AdaptiveSparseGrid& neumann_grid;
    AdaptiveSparseGrid& dirichlet_grid;
    VectorSparseG f_dirichlet;
    VectorSparseG prew;
    MatrixVectorHomogen& matrix;
    VectorSparseG dirichlet;
};
#endif //RUN_PDE_H
