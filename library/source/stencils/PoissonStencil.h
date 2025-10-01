//
// Created by to35jepo on 2/16/23.
//

#ifndef SGRUN_POISSONSTENCIL_H
#define SGRUN_POISSONSTENCIL_H

#include "../extemp/vector.h"
#include "../myMath/myMath.h"
#include "../cells/celldimension.h"
enum TypeMatrixVectorMultiplication  { StencilOnTheFly, StoreLocalStiffnessMatrix};


class PoissonStencil {
public:
    virtual inline void applyStencilOnCell_MPI_OMP(CellDimension& cell,VectorSparseG& input, VectorSparseG& output) {
        cout << " PoissonStencil::applyStencilOnCell_MPI_OMP NOT IMPLEMENTED YET!" << endl;
        exit(1);
    }
    virtual inline void applyStencilOnCell_MPI_OMP_nosymmetry(CellDimension& cell,VectorSparseG& input, VectorSparseG& output) {
        cout << " PoissonStencil::applyStencilOnCell_MPI_OMP_nosymmetry NOT IMPLEMENTED YET!" << endl;
        exit(1);
    }
    virtual inline void applyStencilOnCell(CellDimension& cell,VectorSparseG& input, VectorSparseG& output){
        cout << " PoissonStencil::applyStencilOnCell NOT IMPLEMENTED YET!" << endl;
        exit(1);
    }
    virtual inline void applyStencilOnCell_BinaryIterator(CellDimension& cell,VectorSparseG& input, VectorSparseG& output){
        cout << " PoissonStencil::applyStencilOnCell_BinaryIterator NOT IMPLEMENTED YET!" << endl;
        exit(1);
    }
    void initialize(Depth &T_);
    inline double returnValue(const IndexDimension &Index,const MultiDimCompass &mc) const {

        double factors[DimensionSparseGrid];
        for (int j = 0; j < DimensionSparseGrid; ++j) {
            double h = meshwidth[j];
            Richtung r = mc.getRichtung(j);
            factors[j] = (1.0 / 6.0) * h;
            if(r == Mitte) {
                if((!Index.isAtLeftBoundary(j)) && (!Index.isAtRightBoundary(j)))factors[j] *= 4.0;
                else factors[j] *= 2.0;
            }
        }

        // Stencil == Stiffness
        double sum = 0.0;
        for (int i = 0; i < DimensionSparseGrid; i++) {
            // int du/dxi * dv/dxi  d(x1...xd)
            double prod = 1.0;
            for (int j = 0; j < DimensionSparseGrid; ++j) {
                // integral factor in direction j
                if (i == j) {
                    double h = meshwidth[j];
                    Richtung r = mc.getRichtung(j);
                    // factor: int du/dxi * dv/dxi dxi
                    prod *= (1.0 / h);
                    if (r != Mitte) {
                        prod *= -1.0;
                    } else if((!Index.isAtLeftBoundary(j)) && (!Index.isAtRightBoundary(j))){
                        prod *= 2.0;
                    }
                } else {
                    // factor int du/dxi*dv/dxi dxj = int u*v*xj
                    prod *= factors[j];
                }
            }
            sum += prod;
        }
        return sum;
    }

    TypeMatrixVectorMultiplication getTypeMatrixVectorMultiplication(){return typeMatrixVectorMultiplication;};


private:
    double meshwidth[DimensionSparseGrid];

    TypeMatrixVectorMultiplication typeMatrixVectorMultiplication = StencilOnTheFly;
};


#endif //SGRUN_POISSONSTENCIL_H
