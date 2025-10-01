//
// Created by to35jepo on 2/16/23.
//

#ifndef SGRUN_MASSSTENCIL_H
#define SGRUN_MASSSTENCIL_H


#include "../extemp/vector.h"
#include "../myMath/myMath.h"
#include "../cells/celldimension.h"

class MassStencil {
public:

     inline void initialize(Depth &T_) {

        for (int d = 0; d < DimensionSparseGrid; d++) {
            auto value = double(POW2( T_.at(d)));
            meshwidth[d] = 1.0 / value;
        }

    };
    virtual inline void applyStencilOnCell_MPI_OMP(CellDimension& cell,VectorSparseG& input, VectorSparseG& output) {
        cout << " MassStencil applyStencilOnCell NOT IMPLEMENTED YET!" << endl;
        exit(1);
    }
    virtual inline void applyStencilOnCell(CellDimension& cell,VectorSparseG& input, VectorSparseG& output){
        cout << " MassStencil applyStencilOnCell NOT IMPLEMENTED YET!" << endl;
        exit(1);
       /* for(CellIndexIterator outerIter(&cell); outerIter.goon(); ++outerIter){

            IndexDimension p = outerIter.getIndex();
            CellIndexDirection dirP = outerIter.getCellIndexDirection();


            for(CellIndexIterator innerIter(outerIter); innerIter.goon(); ++innerIter){
                IndexDimension q=innerIter.getIndex();
                CellIndexDirection dirQ = innerIter.getCellIndexDirection();
                unsigned long k;
                if(output.getSparseGrid()->occupied(k,q)){

                    double val = integration(cell, dirP, dirQ);
                    output.setValue(k,val);
                }


            }
        }*/
    }

    inline double returnValue(const IndexDimension &Index, MultiDimCompass &mc) const {
        double factors[DimensionSparseGrid];
        for (int j = 0; j < DimensionSparseGrid; ++j) {
            double h = meshwidth[j];
            Richtung r = mc.getRichtung(j);
            factors[j] = (1.0 / 6.0) * h;;
            if(r == Mitte) {
                if((!Index.isAtLeftBoundary(j)) && (!Index.isAtRightBoundary(j)))factors[j] *= 4;
                else factors[j] *= 2;
            }
            if(r == Links && Index.isAtLeftBoundary(j)) factors[j]=1.0;
            if (r == Rechts && Index.isAtRightBoundary(j)) factors[j]=1.0;

        }



        double prod = 1.0;
        for (int j = 0; j < DimensionSparseGrid; j++) {

            prod *= factors[j];


        }





        return prod;

    }
    // inline double returnValue(IndexDimension &Index, MultiDimCompass &mc) const {
    //     size_t factors[DimensionSparseGrid];
    //     for (int j = 0; j < DimensionSparseGrid; ++j) {
    //         Richtung r = mc.getRichtung(j);
    //         factors[j] = 0;
    //         if(r == Mitte) {
    //             factors[j] += 1;
    //             if((!Index.isAtLeftBoundary(j)) && (!Index.isAtRightBoundary(j)))
    //             factors[j] += 1;
    //         }
    //     }

    //     // Stencil == Stiffness
    //     double sum = 0.0;
    //     for (int i = 0; i < DimensionSparseGrid; i++) {
    //         // int du/dxi * dv/dxi  d(x1...xd)
    //         double prod = 1.0;
    //         for (int j = 0; j < DimensionSparseGrid; ++j) {
    //             // integral factor in direction j

    //             double h = meshwidth[j];

    //             if (i == j) {
    //                 // factor: int du/dxi * dv/dxi dxi
    //                 prod *= (1.0 / h);
    //                 if (factors[j] == 0) {
    //                     prod *= -1;
    //                 } else if(factors[j] == 2){
    //                     prod *= 2.0;
    //                 }
    //             } else {
    //                 // factor int du/dxi*dv/dxi dxj = int u*v*xj
    //                 prod *= (1.0 / 6.0) * h;
    //                 if (factors[j] == 1) {
    //                     prod *= 2.0;
    //                 } else if (factors[j] == 2) {
    //                     prod *= 4.0;
    //                 }

    //             }
    //         }
    //         sum += prod;
    //     }


    //     return sum;


    // }


private:
    double meshwidth[DimensionSparseGrid];

};


#endif //SGRUN_MASSSTENCIL_H
