//
// Created by to35jepo on 2/16/23.
//

#include "PoissonStencil.h"

void PoissonStencil::initialize(Depth &T_) {


    for (int d = 0; d < DimensionSparseGrid; d++) {
        int exp = T_.at(d);
        if(exp == 0) meshwidth[d] = 1.0;
        else {

            int value2 = 1 << exp;
            meshwidth[d] = 1.0 / double(value2);

        }



    }

};