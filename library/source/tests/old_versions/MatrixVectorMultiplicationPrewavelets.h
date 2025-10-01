//
// Created by scherner on 22.04.21.
//

#ifndef GRUN_MATRIXVECTORMULTIPLICATIONPREWAVELETS_H
#define GRUN_MATRIXVECTORMULTIPLICATIONPREWAVELETS_H

#include "../../indices/index.h"
#include "../../sgrid/depth.h"
#include <cmath>
#include "../../extemp/vector.h"
#include "../../extemp/multilevelvector.h"
#include "../../sgrid/multilevelSparseGrid.h"
#include "../../mympi.h"


#include "../../sgoperation/matrix.h"
#include "../../BasisTransformations/BasisTransformations_inhomog.h"
#include "../../BasisTransformations/BasisTransformations.h"



class CasesIterator {
public:
    CasesIterator() { snum = 0; }

    inline bool goon() { return (snum < pow(2, DimensionSparseGrid)); }

    inline void operator++() {
        snum++;
    }


    inline bool *getcase() {
        static bool cases[DimensionSparseGrid];


        for (int i = 0; i < DimensionSparseGrid; ++i) cases[i] = false;

        for (int i = 0; i < DimensionSparseGrid; ++i) {
            // the following condition checks
            // if the ith bit of snum in binary form
            // is 1 or not
            if ((snum & (1 << i)) != 0) {

                cases[i] = true;
            }
        }
        return cases;
    };

    bool *getcase(int snum_);

    int getCaseNumber() { return snum; };

    /* int getCaseNumber(bool* case_){
         int returnvalue;
         returnvalue = 0.0;
         int b = 1;
         for(int d=DimensionSparseGrid-1; d >=0; d--) {
             returnvalue = returnvalue + case_[d] * b;
             b = b * 2;
         }

     }*/





private:

    int snum;
};

enum StencilType {
    Stiffness, Mass
};

//////////////////// Dirichlet Matrix-Vektor-Multiplikation
void MatrixVectorPrewavelet(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, MultiLevelVector &gM);

//Functionen, die man für MaxtrixVectorPrewavelet braucht

//alte Variante ohne MultigridVector
/*void PP_Alg(VectorSparseG &prew, VectorSparseG &u, VectorSparseG &Ax, StencilType type, bool *restrictions);*/


void MatrixVector_Case(VectorSparseG &prew, VectorSparseG &u, VectorSparseG &Ax, StencilType type, bool *restrictions,
                       VectorSparseG &g, VectorSparseG &z, MultiLevelVector &gM, ListOfDepthOrderedSubgrids &list,
                       MultiLevelVector &nodal);

void ApplyStencil(VectorSparseG &x, VectorSparseG &Ax, Depth &T, StencilType type);

double CalcStencilValue(IndexDimension Index, Depth &T, VectorSparseG &u, StencilType type);


IndexDimension GetStencilValueIndex(IndexDimension Index, MultiDimCompass mc, double *stencilvalue, Depth T, int dir,
                                    StencilType type);


void restriction(VectorSparseG &fine, VectorSparseG &coarse, int t, int d);

void ConvertToPrewavelet(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T, SubgridFixedDepth &subgrid);

void ConvertToPrewavelet2(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T);

void ConvertDualNodalToPrewavelet(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T);
////////////// MV Neumann
void MatrixVectorPrewavelet_inhomog(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, MultiLevelVector &gM);

void test_neumann(VectorSparseG &prew, VectorSparseG &Ax, StencilType type, bool *restrictions, MultiLevelVector &gM);


class Stencil {
public:
    Stencil(Depth T_, StencilType type_);

    double returnValue(MultiDimCompass &mc, IndexDimension Index);

private:
    double values[MaxShift<DimensionSparseGrid>::value];
    StencilType type;
    Depth T;

};

void ApplyStencil_inhomog(VectorSparseG &x, VectorSparseG &Ax, Depth &T, StencilType type);

double CalcStencilValue_Boundary(IndexDimension Index, Depth &T, VectorSparseG &u, StencilType type, Stencil stencil);

IndexDimension
GetStencilValueIndex_Boundary(IndexDimension Index, MultiDimCompass mc, double *stencilvalue, Depth T, int dir,
                              StencilType type);

double
GetStencilValue_Boundary(MultiDimCompass mc, Depth T, int dir,
                         StencilType type);

bool getNextIndex(IndexDimension Index, MultiDimCompass mc, Depth T, IndexDimension &NextIndex );

void restriction_neumann(VectorSparseG &fine, VectorSparseG &coarse, int t, int d);

void ConvertToNeumannPrewavelet(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T);

void ConvertToNeumannPrewavelet2(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T, ListOfDepthOrderedSubgrids &list);

////////////////////// volle Gitter

void ApplyStencil(VectorSparseG &x, VectorSparseG &Ax, StencilType type);

void restriction_local(VectorSparseG &fine, VectorSparseG &coarse, Depth Tfine, Depth Tcoarse);

void matrixmult_prew_fullgrid(VectorSparseG &x, VectorSparseG &Ax, StencilType type);






///////////////////////// Alte DualSpace Transformation
/*void ConvertDualSpace(VectorSparseG &neumann, VectorSparseG &dirichlet);*/

#endif //GRUN_MATRIXVECTORMULTIPLICATIONPREWAVELETS_H
