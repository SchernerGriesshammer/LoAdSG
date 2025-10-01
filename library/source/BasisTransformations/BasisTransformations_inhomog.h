//
// Created by scherner on 29.04.21.
//

#ifndef GRUN_BASISTRANSFORMATIONS_INHOMOG_H
#define GRUN_BASISTRANSFORMATIONS_INHOMOG_H

#include "BasisTransformations.h"
#include "../extemp/vector.h"
#include "../sgrid/komponente.h"
#include "../extemp/shift.h"
#include "../sgoperation/matrix.h"
#include "../extemp/multilevelvector.h"
#include "../tests/old_versions/MatrixVectorMultiplicationPrewavelets.h"


void CalcNodalBasis_Neumann(VectorSparseG &u, ListOfDepthOrderedSubgrids &orderedSubgrids);

void CalcNodalBasis_inhomog(VectorSparseG &u, Depth &Tiefe);

void
CalcHierarchicalBasis_Neumann(VectorSparseG &u, const Depth &Tsubgrid, ListOfDepthOrderedSubgrids &orderedSubgrids);

void CalcHierarchicalBasis_Neumann(VectorSparseG &u, ListOfDepthOrderedSubgrids &orderedSubgrids);

void CalcHierbyPrew_Neumann(VectorSparseG &prew, VectorSparseG &hier, ListOfDepthOrderedSubgrids &allgrids,
                            SubgridFixedDepth &subgrid);

void CalcUbyPrew_Neumann(VectorSparseG &prew, VectorSparseG &u);

/**
 * Solve on component for prewavelet coefficients and save them in prew.
 * @param prew
 * @param u
 * @param komponente
 */
void SolveSave_Neumann(VectorSparseG *prew, VectorSparseG *u, ZusammenhangsKomponente &komponente);

void CalcPrew_inhomogen(VectorSparseG &prew, VectorSparseG &u);

/**
 * Welche prewavelet Koeffizienten werden hier genommen?
 *
 * @param prew
 * @param u
 * @param Tiefe
 * @param restrictions
 *
 *
 */
void CalcUbyPrewRestrictions_inhomog(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions);

void CalcUbyPrewRestrictions_inhomog2(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions,
                                      MultiLevelVector &nodal, ListOfDepthOrderedSubgrids &list);

void CalcUbyPrewRestrictions_inhomog3(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions,
                                      MultiLevelVector &nodal);

void CalcUbyPrew_inhomog(VectorSparseG &prew, VectorSparseG &u, Depth &Tiefe);

#endif //GRUN_BASISTRANSFORMATIONS_INHOMOG_H
