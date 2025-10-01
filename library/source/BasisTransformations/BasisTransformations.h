//
// Created by scherner on 12.03.21.
//

#ifndef GRUN_BASISTRANSFORMATIONS_H
#define GRUN_BASISTRANSFORMATIONS_H


#include "../sgrid/komponente.h"
#include "../sgoperation/matrix_operations.h"
#include "../sgoperation/matrix.h"
#include "../extemp/shift.h"
#include "../sgrid/komponente.h"
#include "../sgrid/depth.h"
#include "../tests/old_versions/MatrixVectorMultiplicationPrewavelets.h"


/**
 * berechnet hierarchische Basis auf Tiefen kleiner gleich Tsubgrid
 **/
void CalcHierarchicalBasis(VectorSparseG &vectorSparseG, const Depth &Tsubgrid, ListOfDepthOrderedSubgrids &orderedSubgrids);


/**
 * berechnet hierarchische Basis auf Tiefen kleiner gleich Tsubgrid
 **/
void CalcHierarchicalBasis_NEU(VectorSparseG &vectorSparseG, const Depth &Tsubgrid, std::list<Depth> SortierteTiefen);


void CalcHierarchicalBasisExact(VectorSparseG &u, const Depth &Tsubgrid);



/**
 *
 * Berechnet Prewavelet-Koeffizienten, der Punkte, die in Zusammenhangskomponente enthalten sind und deren Tiefe gleich
 * komponente.Depth() sind.
 *
 * @param prew
 * @param u
 * @param komponente
 */
void GetPrewaveletCoefficients(VectorSparseG *prew, VectorSparseG *u, ZusammenhangsKomponente &komponente);


/**
 * Berechnet nodale Basis.
 *
 * @param u ist in hierarchischer Basis gegeben.
 */
void CalcNodalBasis(VectorSparseG &u);


/**
 * Berechnet nodale Basis.
 *
 * @param u ist in hierarchischer Basis gegeben.
 */
void CalcNodalBasis(VectorSparseG &u, ListOfDepthOrderedSubgrids& orderedSubgrids);


/**
 * Berechnet nodale Basis.
 *
 * @param u ist in hierarchischer Basis gegeben.
 */
void CalcNodalBasis_NEU(VectorSparseG &u, std::list<Depth>& SortierteTiefen);



/**
 * Berechnet nodale Basis auf allen Punkten mit Tiefe kleiner gleich Tiefe
 *
 * @param u ist in hierarchischer Basis gegeben.
 */
void CalcNodalBasis(VectorSparseG &u, const Depth &Tiefe, ListOfDepthOrderedSubgrids &orderedSubgrids);


/**
 * Berechnet nodale Basis auf allen Punkten mit Tiefe kleiner gleich Tiefe,
 * verwendet nur hierarchische Basis auf Punkten mit Tiefe = Depth Tiefe
 *
 * @param hier ist in hierarchischer Basis gegeben.
 */
void CalcNodalBasisExactDepth(VectorSparseG &hier, VectorSparseG &u, const Depth &Tiefe);

/**
 * Calculate hierarchical basis, given by prewavelet coefficients on a fixed subgrid
 *
 * @param prew
 * @param hier
 * @param subgrid
 */
void CalcHierbyPrew(VectorSparseG &prew, VectorSparseG &hier, SubgridFixedDepth &subgrid,
                    ListOfDepthOrderedSubgrids &orderedSubgrids);

void CalcHierbyPrew_NEU(VectorSparseG &prew, VectorSparseG &hier, Depth& T, std::list<Depth> SortierteTiefen);


void CalcHierarchicalBasis(VectorSparseG &u);

void CalcHierarchicalBasis(VectorSparseG &u,  ListOfDepthOrderedSubgrids &orderedSubgrids);

void CalcHierarchicalBasis_NEU(VectorSparseG &u, std::list<Depth>& SortierteTiefen);

void SolvePrewavelet(VectorSparseG &prew, VectorSparseG &u, const Depth &T);





void CalcUbyPrewRestrictions(VectorSparseG &prew, VectorSparseG &u, Depth Tiefe, bool const *restrictions,
                             ListOfDepthOrderedSubgrids &orderedSubgrids);



/// Prewavelet Basis

void calcPrewByNodal(VectorSparseG &prew, VectorSparseG &u);

void calcPrewByNodal_NEU(VectorSparseG &prew, VectorSparseG &u,std::list<Depth>& SortierteTiefen);

void calcPrewByNodal(VectorSparseG &prew, VectorSparseG &u, ListOfDepthOrderedSubgrids& orderedSubgrids);

void calcNodalByPrew(VectorSparseG &prew, VectorSparseG &u);


#endif //GRUN_BASISTRANSFORMATIONS_H
