/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer, Rainer Hartmann
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef INTERFACEMATRICES_H
#define INTERFACEMATRICES_H

/**
 * basic  enum for two possible basis functions on a 1-dimensional interval 
 */
enum BasisFunctionType { leftBasis = 0, rightBasis = 1, gradLeftBasis = 2, gradRightBasis = 3 };


/**
 * interface for calculation of local stiffness matrices
 */
template <size_t DIM>
class InterfaceLocalStiffnessMatrices {
   public:
      /**
       * @param p_left  coordinate of left  boundary of cell
       * @param p_right coordinate of right boundary of cell
       * @param u,v:    vector describing type of basis function 
       * @return a(u,v) with respect to the integration on the cell
       */         
      virtual double stencil_integration(double p_left[DIM], double p_right[DIM],
				                                 BasisFunctionType u[DIM], BasisFunctionType v[DIM]) const = 0;
};


#endif  // INTERFACEMATRICES_H

