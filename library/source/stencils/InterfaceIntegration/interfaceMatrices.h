/**********************************************************************************
 * Copyright 2016 Christoph Pflaum, Rainer Hartmann
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **********************************************************************************/

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

