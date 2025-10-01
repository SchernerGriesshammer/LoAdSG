/**********************************************************************************
 * Copyright 2010 Christoph Pflaum 
 * 		Department Informatik Lehrstuhl 10 - Systemsimulation
 *		Friedrich-Alexander Universität Erlangen-Nürnberg
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
 
// ------------------------------------------------------------
//
// vector.h
//
// ------------------------------------------------------------

#ifndef COORD_H
#define COORD_H

/*#include "../abbrevi.h"
#include "../myAssert.h"*/
#include "../indices/index.h"
#include "../primes/prime.h"
#include "../sgrid/sparseGrid.h"
#include "extempAlg.h"

//test NOW
#include <iostream>
using namespace std;

//////////////////////////////////////
// 1. coordinates
// 2. Function applied to coordinates
// 3. test Functions
//////////////////////////////////////


//////////////////////////////////////
// 1. coordinates
//////////////////////////////////////

/** \addtogroup ApplicationGroup **/
/* @{ */
class Coordinate : public ExprSparseG<Coordinate> {		
	public:			
		Coordinate(AdaptiveSparseGrid_Base& sg, int dimension_);
		
		//inline double getValue(double* data, const IndexDimension& I) const;
       inline double getValue(int i, const IndexDimension& I)const {
        
                  return I.coordinate(dimension);};   

		ExpressionDescription getDescription() const {  return ExpressionDescription(true); }		
		AdaptiveSparseGrid_Base* getSparseGrid() const { return  sparseGrid; }
	private:		
		int dimension;
		AdaptiveSparseGrid_Base* sparseGrid;		
};
//@}
/*double Coordinate::getValue(double* data, const IndexDimension& I) const {
  return I.coordinate(dimension);
}*/
/* @} */



#endif // COORD_H
