/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
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
