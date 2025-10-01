/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer, Rainer Hartmann
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef INTERATORBASISFUNCTION_H
#define INTERATORBASISFUNCTION_H


/***
 *  iterator of all possible basis function in a d-dimensional cell
 *  2^d possibilities
 *  Use iterator in the following form:
 * 
 *  for(IteratorBasisFunction<Dim> iterU;iterU.hasNext();iterU.next()) {
 *	 .....  iterU.getBasisTypeCoord(d);
 *	 .....  iterU.getBasisTypeNum(d);
 *  }
 * 
 */
template <size_t DIM>
class IteratorBasisFunction {
   public: 
      IteratorBasisFunction()                    { value = 0; 	maxValue=1; 	maxValue = maxValue << DIM;   }
      void next()                                { ++value; }
      bool hasNext()                             { return  value<maxValue; }
      BasisFunctionType getBasisTypeCoord(int d) { return  BasisFunctionType((value>>d)&1); }
      int               getBasisTypeNum(int d)   { return  ((value>>d)&1); }
   private:
      long value;
      long maxValue;
};





#endif  // INTERATORBASISFUNCTION_H

