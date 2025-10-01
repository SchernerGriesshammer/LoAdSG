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

