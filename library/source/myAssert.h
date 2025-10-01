/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
#ifndef MYASSERT_H
#define MYASSERT_H


#include <assert.h>
#include "abbrevi.h"





inline void myAssert(bool thisShouldhold) {
#ifdef ASSERT
   assert(thisShouldhold);
#endif
}

inline void assertDimension(const int d) {
#ifdef ASSERT
  assert(d>=0);
  assert(d<DimensionSparseGrid);
#endif
}



#endif  // MYASSERT_H

