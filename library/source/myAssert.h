  
/**********************************************************************************
 * Copyright 2015 Christoph Pflaum 
 *              Department Informatik Lehrstuhl 10 - Systemsimulation
 *              Friedrich-Alexander Universität Erlangen-Nürnberg
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

