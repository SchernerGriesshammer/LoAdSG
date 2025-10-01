  
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

#ifndef MYSPARSEMATH_H
#define MYSPARSEMATH_H

inline double ABS(double x) {
    if (x >= 0.0) return x;
    return -x;
}

inline int POW(int base, unsigned int exponent) {
    int returnvalue = 1;
    for (int j = 0; j < exponent; j++) {
        returnvalue = returnvalue * base;
    }
    return returnvalue;
}

inline int POW2(unsigned int exponent) {
    unsigned int one=1;
    return one << exponent;
}


/**
 * Converts an integer to a binary and stores the resulting binary as a boolean.
 *
 * @param num
 * @param binary
 */
inline void intToBinary(int num, bool binary[], int dim) {
    for (int i = dim - 1; i >= 0; --i) {
        binary[i] = num & 1;
        num >>= 1;
    }
}
#endif
