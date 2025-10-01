/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

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
