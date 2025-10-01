/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef PRIME_H
#define PRIME_H

class PrimeNumbers {
public:
    static unsigned int getPrimeForHash(int d) { return numberDim[d]; }

    static unsigned long getNextPrime(unsigned long n) { return 199831; } // { return 11;  }
private:
    static const unsigned int numberDim[20];
};

#endif // PRIME_H