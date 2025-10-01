
  
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