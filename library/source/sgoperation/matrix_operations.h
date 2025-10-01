

#ifndef PREW_MATRIX_H
#define PREW_MATRIX_H



#include "../abbrevi.h"
#include "../myAssert.h"
#include "../indices/index.h"
#include "../primes/prime.h"

#include "../sgrid/sparseGrid.h"


/*#include "../extemp/vector.h"
#include "../extemp/shift.h"
#include "../extemp/meshExtemp.h"*/



#include <iostream>
#include <cmath>

using namespace std;


void LR(double *L, double *R, int I);

void forward(double *y, double *lambda, double *L, int I);

void backward(double *c, double *y, double *R, int I);


#endif
