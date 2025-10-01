//
// Created by scherner on 01.03.21.
//

#ifndef GRUN_MATRIX_H
#define GRUN_MATRIX_H

//#include "prewavelet_matrix.h"
#include "../sgrid/depth.h"
#include "../sgrid/sparseGrid.h"

class PrewaveletMatrixHomogen {
public:
    PrewaveletMatrixHomogen(int size_, Depth T_, AdaptiveSparseGrid_Base *grid_);

    ~PrewaveletMatrixHomogen();

    static int convert(int i, int j, int I);

    void create(IndexDimension StartIndex, int dir);

    void solve(double *x, double *b);


    void createLR();

    void print();


protected:
    double *M;
    double *R;
    double *L;
    bool LRcreated;
    bool Mcreated;
    int size;
    int size_R;
    AdaptiveSparseGrid_Base *grid;

    Depth T;
    int k;
};


class PrewaveletMatrixInhomogen {
public:
    PrewaveletMatrixInhomogen(int size_, Depth T_, AdaptiveSparseGrid_Base *grid_);

    ~PrewaveletMatrixInhomogen();

    static int convert(int i, int j, int I);

    void create(IndexDimension StartIndex, int dir);

    void solve(double *x, double *b);


    void createLR();

    void print();

    void printL();


protected:
    double *M;
    double *R;
    double *L;
    bool LRcreated;
    bool Mcreated;
    int size;
    int size_R;
    AdaptiveSparseGrid_Base *grid;

    Depth T;
    int k;
};


#endif //GRUN_MATRIX_H
