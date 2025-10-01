//
// Created by scherner on 12.03.21.
//
#include "matrix.h"
#include "matrix_operations.h"

PrewaveletMatrixHomogen::PrewaveletMatrixHomogen(int size_, Depth T_, AdaptiveSparseGrid_Base *grid_) : T(T_), grid(grid_) {
    k = -1;
    size = size_;


    if (size < 6 && size >= 1) {
        Mcreated = true;
        size_R = size * size;
        M = new double[size * size];
        L = new double[size * size];
        R = new double[size * size];
        for (int i = 0; i < size * size; i++) {
            M[i] = 0.0;
            L[i] = 0.0;
            R[i] = 0.0;
        }
    } else {
        Mcreated = true;
        size_R = size * 5;
        M = new double[size * 5];
        L = new double[size * 5];
        R = new double[size * 5];
        for (int i = 0; i < size * 5; i++) {
            M[i] = 0.0;
            L[i] = 0.0;
            R[i] = 0.0;
        }
    }
    LRcreated = false;
    //Mcreated = false;

}

PrewaveletMatrixHomogen::~PrewaveletMatrixHomogen() {

    if (Mcreated) {
        delete[] M;
        delete[] L;
        delete[] R;
    }
}

int PrewaveletMatrixHomogen::convert(int i, int j, int I) {

    int savei = i;
    if (I < 6) { return i * I + j; }
    if (j > 2 && j < I - 2) {
        i = i - ((j - 3) + 1);
    }
    if (j == I - 2) i = i - j + 3;
    if (j == I - 1) i = i - j + 4;

    if (I > 5) {
        if (i > 4 || i < 0) {
            cout << "error" << endl;
            cout << " i " << i << " j  " << j << endl;
            cout << "save i  " << savei << "   I " << I << endl;
        }
    }


    return i * I + j;


}


void PrewaveletMatrixHomogen::create(IndexDimension StartIndex, int dir) {

    int depth = T.at(dir);


    unsigned long i;
    k = -1;
    while (grid->occupied(i, StartIndex) && StartIndex.getDepth(dir) != 0 && grid->workonindex(i)) {
        k++;


        //for (int k = 0; k < size; k++) {
        if (StartIndex.getDepth(dir) == depth) {
            if (StartIndex.nextLeft(dir, depth).isAtLeftBoundary(dir)) {
                M[convert(k, k, size)] = 0.9;
                if (StartIndex.getDepth(dir) == 1 && depth == 1) {
                    M[convert(k, k, size)] = 1.0;
                    break;
                }
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)) && grid->workonindex(l))
                    M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)) && grid->workonindex(l))
                    M[convert(k + 2, k, size)] = 0.1;
            } else if (StartIndex.nextRight(dir, depth).isAtRightBoundary(dir)) {
                M[convert(k, k, size)] = 0.9;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth)) && grid->workonindex(l)) {

                    M[convert(k - 1, k, size)] = -0.6;

                }
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)) && grid->workonindex(l))
                    M[convert(k - 2, k, size)] = 0.1;

            } else {
                M[convert(k, k, size)] = 1.0;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth)) && grid->workonindex(l))
                    M[convert(k - 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)) && grid->workonindex(l))
                    M[convert(k - 2, k, size)] = 0.1;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)) && grid->workonindex(l))
                    M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)) && grid->workonindex(l))
                    M[convert(k + 2, k, size)] = 0.1;

            }
        }

        if (StartIndex.getDepth(dir) < depth) {

            M[convert(k, k, size)] = 1.0;
            long unsigned int l;
            if (grid->occupied(l, StartIndex.nextLeft(dir, depth)) && grid->workonindex(l)) {
                M[convert(k - 1, k, size)] = 0.5;
            }
            if (grid->occupied(l, StartIndex.nextRight(dir, depth)) && grid->workonindex(l)) {
                M[convert(k + 1, k, size)] = 0.5;
            }


        }
        StartIndex = StartIndex.nextRight(dir, depth);


    }
};


void PrewaveletMatrixHomogen::solve(double *x, double *b) {
    if (size < 2) {
        x[0] = b[0] / M[0];
        return;
    }
    if (!LRcreated) createLR();
    double y[size];
    for (int i = 0; i < size; i++) {
        y[i] = b[i];

    }


    forward(y, b, L, size);
    backward(x, y, R, size);
};


void PrewaveletMatrixHomogen::createLR() {


    if (size >= 6) {

        for (int i = 0; i < size * 5; i++) {
            L[i] = 0.0;
            R[i] = 0.0;
        }
    } else {

        for (int i = 0; i < size * size; i++) {

            L[i] = 0.0;
            R[i] = 0.0;
        }
    }


    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (abs(i - j) < 3) {
                R[convert(i, j, size)] = M[convert(i, j, size)];

            }

        }

    }

    LR(L, R, size);
};

void PrewaveletMatrixHomogen::print() {
    int I = size;
    for (int i = 0; i < I; i++) {
        for (int j = 0; j < I; j++) {
            if (i <= j + 2 && j <= i + 2) {
                cout << "\t" << M[convert(i, j, I)];
            } else { cout << "\t"; }
        }
        cout << "\n";
    }
    cout << endl;

};


/*void PrewaveletMatrix_Neumann::create_neumann(IndexDimension StartIndex, int dir) {

    int depth = T.at(dir);
    unsigned long i;
    k = -1;
    while (grid->occupied(i, StartIndex)) {
        k++;

        //for (int k = 0; k < size; k++) {
        if (StartIndex.getDepth(dir) == depth) {
            if (StartIndex.nextLeft(dir, depth).isAtLeftBoundary(dir)) {
                M[convert(k - 1, k, size)] = -1.2;
                M[convert(k, k, size)] = 1.1;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)))
                    M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)))
                    M[convert(k + 2, k, size)] = 0.1;
            } else if (StartIndex.nextRight(dir, depth).isAtRightBoundary(dir)) {
                M[convert(k, k, size)] = 1.1;
                M[convert(k + 1, k, size)] = -1.2;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth))) {

                    M[convert(k - 1, k, size)] = -0.6;

                }
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)))
                    M[convert(k - 2, k, size)] = 0.1;

            } else {
                M[convert(k, k, size)] = 1.0;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth)))
                    M[convert(k - 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)))
                    M[convert(k - 2, k, size)] = 0.1;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)))
                    M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)))
                    M[convert(k + 2, k, size)] = 0.1;

            }
        }

        if (StartIndex.getDepth(dir) < depth) {

            M[convert(k, k, size)] = 1.0;
            long unsigned int l;
            if (grid->occupied(l, StartIndex.nextLeft(dir, depth)) && !StartIndex.isAtLeftBoundary(dir)) {
                M[convert(k - 1, k, size)] = 0.5;
            }
            if (grid->occupied(l, StartIndex.nextRight(dir, depth)) && !StartIndex.isAtRightBoundary(dir)) {
                M[convert(k + 1, k, size)] = 0.5;
            }


        }
        if (!StartIndex.isAtRightBoundary(dir)) StartIndex = StartIndex.nextRight(dir, depth);
        else break;


    }
}*/


//////////////////////////////////////////////////

PrewaveletMatrixInhomogen::PrewaveletMatrixInhomogen(int size_, Depth T_, AdaptiveSparseGrid_Base *grid_) : T(T_),
                                                                                                            grid(grid_) {
    k = -1;
    size = size_;

    if (size < 6 && size >= 1) {
        Mcreated = true;
        size_R = size * size;
        M = new double[size * size];
        L = new double[size * size];
        R = new double[size * size];
        for (int i = 0; i < size * size; i++) {
            M[i] = 0.0;
            L[i] = 0.0;
            R[i] = 0.0;
        }
    } else {
        Mcreated = true;
        size_R = size * 5;
        M = new double[size * 5];
        L = new double[size * 5];
        R = new double[size * 5];
        for (int i = 0; i < size * 5; i++) {
            M[i] = 0.0;
            L[i] = 0.0;
            R[i] = 0.0;
        }
    }
    LRcreated = false;
    //Mcreated = false;

}

PrewaveletMatrixInhomogen::~PrewaveletMatrixInhomogen() {

    if (Mcreated) {

        delete[] M;
        delete[] L;
        delete[] R;
    }
}

int PrewaveletMatrixInhomogen::convert(int i, int j, int I) {

    int savei = i;
    int savej = j;
    if (I < 6) { return i * I + j; }
    if (j > 2 && j < I - 2) {

        i = i - ((j - 3) + 1);
    }
    if (j == I - 2) {


        i = i - j + 3;

    }
    if (j == I - 1) {

        i = i - j + 4;
    }

    if (I > 5) {
        if (i > 4 || i < 0) {
            cout << "error" << endl;
            cout << " i " << i << " j  " << j << endl;
            cout << "save i  " << savei << endl;
            cout << "save j  " << savej << endl;
            cout << "I " << I << endl;

            throw std::invalid_argument("received negative value");
        }
    }


    return i * I + j;


}


void PrewaveletMatrixInhomogen::create(IndexDimension StartIndex, int dir) {


    int depth = T.at(dir);
    unsigned long i;
    k = -1;

    while (grid->occupied(i, StartIndex)) {
        if (!grid->workonindex(i)) break;

        k++;


        //for (int k = 0; k < size; k++) {
        if (StartIndex.getDepth(dir) == depth) {


            if (StartIndex.getDepth(dir) == 1 && depth == 1) {
                M[convert(0, 0, size)] = 1;
                M[convert(0, 1, size)] = -1;
                M[convert(1, 1, size)] = 1;
                M[convert(1, 0, size)] = 0.5;
                M[convert(1, 2, size)] = 0.5;
                M[convert(2, 1, size)] = -1;
                M[convert(2, 2, size)] = 1;
                break;
            }

            if (StartIndex.getDepth(dir) == 0 && depth == 0) {
                M[convert(k, k, size)] = 1.0;
                break;
            }
            if (StartIndex.nextLeft(dir, depth).isAtLeftBoundary(dir)) {

                if (k == 0)k++;
                M[convert(k - 1, k, size)] = -1.2;

                M[convert(k, k, size)] = 1.1;

                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)))
                    M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k + 2, k, size)] = 0.1;
            } else if (StartIndex.nextRight(dir, depth).isAtRightBoundary(dir)) {
                M[convert(k, k, size)] = 1.1;
                if (StartIndex.getDepth(dir) <= 1 && depth <= 1) {
                    M[convert(k, k, size)] = 1.0;
                    break;
                }


                M[convert(k + 1, k, size)] = -1.2;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth))) {
                    if (grid->workonindex(l))
                        M[convert(k - 1, k, size)] = -0.6;

                }
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k - 2, k, size)] = 0.1;

            } else {
                M[convert(k, k, size)] = 1.0;
                long unsigned int l;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k - 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextLeft(dir, depth).nextLeft(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k - 2, k, size)] = 0.1;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k + 1, k, size)] = -0.6;
                if (grid->occupied(l, StartIndex.nextRight(dir, depth).nextRight(dir, depth)))
                    if (grid->workonindex(l))
                        M[convert(k + 2, k, size)] = 0.1;

            }
        }

        if (StartIndex.getDepth(dir) < depth) {

            M[convert(k, k, size)] = 1.0;
            long unsigned int l;
            if (grid->occupied(l, StartIndex.nextLeft(dir, depth)) && !StartIndex.isAtLeftBoundary(dir)) {
                if (grid->workonindex(l))
                    M[convert(k - 1, k, size)] = 0.5;
            }
            if (grid->occupied(l, StartIndex.nextRight(dir, depth)) && !StartIndex.isAtRightBoundary(dir)) {
                if (grid->workonindex(l))
                    M[convert(k + 1, k, size)] = 0.5;
            }


        }
        if (!StartIndex.isAtRightBoundary(dir)) StartIndex = StartIndex.nextRight(dir, depth);
        else break;


    }
};


void PrewaveletMatrixInhomogen::solve(double *x, double *b) {
    if (size < 2) {
        x[0] = b[0] / M[0];
        return;
    }
    if (!LRcreated) createLR();
    double y[size];
    for (int i = 0; i < size; i++) {
        y[i] = b[i];

    }


    forward(y, b, L, size);
    backward(x, y, R, size);
};


void PrewaveletMatrixInhomogen::createLR() {


    if (size >= 6) {

        for (int i = 0; i < size * 5; i++) {
            L[i] = 0.0;
            R[i] = 0.0;
        }
    } else {

        for (int i = 0; i < size * size; i++) {

            L[i] = 0.0;
            R[i] = 0.0;
        }
    }


    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (abs(i - j) < 3) {
                R[convert(i, j, size)] = M[convert(i, j, size)];

            }

        }

    }

    LR(L, R, size);
};

void PrewaveletMatrixInhomogen::print() {
    int I = size;
    for (int i = 0; i < I; i++) {
        for (int j = 0; j < I; j++) {
            if (i <= j + 2 && j <= i + 2) {
                cout << "\t" << M[convert(i, j, I)];
            } else { cout << "\t"; }
        }
        cout << "\n";
    }
    cout << endl;

};

void PrewaveletMatrixInhomogen::printL() {
    int I = size;
    for (int i = 0; i < I; i++) {
        for (int j = 0; j < I; j++) {
            if (i <= j + 2 && j <= i + 2) {
                cout << "\t" << L[convert(i, j, I)];
            } else { cout << "\t"; }
        }
        cout << "\n";
    }
    cout << endl;

};

