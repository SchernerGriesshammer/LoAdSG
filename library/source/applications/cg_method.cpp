//
// Created by scherner on 26.07.21.
//

#include "cg_method.h"


#include "../stencils/PoissonStencil.h"
#include "../stencils/Stencil.h"

bool
CG::solveNeumann(int maxIteration, double eps, VectorSparseG &x,  VectorSparseG &f,
                 int *iterations, MatrixVectorInhomogen& matrix) {

    cout << "Start CG " << endl;
    double tau;
    double delta;
    double delta_prime;
    double beta;



    AdaptiveSparseGrid *grid = x.getSparseGrid();

    VectorSparseG z(*grid),r(*grid),g(*grid),d(*grid);

    VectorSparseG vec(*grid);
    VectorSparseG ddelta(*grid);
    VectorSparseG xt(*grid);
    VectorSparseG h(*grid);

    double j =1.0;





    // solves  Ax = f;




    // r = A*x - f;


    //MatrixVectorPrewavelet_inhomog(x, vec, Stiffness,gM);

    PoissonStencil poissonStencil;
    matrix.multiplication<PoissonStencil>(x, vec, poissonStencil);


    //matrixH.multiplication<StiffnessStencilH>(x,vec,stencilH);
    x = 0.0;
    OperatorT::apply(vec, r);
    vec = 0.0;

    // ftilde = T*f;
    OperatorT::apply(f, vec);

    // r = r - f = Ax - f;
    r = vec - r;
    h = j * r;
    ddelta = h;


    delta = product(r, h);


    cout << "Solve cg" << endl;
    int k = 0;

    for (int i = 1; i <= maxIteration && delta > eps; ++i) {
        grid->WorkOnHangingNodes = false;
        k = i;

        // T*A*T*d=T*g;

/*
        OperatorTtrans(ddelta,d);
        CalcUbyPrew_Neumann(d,test);
        for(unsigned long k=0; k< grid->getMaximalOccupiedSecondTable(); k++){
            IndexDimension I = grid->getIndexOfTable(k);
            if(I.isAtBoundary() && abs(test.getValue(I))>0.01){
                cout << "error in Iteration " << i << endl;
                I.PrintCoord();
                cout << endl;
                cout << abs(test.getValue(I)) << endl;


           }
        }
*/

        g = 0.0;

        OperatorT::applyTranspose(ddelta, d);

        matrix.multiplication<PoissonStencil>(d,vec,poissonStencil);

        OperatorT::apply(vec, g);


        tau = double(delta / product(ddelta, g));

        r = r - (tau * g);
        x = x + (tau * ddelta);


        delta_prime = product(r, r);
        if (delta_prime < eps) break;
        h = j * r;

        delta_prime = product(r, h);


        beta = delta_prime / delta;
        delta = delta_prime;
        ddelta = h + beta * ddelta;


    }

    cout << "CG finished after " << k << " iterations and residuum delta = " << delta_prime << endl;
    *iterations = k;
    vec = 0.0;
    OperatorT::applyTranspose(x, vec);
    x = 0.0;
    x = vec;
    return true;
}



double CG::error[MaximumDepth][maxIteration] = {};

int CG::count = 0;


void OperatorT::apply(VectorSparseG &input, VectorSparseG &output) {
    AdaptiveSparseGrid_Base *grid = input.getSparseGrid();
    VectorSparseG input_cp(*input.getSparseGrid());
    input_cp = input;
    output = input;


    for (int d = 0; d < DimensionSparseGrid; d++) {
        for (unsigned long i = 0; i < grid->getMaximalOccupiedSecondTable(); i++) {
            if (grid->getSecondTable()[i] != 0) {
                IndexDimension Index = grid->getIndexOfTable(i);
                // if Index is not at Boundary
                if (Index.getIndex(d) != 0 && Index.getIndex(d) != 1) {

                    if (Index.isLinksRandNah(d) && Index.isRechtsRandNah(d)) {

                        double value = input.getValue(i) + input.getValue(Index.nextLeft(d)) +
                                       input.getValue(Index.nextRight(d));

                        output.setValue(i,value);

                    }

                    if (Index.isLinksRandNah(d) && !Index.isRechtsRandNah(d)) {
/*                        ShiftOperator left(d, Left);
                        output = input + 1.2 * left(input) | Index;*/
                        double value = input.getValue(i) + 1.2 * input.getValue(Index.nextLeft(d));
                        //output = value | Index;
                        output.setValue(i,value);
                    }

                    if (Index.isRechtsRandNah(d) && !Index.isLinksRandNah(d)) {
                        /*                ShiftOperator right(d, Right);
                                        output = input + 1.2 * right(input) | Index;*/

                        double value = input.getValue(i) + 1.2 * input.getValue(Index.nextRight(d));
                        //output = value | Index;
                        output.setValue(i,value);
                    }


                } else {
                    //output = 0.0 | Index;
                    output.setValue(i,0.0);
                }

            }

        }
        input = output;


    }
    input = input_cp;
}

void OperatorT::applyTranspose(VectorSparseG &input, VectorSparseG &output) {
    VectorSparseG input_cp(*input.getSparseGrid());
    input_cp = input;
    output = 0.0;
    AdaptiveSparseGrid_Base *grid = input.getSparseGrid();
    for (int d = 0; d < DimensionSparseGrid; d++) {
        for (unsigned long i = 0; i < grid->getMaximalOccupiedSecondTable(); i++) {
            if (grid->getSecondTable()[i] != 0) {
                IndexDimension Index = grid->getIndexOfTable(i);
                if (Index.getIndex(d) != 0 && Index.getIndex(d) != 1) {

                    output = input.getValue(Index) | Index;
                    if (Index.isLinksRandNah(d) && Index.isRechtsRandNah(d)) {
                        IndexDimension IndexLeft = Index.nextLeft(d);
                        IndexDimension IndexRight = Index.nextRight(d);

                        double value = input.getValue(Index);
                        //output = output + input + value | IndexLeft;
                        output = output.getValue(IndexLeft) + input.getValue(IndexLeft) + input.getValue(Index) |
                                 IndexLeft;
                        //output = output + input + value | IndexRight;
                        output = output.getValue(IndexRight) + input.getValue(IndexRight) + value | IndexRight;
                    }
                    if (Index.isRechtsRandNah(d) && !Index.isLinksRandNah(d)) {

                        IndexDimension IndexRight = Index.nextRight(d);
                        double value = input.getValue(Index);
                        //output = output + input + 1.2 * value | IndexRight;
                        output = output.getValue(IndexRight) + input.getValue(IndexRight) + 1.2 * value | IndexRight;
                    }
                    if (!Index.isRechtsRandNah(d) && Index.isLinksRandNah(d)) {

                        IndexDimension IndexLeft = Index.nextLeft(d);
                        double value = input.getValue(Index);
                        //output = output + input + 1.2 * value | IndexLeft;
                        output = output.getValue(IndexLeft) + input.getValue(IndexLeft) + 1.2 * value | IndexLeft;
                    }


                }

            }
        }

        input = output;
    }
    input = input_cp;

}



