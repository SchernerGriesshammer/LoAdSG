//
// Created by to35jepo on 9/19/23.
//

#include "condition.h"
#include "../extemp/vector.h"
#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "../stencils/Stencil.h"
#include "cg_method.h"
#include "norms.h"
#include "H1_Norm.h"


double potenzmethode(int k,AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, bool conditioned){

    VectorSparseG h(grid),xneu(grid),z(grid), xold(grid);


    MatrixVectorHomogen matrix(grid,mgrid,1,0);
    double rho;
    double rho_Z;
    double rho_N;
    //ADD<Poisson, VariableCoefficientMass> stencil(&grid);
    Poisson stencil(grid);
    Preconditioning<Poisson> P(z,stencil);
    IndexDimension centerPoint;
    xold=1.0;
    double rho_old=0.0;
    for(int i = 0;  i<k; i++) {




        h = xold;



        matrix.multiplication(h,xneu,stencil);





        rho_Z = product(xold, xneu);
        rho_N= product(xneu,xneu);
        rho =rho_Z;

        if(abs(rho-rho_old)<1e-15)break;
        rho_old=rho;

        double sum = 0.0;

        sum = 1.0 /sqrt(rho_N);
        xold = sum * xneu;

    }



    return rho;
}

double potenzmethode_min(int k,AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid,bool conditioned, double lambda){

    VectorSparseG h(grid),xneu(grid),z(grid), xold(grid),l(grid),l2(grid),xneu_v(grid);


    MatrixVectorHomogen matrix(grid,mgrid,1,0);
    double rho;
    double rho_Z;
    double rho_N;
    //ADD<Poisson, VariableCoefficientMass> stencil(&grid);
    Poisson stencil(grid);
    HelmHoltz helmHoltz(grid);
    Preconditioning<Poisson> P(z,stencil);

    IndexDimension centerPoint;
    xold.setValue(centerPoint,1.0);

    double rho_old=0.0;
    for(int i = 0;  i<k; i++) {




        h=xold;
        l=xold;


        xneu =0.0;

        matrix.multiplication(h,xneu,stencil);
        //P.apply_inverse(&xneu);




        l2 = lambda*l;
        xneu = l2 - xneu;






        rho_Z = product(xold, xneu);
        rho_N= product(xneu,xneu);
        rho =rho_Z;
        if(abs(rho_old-rho)<1e-15)break;
        rho_old=rho;


        double sum;
        sum = rho_N;

        sum = 1.0 /sqrt(sum);
        xold = sum * xneu;

    }

    xneu.PrintDouble(3);
    xold.PrintDouble(3);


    rho = -rho+lambda;


    return rho;
}


double potenzmethode_invers(int k, AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, VectorSparseG& emin){
    //Finde Eigenwert zu Ax=lambda*Mx <=> inv(M)Ax = lambda x

    VectorSparseG x(grid),xneu(grid),diff(grid),xv(grid);
    MatrixVectorHomogen matrix(grid,mgrid,1,0);

    // rho = Raiyleigh-Quotient
    double rho=0.0;
    double rho_old=1.0;


    // Setze Startwert für x
    x=1.0;

    // Input für das CG-Verfahren
    int iterations = 100;
    double time;

    // Definiere Stencil für die Matrizen
    Poisson A(grid);
    HelmHoltz M(grid);

    int i;
    for(i = 0;  i<k; i++) {

        xneu=0.0;

        // Berechne M*x
        matrix.multiplication(x,xv,M);

        // Löse Ax_neu = Mx
        CG::solveHomogen<Poisson>(1e-15, xneu,xv, &iterations, diff, matrix, A, &time);

        // Berechne Rayleigh-Quotienten
        rho = product(x,xneu)/ product(x,x);


        // Setze x = xneu/c;
        L2 l2(grid,mgrid);

        double c = product(xneu,xneu);
        x = xneu/sqrt(c);



        //Abbruch
        if(abs(rho_old-rho)<1e-15) break;
        rho_old=rho;

    }


    emin = x;

    cout << "Abbruch nach " << i << " Iterationen " << endl;
    cout << "Lambda min =  " << 1.0/rho << endl;
    cout << " abs(rho_old-rho) = " << abs(rho_old-rho) << endl;
    return 1.0/rho;

}




double condition(int k, AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid&mgrid, bool conditioned){
    double lmax = potenzmethode_min(k,grid,mgrid,conditioned,0.0);
    double lmin= potenzmethode_min(k,grid,mgrid,conditioned,lmax);
    //double lmin2  = potenzmethode_invers(k,grid,mgrid);
    cout <<"lmax: "  << lmax << " lmin: " << lmin <<  endl;
    return lmax/lmin;
}
