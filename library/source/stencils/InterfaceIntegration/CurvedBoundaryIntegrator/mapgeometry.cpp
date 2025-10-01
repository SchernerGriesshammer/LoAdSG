
#include <iostream>
#include <fstream>

//test NOW
using namespace std;

#include "../interfaceMatrices.h"
#include "../constantIntegrators.h"
#include "../interatorBasisFunction.h"
#include "curvedPoissonIntegrator.h"
#include "SGIntegrator.h"

#define Dim 3

double solution(D3vector& x_){
    double x = x_[0],y = x_[1], z = x_[2];

    return(x*y*z*(z-0.5)*cos(0.5*M_PI*(x*x+y*y))*sin(0.5*M_PI*(x*x+y*y)));
}
//Transformation functions and their derivatives
D3vector transform_SD(double t){
    return D3vector(1.0+t,0.0,0.0);
}

D3vector transDerivs_SD(double t){
    return D3vector(1.0,0.0,0.0);
}

D3vector transform_ND(double t){
    return D3vector(0.0,1.0+t,0.0);
}

D3vector transDerivs_ND(double t){
    return D3vector(0.0,1.0,0.0);
}

D3vector transform_ST(double t){
    return D3vector(1.0+t,0.0,0.5);
}

D3vector transDerivs_ST(double t){
    return D3vector(1.0,0.0,0.0);
}

D3vector transform_NT(double t){
    return D3vector(0.0,1.0+t,0.5);
}

D3vector transDerivs_NT(double t){
    return D3vector(0.0,1.0,0.0);
}


D3vector transform_WD(double t){
    return D3vector(cos(t*0.5*M_PI), sin(t*0.5*M_PI), 0.0 );
}

D3vector transDerivs_WD(double t){
    return D3vector(-0.5*M_PI*sin(t*0.5*M_PI),0.5*M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_ED(double t){
    return D3vector(2.0*cos(t*0.5*M_PI), 2.0*sin(t*0.5*M_PI), 0.0 );
}

D3vector transDerivs_ED(double t){
    return D3vector(-1.0*M_PI*sin(t*0.5*M_PI),M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_WT(double t){
    return D3vector(cos(t*0.5*M_PI), sin(t*0.5*M_PI), 0.5 );
}

D3vector transDerivs_WT(double t){
    return D3vector(-0.5*M_PI*sin(t*0.5*M_PI),0.5*M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_ET(double t){
    return D3vector(2.0*cos(t*0.5*M_PI), 2.0*sin(t*0.5*M_PI), 0.5);
}

D3vector transDerivs_ET(double t){
    return D3vector(-1.0*M_PI*sin(t*0.5*M_PI),M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_SW(double t){
    return D3vector(1.0,0.0,0.5*t);
}

D3vector transDerivs_SW(double t){
    return D3vector(0.0,0.0,0.5);
}


D3vector transform_SE(double t){
    return D3vector(2.0,0.0,0.5*t);
}

D3vector transDerivs_SE(double t){
    return D3vector(0.0,0.0,0.5);
}


D3vector transform_NW(double t){
    return D3vector(0.0,1.0,0.5*t);
}

D3vector transDerivs_NW(double t){
    return D3vector(0.0,0.0,0.5);
}

D3vector transform_NE(double t){
    return D3vector(0.0,2.0,0.5*t);
}

D3vector transDerivs_NE(double t){
    return D3vector(0.0,0.0,0.5);
}

int main(int argc, char *argv[]){
    cout.precision(5);
    cout.setf(std::ios::fixed,std::ios::floatfield);


    // Example: POissin Curved
    //////////////////////////////////
    D3vector corners[8];
    D3vector (*transformEdge[12])(double);
    D3vector (*transDerivs[12])(double);

    // this describes a small ring segment of radii r, R and length L:
    double r = 1.0;
    double R = 2.0;
    double L = 0.5;

    corners[WSDdir3D] = D3vector ( r,0.0,0.0 );
    corners[ESDdir3D] = D3vector ( R,0.0,0.0 );
    corners[WNDdir3D] = D3vector ( 0.0,r,0.0 );
    corners[ENDdir3D] = D3vector ( 0.0,R,0.0 );
    corners[WSTdir3D] = D3vector ( r,0.0,L );
    corners[ESTdir3D] = D3vector ( R,0.0,L );
    corners[WNTdir3D] = D3vector ( 0.0,r,L );
    corners[ENTdir3D] = D3vector ( 0.0,R,L );

    /*transformEdge[SDed] = transform_const;
    transformEdge[NDed] = transform_const;
    transformEdge[STed] = transform_const;
    transformEdge[NTed] = transform_const;
    transformEdge[WDed] = transform_cyl_NE_inner;
    transformEdge[EDed] = transform_cyl_NE;
    transformEdge[WTed] = transform_cyl_NE_inner;
    transformEdge[ETed] = transform_cyl_NE;
    transformEdge[SWed] = transform_const;
    transformEdge[SEed] = transform_const;
    transformEdge[NWed] = transform_const;
    transformEdge[NEed] = transform_const;*/

    transformEdge[SDed] = transform_SD;
    transformEdge[NDed] = transform_ND;
    transformEdge[STed] = transform_ST;
    transformEdge[NTed] = transform_NT;
    transformEdge[WDed] = transform_WD;
    transformEdge[EDed] = transform_ED;
    transformEdge[WTed] = transform_WT;
    transformEdge[ETed] = transform_ET;
    transformEdge[SWed] = transform_SW;
    transformEdge[SEed] = transform_SE;
    transformEdge[NWed] = transform_NW;
    transformEdge[NEed] = transform_NE;

    transDerivs[SDed] = transDerivs_SD;
    transDerivs[NDed] = transDerivs_ND;
    transDerivs[STed] = transDerivs_ST;
    transDerivs[NTed] = transDerivs_NT;
    transDerivs[WDed] = transDerivs_WD;
    transDerivs[EDed] = transDerivs_ED;
    transDerivs[WTed] = transDerivs_WT;
    transDerivs[ETed] = transDerivs_ET;
    transDerivs[SWed] = transDerivs_SW;
    transDerivs[SEed] = transDerivs_SE;
    transDerivs[NWed] = transDerivs_NW;
    transDerivs[NEed] = transDerivs_NE;


    IntegratorPoissonCurved ipc(transformEdge,transDerivs,corners);

    //Number of elements = N*N*N
    int N = atoi(argv[1]);

    double dh = 1.0/(double)N;

    ofstream outfile;
    outfile.open("geometry.txt");

    D3vector x(0.0,0.0,0);
    D3vector x_trans(0,0,0);

    for(int k = 0;k<=N;k++){
        for(int j=0;j<=N;j++){
            for(int i=0;i<=N;i++){
             x[0] = (i)*dh;
             x[1] = (j)*dh;
             x[2] = (k)*dh;

             x_trans = ipc.Map2(x[0],x[1],x[2]);

             outfile<<x_trans[0]<<" "<<x_trans[1]<<" "<<x_trans[2]<<std::endl;
            }

        }
    }

}
