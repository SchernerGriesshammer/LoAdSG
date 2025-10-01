// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>

//test NOW
using namespace std;

#include "sour"
#include "../constantIntegrators.h"
#include "../interatorBasisFunction.h"
#include "helmIntegrator.h"

#define Dim 3

double integral_analytical( double p_left[Dim], double p_right[Dim], BasisFunctionType u[Dim], BasisFunctionType v[Dim] )
{

    double integral = 1.0;

    for( size_t i = 0; i < Dim; ++i )
    {

        if( u[i] == v[i] )
        {
            integral *= ( (p_right[i] - p_left[i]) / 3.0 );
        }
        else
        {
            integral *= ( (p_right[i] - p_left[i]) / 6.0 );
        }

    }

    return integral;

}

static double var_coeff( const std::array<double, Dim>& x )
{
    return 1.0;
}

int main(int argc, char** argv)
{
    cout.precision(5);
    cout.setf(std::ios::fixed,std::ios::floatfield);

    IntegratorHelm<double (*) (const std::array<double, Dim>&), Dim> localStiffnessHelm(&var_coeff, 1e-8, 0, 40, 1.0e8 );

    double p_left[Dim];
    double p_right[Dim];

    p_left[0] = -0.1;
    p_right[0] = 0.3;

    for(int i=1; i<Dim; ++i)
    {
        p_left[i]  = 2.0 * p_left[i-1];
        p_right[i] = 3.0 * p_right[i-1];
    }

    BasisFunctionType u[Dim];
    BasisFunctionType v[Dim];

    for(IteratorBasisFunction<Dim> iterU; iterU.hasNext(); iterU.next())
    {
        cout << endl;
        cout << " u = (";
        for(int d=0; d<Dim; ++d)
        {
            cout << iterU.getBasisTypeNum(d);
            if(d<Dim-1)
            {
                cout << ", ";
            }
        }
        cout << ")" << endl;
        cout << " ----------- " << endl;
        for(IteratorBasisFunction<Dim> iterV; iterV.hasNext(); iterV.next())
        {
            cout << " v = (";
            for(int d=0; d<Dim; ++d)
            {
                cout << iterV.getBasisTypeNum(d);
                if(d<Dim-1)
                {
                    cout << ", ";
                }
            }

            for(int d=0; d<Dim; ++d)
            {
                u[d] =iterU.getBasisTypeCoord(d);
                v[d] =iterV.getBasisTypeCoord(d);
            }

            cout << "), Error: "
                 << std::abs(localStiffnessHelm.stencil_integration(p_left,p_right,u,v) - integral_analytical(p_left,p_right,u,v))
                 << endl;
        }
    }













    cout << " Ende! " << endl;
}


