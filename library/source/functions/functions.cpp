//
// Created by to35jepo on 8/1/24.
//

#include "functions.h"
#include "../indices/index.h"
#include <cmath>
#include <iostream>

double schroedinger2D(double a, double b, double x, double y) {
    // Compute the exponential term
    double exp_term = exp(a * (pow(x - b, 2) + pow(y - b, 2)));

    // Compute the dot product of vecx and vecy
    double dot_product = (x - y) * (x - y);

    // Compute the polynomial terms
    double x_term = pow(x, 3) * pow(x - 1.0, 3);
    double y_term = pow(y, 3) * pow(y - 1.0, 3);

    // Compute the final function value
    double f = exp_term * dot_product * x_term * y_term;

    return f;
}

double lap_schroedinger2D(double a, double b, double x, double y) {
    double exp_part = std::exp(a * (std::pow(x - 0.5, 2) + std::pow(y - 0.5, 2)));
    double x1 = x - 1.0;
    double y1 = y - 1.0;
    double x0_5 = x - 0.5;
    double y0_5 = y - 0.5;
    double common_term = std::pow(y1, 3) * std::pow(y, 3) * std::pow(x - y, 2);

    double term1 = x * x * x * ((std::pow(x1, 3) * (4.0 * a * a * x0_5 * x0_5 * exp_part + 2.0 * a * exp_part)) + 12.0 * a * x0_5 * x1 * x1 * exp_part + 6.0 * x1 * exp_part);
    double term2 = 6.0 * x * x * (2.0 * a * x0_5 * std::pow(x1, 3) * exp_part + 3.0 * x1 * x1 * exp_part);
    double term3 = 6.0 * x * std::pow(x1, 3) * exp_part;
    double term4 = 2.0 * std::pow(x1, 3) * x * x * x * std::pow(y1, 3) * y * y * y * exp_part;
    double term5 = 4.0 * std::pow(y1, 3) * std::pow(y, 3) * (x - y) * (2.0 * a * x0_5 * std::pow(x, 3) * std::pow(x1, 3) * exp_part + 3.0 * std::pow(x, 3) * x1 * x1 * exp_part + 3.0 * x * x * std::pow(x1, 3) * exp_part);

   double d2dx =  common_term * (term1 + term2 + term3) + term4 + term5;

    double common_term1 = std::pow(x1, 3) * std::pow(x, 3) * std::pow(y1, 3) * std::pow(y, 3) * std::pow(x - y, 2);
    term1 = common_term1 * (4.0 * a * a * y0_5 * y0_5 * exp_part + 2.0 * a * exp_part);

    double common_term2 = std::pow(x1, 3) * std::pow(x, 3) * exp_part;
    term2 = common_term2 * (
            std::pow(y, 3) * (-12.0 * std::pow(y1, 2) * (x - y) + 6.0 * y1 * std::pow(x - y, 2) + 2.0 * std::pow(y1, 3)) +
            6.0 * std::pow(y, 2) * (3.0 * std::pow(y1, 2) * std::pow(x - y, 2) - 2.0 * std::pow(y1, 3) * (x - y)) +
            6.0 * y * std::pow(y1, 3) * std::pow(x - y, 2)
    );

    double common_term3 = 4.0 * a * std::pow(x1, 3) * std::pow(x, 3) * y0_5 * exp_part;
    term3 = common_term3 * (
            -2.0 * std::pow(y, 3) * std::pow(y1, 3) * (x - y) +
            3.0 * std::pow(y, 3) * std::pow(y1, 2) * std::pow(x - y, 2) +
            3.0 * std::pow(y, 2) * std::pow(y1, 3) * std::pow(x - y, 2)
    );

    double d2dy = term1 + term2 + term3;

    return  d2dx+d2dy;
}




double lap_schroedinger4D( double a, double b,double x, double y, double u, double v) {
    double t1 =pow(u, 3.0)*pow(x, 3.0)*pow(y, 3.0)*pow(u - 1.0, 3.0)*pow(x - 1.0, 3.0)*pow(y - 1.0, 3.0)*(12.0*a*pow(v, 2.0)*pow(-b + v, 1.0)*pow(v - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*a*pow(v, 3.0)*pow(-b + v, 1.0)*pow(v - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 8.0*a*pow(v, 3.0)*pow(-b + v, 1.0)*pow(v - 1.0, 3.0)*(v - y) + a*pow(v, 3.0)*pow(v - 1.0, 3.0)*(4.0*a*pow(-b + v, 2.0) + 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 6.0*pow(v, 1.0)*pow(v - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) + 18.0*pow(v, 2.0)*pow(v - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*pow(v, 2.0)*pow(v - 1.0, 3.0)*(v - y) + 6.0*pow(v, 3.0)*pow(v - 1.0, 1.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*pow(v, 3.0)*pow(v - 1.0, 2.0)*(v - y) + 2*pow(v, 3.0)*pow(v - 1.0, 3.0))*exp(a*(pow(-b + u, 2.0) + pow(-b + v, 2.0) + pow(-b + x, 2.0) + pow(-b + y, 2.0)));
    double t2 = pow(v, 3.0)*pow(x, 3.0)*pow(y, 3.0)*pow(v - 1.0, 3.0)*pow(x - 1.0, 3.0)*pow(y - 1.0, 3.0)*(12.0*a*pow(u, 2.0)*pow(-b + u, 1.0)*pow(u - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*a*pow(u, 3.0)*pow(-b + u, 1.0)*pow(u - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 8.0*a*pow(u, 3.0)*pow(-b + u, 1.0)*pow(u - 1.0, 3.0)*(u - x) + a*pow(u, 3.0)*pow(u - 1.0, 3.0)*(4.0*a*pow(-b + u, 2.0) + 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 6.0*pow(u, 1.0)*pow(u - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) + 18.0*pow(u, 2.0)*pow(u - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*pow(u, 2.0)*pow(u - 1.0, 3.0)*(u - x) + 6.0*pow(u, 3.0)*pow(u - 1.0, 1.0)*(pow(u - x, 2) + pow(v - y, 2)) + 12.0*pow(u, 3.0)*pow(u - 1.0, 2.0)*(u - x) + 2*pow(u, 3.0)*pow(u - 1.0, 3.0))*exp(a*(pow(-b + u, 2.0) + pow(-b + v, 2.0) + pow(-b + x, 2.0) + pow(-b + y, 2.0)));
    double t3 = pow(u, 3.0)*pow(v, 3.0)*pow(y, 3.0)*pow(u - 1.0, 3.0)*pow(v - 1.0, 3.0)*pow(y - 1.0, 3.0)*(12.0*a*pow(x, 2.0)*pow(-b + x, 1.0)*pow(x - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) - 8.0*a*pow(x, 3.0)*pow(-b + x, 1.0)*(u - x)*pow(x - 1.0, 3.0) + 12.0*a*pow(x, 3.0)*pow(-b + x, 1.0)*pow(x - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + a*pow(x, 3.0)*pow(x - 1.0, 3.0)*(4.0*a*pow(-b + x, 2.0) + 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 6.0*pow(x, 1.0)*pow(x - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) - 12.0*pow(x, 2.0)*(u - x)*pow(x - 1.0, 3.0) + 18.0*pow(x, 2.0)*pow(x - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) - 12.0*pow(x, 3.0)*(u - x)*pow(x - 1.0, 2.0) + 6.0*pow(x, 3.0)*pow(x - 1.0, 1.0)*(pow(u - x, 2) + pow(v - y, 2)) + 2*pow(x, 3.0)*pow(x - 1.0, 3.0))*exp(a*(pow(-b + u, 2.0) + pow(-b + v, 2.0) + pow(-b + x, 2.0) + pow(-b + y, 2.0)));
double t4 = pow(u, 3.0)*pow(v, 3.0)*pow(x, 3.0)*pow(u - 1.0, 3.0)*pow(v - 1.0, 3.0)*pow(x - 1.0, 3.0)*(12.0*a*pow(y, 2.0)*pow(-b + y, 1.0)*pow(y - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) - 8.0*a*pow(y, 3.0)*pow(-b + y, 1.0)*(v - y)*pow(y - 1.0, 3.0) + 12.0*a*pow(y, 3.0)*pow(-b + y, 1.0)*pow(y - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + a*pow(y, 3.0)*pow(y - 1.0, 3.0)*(4.0*a*pow(-b + y, 2.0) + 2.0)*(pow(u - x, 2) + pow(v - y, 2)) + 6.0*pow(y, 1.0)*pow(y - 1.0, 3.0)*(pow(u - x, 2) + pow(v - y, 2)) - 12.0*pow(y, 2.0)*(v - y)*pow(y - 1.0, 3.0) + 18.0*pow(y, 2.0)*pow(y - 1.0, 2.0)*(pow(u - x, 2) + pow(v - y, 2)) - 12.0*pow(y, 3.0)*(v - y)*pow(y - 1.0, 2.0) + 6.0*pow(y, 3.0)*pow(y - 1.0, 1.0)*(pow(u - x, 2) + pow(v - y, 2)) + 2*pow(y, 3.0)*pow(y - 1.0, 3.0))*exp(a*(pow(-b + u, 2.0) + pow(-b + v, 2.0) + pow(-b + x, 2.0) + pow(-b + y, 2.0)));


return t1+t2+t3+t4;

}

double schroedinger4D(double a,double b,double x1, double x2, double y1, double y2) {
    // Calculate the exponential part
    double expPart = exp(a * (pow(x1 - b, 2) + pow(y1 - b, 2) + pow(x2 - b, 2) + pow(y2 - b, 2)));

    // Calculate the dot product
    double dotProduct = (x1-y1)*(x1-y1)+(x2-y2)*(x2-y2);

    // Calculate the polynomial part
    double polyPart = pow(x1, 3) * pow(x1 - 1.0, 3) * pow(y1, 3) * pow(y1 - 1.0, 3) * x2 * pow(x2, 3) * pow(x2 - 1.0, 3) * pow(y2, 3) * pow(y2 - 1.0, 3);

    // Combine all parts
    double result =expPart*dotProduct* polyPart;

    return result;
}

double l2norm_minusVal(double *x, double z) {
    double result=0.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        result += (x[d]-z)*(x[d]-z);
    }
    return  result;
}

double f(double *x, double z, double alpha) {
    double l2_value = l2norm_minusVal(x,z);
    double exp_value = exp(alpha*l2_value);
    return  exp_value;
}

double fp(double *x, double z, double alpha, int dir) {

    double result = 2.0*alpha*(x[dir]-z)*f(x,z,alpha);
    return  result;

}

double fpp(double *x, double z, double alpha, int dir) {

    double result = 2.0*alpha*(1.0+2.0*alpha*(x[dir]-z)*(x[dir]-z))*f(x,z,alpha);
    return  result;

}

double g(double *x) {
    double result=1.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        result*=x[d]*x[d]*x[d]*(x[d]-1)*(x[d]-1)*(x[d]-1);
    }
    return result;
}

double gp(double *x, int dir) {
    double result=1.0;
    for(int d=0;d<DimensionSparseGrid;d++){
        if(d!=dir){
            result*=x[d]*x[d]*x[d]*(x[d]-1)*(x[d]-1)*(x[d]-1);
        }
    }

    double result2=0.0;

    result2 += x[dir]*x[dir]*x[dir]*x[dir]*x[dir]*6.0;
    result2 += x[dir]*x[dir]*x[dir]*x[dir]*(-15.0);
    result2 += x[dir]*x[dir]*x[dir]*12.0;
    result2 +=x[dir]*x[dir]*(-3.0);

    return result2*result;
}

double gpp(double *x, int dir) {
    double result=1.0;
    for(int d=0;d<DimensionSparseGrid;d++){
        if(d!=dir){
            result*=x[d]*x[d]*x[d]*(x[d]-1)*(x[d]-1)*(x[d]-1);
        }
    }

    double result2=0.0;

    result2 += x[dir]*x[dir]*x[dir]*x[dir]*30.0;
    result2 +=x[dir]*x[dir]*x[dir]*(-60.0);
    result2 +=x[dir]*x[dir]*36.0;
    result2 +=x[dir]*(-6.0);

    return result2*result;
}

double fg_pp(double *x, double z, double alpha, int dir) {
    double result=0.0;
    //f''g
    result+= fpp(x,z,alpha,dir)*g(x);

    //2f'g'
    result+= 2.0*fp(x,z,alpha,dir)*gp(x,dir);

    //g''f
    result+= f(x,z,alpha)*gpp(x,dir);

    return  result;
}

double gv_pp(double *x, double z, double alpha, int dir) {
    double result=0.0;
    //f''g
    result+= gpp(x,dir)*v(x);


    //2f'g'
    result+= 2.0*gp(x,dir)*vp(x,dir);


    //g''f
    result+= vpp(x,dir)*g(x);

    return  result;
}


double vpp(double *x, int dir) {


    return 2.0;}

double vp(double *x, int dir) {
    double  value;
    int dim2=DimensionSparseGrid/2;

    if(dir<dim2)value = 2.0*x[dir]-2.0*x[dir+dim2];
    else value =2.0*x[dir]-2.0*x[dir-dim2];

    return value;
}

double v(double *x) {
    double  value=0.0;
    int dim2=DimensionSparseGrid/2;

    for(int d=0; d<dim2;d++){
        value = value + ((x[d]-x[d+dim2])*(x[d]-x[d+dim2]));
    }

    return value;
}

double fg_p(double *x, double z, double alpha, int dir) {
    double result = 0.0;

    //f'g
    result += fp(x,z,alpha,dir)*g(x);


    //g'f
    result += gp(x,dir)*f(x,z,alpha);

    return result;

}

double fgv_pp(double *x, double z, double alpha, int dir) {
    double result=0.0;
    //fg(v'')
    result+= f(x,z,alpha)*g(x)*vpp(x,dir);

    //2(fg)'v'
    result+= 2.0*fg_p(x,z,alpha,dir)*vp(x,dir);

    //(gf)''v
    result+= fg_pp(x,z,alpha,dir)*v(x);

    return  result;
}

double laplacian_fg(double *x, double z, double alpha) {
    double result = 0.0;

    for(int d=0; d<DimensionSparseGrid; d++){
        result +=fg_pp(x,z,alpha,d);

    }
    return  result;
}

double laplacian_gv(double *x, double z, double alpha) {
    double result = 0.0;

    for(int d=0; d<DimensionSparseGrid; d++){
        result +=gv_pp(x,z,alpha,d);

    }
    return  result;
}


double laplacian_fgv(double *x, double z, double alpha) {
    double result = 0.0;

    for(int d=0; d<DimensionSparseGrid; d++){
        result +=fgv_pp(x,z,alpha,d);

    }
    return  result;
}

double sin_f(double *x){
    double val = 1.0;
    for(int d=0; d<DimensionSparseGrid; d++){
        val*=sin(M_PI*x[d]);
    }

    return val;
}
double laplacian_sin_f(double * x){
    auto dim = double(DimensionSparseGrid);
    return dim*M_PI*M_PI*sin_f(x);
};

