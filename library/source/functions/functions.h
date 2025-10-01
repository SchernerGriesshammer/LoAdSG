//
// Created by to35jepo on 8/1/24.
//

#ifndef RUN_FUNCTIONS_H
#define RUN_FUNCTIONS_H


//example of functions for integration



double schroedinger2D(double a, double b, double x, double y);

double lap_schroedinger2D(double a, double b, double x, double y);


double lap_schroedinger4D(double a, double b,double x1, double x2, double y1, double y2);


double schroedinger4D(double a, double b,double x1, double x2, double y1, double y2);



double l2norm_minusVal(double* x, double z);



double f(double* x, double z, double alpha);
double fp(double* x, double z, double alpha,int dir);
double fpp(double* x, double z, double alpha,int dir);



double g(double* x);

double gp(double *x, int dir);


double gpp(double *x, int dir);

double fg_pp(double*x,double z,double alpha, int dir);

double vpp(double* x,int dir);;

double vp(double* x, int dir);

double v(double* x);

double fg_p(double* x,double z,double alpha,int dir);
double fgv_pp(double*x,double z,double alpha, int dir);

double laplacian_fg(double* x, double z, double alpha);

double laplacian_gv(double* x, double z, double alpha);

double laplacian_fgv(double* x, double z, double alpha);

double sin_f(double* x);

double laplacian_sin_f(double * x);

#endif //RUN_FUNCTIONS_H
