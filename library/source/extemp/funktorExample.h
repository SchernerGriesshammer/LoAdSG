#ifndef FUNKTOREXAMPLE_SG_H
#define FUNKTOREXAMPLE_SG_H

#include "extempAlg.h"
#include "funktor.h"
#include <cmath>


/** \defgroup ExampleGroup  ''Often used examples''    
 **/
/* @{  */ 
/**
 * examples of functors which are often used
 * 
 **/
class SinFunctor {
public:
    SinFunctor() = default;

    static double evaluate(double x) {
        return sin(x);
    }
};


class Atan2Functor {
public:
    Atan2Functor() = default;

    static double evaluate(double y, double x) {
        return atan2(y, x);
    }
};

class CosFunctor {
public:
    CosFunctor() = default;

    static double evaluate(double x) {
        return cos(x);
    }
};


class SinHFunctor {
public:
    SinHFunctor() = default;

    static double evaluate(double x) {
        return sinh(x);
    }
};


class LogFunctor {
public:
    LogFunctor() = default;

    static double evaluate(double x) {
        return log(x);
    }
};

class CosHFunctor {
public:
    CosHFunctor() = default;

    static double evaluate(double x) {
        return cosh(x);
    }
};


class ExpFunctor {
public:
    ExpFunctor() = default;

    static double evaluate(double x) {
        //return exp(x);
        return exp(x);

    }
};


class SqrtFunctor {
public:
    SqrtFunctor() = default;

    static double evaluate(double x) {
        return sqrt(x);
    }
};

class ArcTanFunctor {
public:
    ArcTanFunctor() = default;

    static double evaluate(double x) {

        return atan(x);
    }
};

class PowerFunctor {
public:
    PowerFunctor() {
        a = 1.0;
    }

    double evaluate(double x) {
        return pow(x, a);
    }

    double a;
};




class CountFunctor {
public:
    CountFunctor() { count = 0; }

    double evaluate(double x) {
        count++;
        return (double) count;
    }

    inline void resetCount() { count = 0; }

private:
    unsigned long count;
};







/* @} */
  












#endif
