#ifndef FUNKTOR_SG_H
#define FUNKTOR_SG_H

/*#include "extempAlg.h"
#include <math.h>  */

/**
 * Implementation of Functors, see Example in main.cc. (Intern: Welche Funktoren/Funktionen benötigen wir?)
 * 
 * 
 * 
 * 
 **/
template <class A, class Func>
class Exp_Functor : public ExprSparseG<Exp_Functor<A, Func> > {
    const A &a_;
    Func *functor_;
public:
    inline Exp_Functor(const A &a, Func *functor)
            : a_(a), functor_(functor) {}

    inline double getValue(int i, const IndexDimension &I) const {
        return functor_->evaluate(a_.getValue(i, I));
    };


    ExpressionDescription getDescription() const { return a_.getDescription(); }
};

template<class A, class B, class Func>
class Atan2_Functor : public ExprSparseG<Atan2_Functor<A, B, Func> > {
    const A &a_;
    const B &b_;
    Func *functor_;
public:
    inline Atan2_Functor(const A &a, const B &b, Func *functor)
            : a_(a), b_(b), functor_(functor) {}

    inline double getValue(int i, const IndexDimension &I) const {
        return functor_->evaluate(a_.getValue(i, I), b_.getValue(i, I));
    };




    ExpressionDescription getDescription() const { return a_.getDescription(); }
};






#include <cmath>



template<class Func>
class Functor {
public:
    Functor(Func &functor) : functor_(&functor) {};

    template<class A>
    inline Exp_Functor<A, Func>
    operator()(const ExprSparseG<A> &a) const {
        return Exp_Functor<A, Func>(a, functor_);
    }

    template<class A, class B>
    inline Atan2_Functor<A, B, Func>
    operator()(const ExprSparseG<A> &a, const ExprSparseG<B> &b) const {
        return Atan2_Functor<A, B, Func>(a, b, functor_);
    }



private:
    Func *functor_;
};



/* @} */


#endif
