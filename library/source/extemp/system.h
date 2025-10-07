// ------------------------------------------------------------
//
// shift.h
//
// ------------------------------------------------------------
#ifndef SYSTEM_H
#define SYSTEM_H


#include "../myAssert.h"
#include "extempAlg.h"
#include "../sgrid/sparseGrid.h"
#include "../indices/index.h"
#include "vector.h"





////////////////////////////////////////
// 1. class ShiftExpr
// 2. system of equation operators
////////////////////////////////////////


template<class B>
class DWrapSim {
 private:
  // variable
  double* dataTableVector;   

  // expression
  B bo;
 public:
  DWrapSim(const VectorSparseG& a,  const ExprSparseG<B>& b) : dataTableVector(a.dataTableVector), bo(b) {  }
  void RunValue(int i, const IndexDimension& I) const {
       dataTableVector[i]=bo.getValue(i, I);
       // geht noch icht!!
  };


  ExpressionDescription getDescription() const {
        return ExpressionDescription(bo.getDescription()); 
  }	
  
};


/* DExprAnd */
template<class A, class B>
class ExprAndSparseG : public ExprSparseG<ExprAndSparseG<A, B> > {
  const A a_;
  DWrapSim<B> b_;
 public:
    ExprAndSparseG(const A& a, DWrapSim<B>& b)  : a_(a), b_(b) {}
	 
    double getValue(int i, const IndexDimension& I) const { 
        b_.RunValue(i,I);
        return a_.getValue(i,I);  
    }
 
    ExpressionDescription getDescription() const {
        return ExpressionDescription(a_.getDescription(),b_.getDescription()); 
    }		 
    AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		    
    
    // at the and of an iteration of an expression
    //void clean() const { a_.clean(); b_.clean(); };
};

/*
// DExprAndSim
template<class A, class B>
class DExprAndSim : public DWrapSim<DExprAndSim<A, B> > {
  A a_;
  B b_;
 public:
    DExprAndSim(const A& a, const B& b)  : a_(a), b_(b) {}
	 
    void RunValue(double* data, const IndexDimension& I) const { 
        b_.RunValue(data,I);
        a_.RunValue(data,I);  
    }
    void RunDataTableVector(int i) const { 
        b_.RunDataTableVector(i);
        a_.RunDataTableVector(i);  
    }   
    ExpressionDescription getDescription() const {
        return ExpressionDescription(a_.getDescription(),b_.getDescription()); 
    }		 
    AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		    
    
    // at the and of an iteration of an expression
    //void clean() const { a_.clean(); b_.clean(); };
};
*/
 
class DExprLIT : public ExprSparseG<DExprLIT> {		
    public:			
	DExprLIT(const double a_): a(a_) {}
		
        double getValue(int i, const IndexDimension& I) const { return a; }
    
        ExpressionDescription getDescription() const { return ExpressionDescription(false); }
        AdaptiveSparseGrid* getSparseGrid() const { return NULL; }	//  Das ist natürlich komisch so! Riccarda: Ist das ein Problem???				
    private:	
	double a;
};	
		

///////////////////////////////////////////////////
// system of equation operators
///////////////////////////////////////////////////


/* @{ */  

template<class A, class B>
ExprAndSparseG< A, B >
operator&(const ExprSparseG<A>& a, DWrapSim<B> b) {
  return ExprAndSparseG< A, B > (a,b);
}


/*
template<class A, class B>
ExprAndSparseG< A, B >
operator&(A a, DWrapSim<B> b) {
  return ExprAndSparseG< A, B > (a,b);
}
*/

/*
template<class A, class B>
DExprAndSim< A, B >
operator&(const  DWrapSim<A>& a, const DWrapSim<B>& b) {
  return DExprAndSim< A, B > (a,b);
}
*/

/*
template<class B>
ExprAndSparseG< VectorSparseG, B >
operator&(VectorSparseG& a, const DWrapSim<B>& b) {
  return ExprAndSparseG< VectorSparseG, B > (a,b);
}
*/

template<class B>
ExprAndSparseG< DExprLIT, B >
operator&(double a, const DWrapSim<B>& b) {
  return ExprAndSparseG< DExprLIT, B > (DExprLIT(a),b);
}


/*
template<class B>
DExpr<DExprAnd<Local_var,  DWrapSim<B> > >
operator&(Local_var &a, const DWrapSim<B>& b)
{
  typedef DExprAnd<Local_var,  DWrapSim<B> >  ExprT;
  return DExpr<ExprT>(ExprT(a,b));
}
*/

template<class B> 
DWrapSim<B>
operator== (VectorSparseG& a, const ExprSparseG<B>& b_) {
  return DWrapSim<B>(a,b_); 
}


/*
DWrapSim<DExprSimLit>
operator== (Variable& a, double value);


template<class B> 
DWrapSim<DExprSimloc<B> > 
operator== (const Local_var& a, const  B& b_)
{
  typedef DWrapSim<DExprSimloc<B> > ExprT;
  return ExprT(DExprSimloc<B>(a,b_)); 
}

DWrapSim<DExprSimlocLit>
operator== (const Local_var& a, double x);

DWrapSim<DExprSimloc<DExprVAR> >
operator== (const Local_var& a, Variable& b);
*/


//@}


#endif
