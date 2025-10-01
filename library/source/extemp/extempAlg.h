/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef EXTEMP_SPARSE_double_H
#define EXTEMP_SPARSE_double_H
#include "../sgrid/sparseGrid.h"
#include "../indices/index.h"

// ExprSparseG makes constructor "_a(a)" possible
// Wrapper Class
// SiWir Skript -> S.63 
template <class A> 
struct ExprSparseG {
    // operator () is only defined for const class A& element und returns const
	inline operator const A&() const {
		return *static_cast<const A*> ( this );
	}
};

	

//----------------------------------------------------------------------------

/*
Exp_Function_Coord
 */
class ExpressionDescription {
  public:
    ExpressionDescription(bool needsIndex_) : needsIndex(needsIndex_) {};
    ExpressionDescription(const ExpressionDescription& a, const ExpressionDescription& b) {
                          needsIndex = a.needsIndex || b.needsIndex;  };
    bool isIndexNeeded() { return needsIndex; };
  private:  
    bool needsIndex;
};
  
//----------------------------------------------------------------------------
template <class A, class B>
class AddSparseG : public ExprSparseG<AddSparseG<A, B> > {
		const A& a_;
		const B& b_;
		
	public:
		inline AddSparseG ( const A& a, const B& b ) : a_ ( a ), b_ ( b ) {};

		//inline double getValue(double* data, const IndexDimension& I) const { 
		//    return a_.getValue(data,I) + b_.getValue(data,I);  
		//}
		
		
		inline double getValue(int i, const IndexDimension& I) const { 
		    return (a_.getValue(i,I) + b_.getValue(i,I));  
		};
		
	
		
		ExpressionDescription getDescription() const {
		    return ExpressionDescription(a_.getDescription(),b_.getDescription()); 
		};
		AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); };
};

//----------------------------------------------------------------------------
template <class A, class B>
inline AddSparseG<A, B> operator+ ( const ExprSparseG<A>& a, const ExprSparseG<B>& b ) {
	return AddSparseG<A, B> ( a, b );
}






//----------------------------------------------------------------------------
template <class A>
class CAddSparseG : public ExprSparseG<CAddSparseG<A> > {
		const A& a_;
		double b_;
		
	public:
		inline CAddSparseG(const A& a, double b ) : a_ ( a ), b_ ( b ) {}


		inline double getValue(int i, const IndexDimension& I) const { 
		    return a_.getValue(i,I) + b_;  
		}
		
	
	
        
		ExpressionDescription getDescription() const {
		    return a_.getDescription(); 
		}		
		AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		
};

//----------------------------------------------------------------------------
template <class A>
inline CAddSparseG<A>
operator+(const ExprSparseG<A>& a, double b ) {
	return CAddSparseG<A> ( a, b );
}

//----------------------------------------------------------------------------
template <class A>
inline CAddSparseG<A>
operator+(double b, ExprSparseG<A>& a ) {
	return CAddSparseG<A> ( a, b );
}

//----------------------------------------------------------------------------
template <class A>
inline CAddSparseG<A>
operator-(ExprSparseG<A>& a, double b ) {
	return CAddSparseG<A> ( a, -1.0*b );
}


//----------------------------------------------------------------------------
template <class A, class B>

class SubSparseG : public ExprSparseG<SubSparseG<A, B> > {
		const A& a_;
		const B& b_;
		
	public:
		 inline SubSparseG ( const A& a, const B& b ) : a_ ( a ), b_ ( b ) {}
		
		 inline double getValue(int i, const IndexDimension& I) const { 
		    return a_.getValue(i,I) - b_.getValue(i,I);  
		 }

        
		 ExpressionDescription getDescription() const {
		    return ExpressionDescription(a_.getDescription(),b_.getDescription()); 
		 }		 
		 AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		 
};

//----------------------------------------------------------------------------
template <class A, class B>
inline SubSparseG<A, B> operator- ( const ExprSparseG<A>& a, const ExprSparseG<B>& b ) {
	return SubSparseG<A, B> ( a, b );
}

//----------------------------------------------------------------------------
template < class B>
class CSubSparseG : public ExprSparseG<CSubSparseG< B> > {
		double a_;
		const B& b_;
		
	public:
 		 inline CSubSparseG ( double a, const B& b ) : a_ ( a ), b_ ( b ) {}

		 inline double getValue(int i, const IndexDimension& I) const { 
		    return a_ - b_.getValue(i,I);  }

		 ExpressionDescription getDescription() const {
		    return b_.getDescription(); 
		 }		 
		 AdaptiveSparseGrid* getSparseGrid() const { return  b_.getSparseGrid(); }		 
};

//----------------------------------------------------------------------------
template <class B>
inline CSubSparseG< B> operator- ( double a, const ExprSparseG<B>& b ) {
	return CSubSparseG< B> ( a, b );
}


//----------------------------------------------------------------------------
template <class A>

class MinusSparseG : public ExprSparseG<MinusSparseG<A> > {
		const A& a_;
		
	public:
		inline MinusSparseG ( const A& a ) : a_ ( a ) {}
		
		
		 inline double getValue(int i, const IndexDimension& I) const { 
		   return -1.0*a_.getValue(i,I);  
		 }

		 ExpressionDescription getDescription() const {
		    return a_.getDescription(); 
		 }		 
		 AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		 
};

//----------------------------------------------------------------------------
template <class A>
inline MinusSparseG<A> operator- ( const ExprSparseG<A>& a ) {
	return MinusSparseG<A> ( a );
}


//----------------------------------------------------------------------------
template <class A>

class TimesSparseG : public ExprSparseG<TimesSparseG<A> > {
		const A& a_;
		double b_;
		
	public:
		inline TimesSparseG ( const A& a, const double& b ) : a_ ( a ), b_ ( b ) {}
		

		inline double getValue(int i, const IndexDimension& I) const { 
		    return a_.getValue(i,I) * b_;  
		}
		
	
		ExpressionDescription getDescription() const {
		    return a_.getDescription(); 
		}		
		AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		
};


//----------------------------------------------------------------------------
template <class A>
inline TimesSparseG<A>
operator* ( const ExprSparseG<A>& a, double b ) {
	return TimesSparseG<A> ( a, b );
}

//----------------------------------------------------------------------------
template <class A>
inline TimesSparseG<A>
operator* ( double b, const ExprSparseG<A>& a ) {
	return TimesSparseG<A> ( a, b );
}


//----------------------------------------------------------------------------
template < class B>

class DivSparseG : public ExprSparseG<DivSparseG< B> > {
		const double a_;
		const B& b_;
		
	public:
		inline DivSparseG ( const double& a, const B& b ) : a_ ( a ), b_ ( b ) {}

		inline double getValue(int i, const IndexDimension& I) const { 
		   return a_ / b_.getValue(i,I);  
		}
		

		ExpressionDescription getDescription() const {
		    return b_.getDescription(); 
		}		
		AdaptiveSparseGrid* getSparseGrid() const { return  b_.getSparseGrid(); }		
};

//----------------------------------------------------------------------------
template <class A>
inline TimesSparseG<A> operator/ ( const ExprSparseG<A>& a, double b ) {
	return TimesSparseG<A> ( a, 1.0 / b );
}

//----------------------------------------------------------------------------

template <class B>
inline DivSparseG< B> operator/ ( double a, const ExprSparseG<B>& b ) {
	return DivSparseG< B> ( a, b );
}


//----------------------------------------------------------------------------
template <class A, class B>

class VTimesSparseG : public ExprSparseG<VTimesSparseG<A, B> > {
    const A &a_;
    const B &b_;

public:
    inline VTimesSparseG(const A &a, const B &b) : a_(a), b_(b) {}

    inline double getValue(int i, const IndexDimension &I) const {
        if (a_.getValue(i, I) == 0 || b_.getValue(i, I) == 0) return 0.0;
        return a_.getValue(i, I) * b_.getValue(i, I);
    };


    ExpressionDescription getDescription() const {
        return ExpressionDescription(a_.getDescription(), b_.getDescription());
    };

    AdaptiveSparseGrid *getSparseGrid() const { return a_.getSparseGrid(); }
};


//----------------------------------------------------------------------------
template <class A, class B>
inline VTimesSparseG<A, B> operator* ( const ExprSparseG<A>& a, const ExprSparseG<B>& b ) {
	return VTimesSparseG<A, B> ( a, b );
}


//----------------------------------------------------------------------------
template <class A, class B>

class VDivSparseG : public ExprSparseG<VDivSparseG<A, B> > {
		const A& a_;
		const B& b_;
		
	public:
		inline VDivSparseG ( const A& a, const B& b ) : a_ ( a ), b_ ( b ) {}

		 inline double getValue(int i, const IndexDimension& I) const { 
		    return a_.getValue(i,I) / b_.getValue(i,I);  
		 };
		 

		 ExpressionDescription getDescription() const {
		    return ExpressionDescription(a_.getDescription(),b_.getDescription()); 
		 }		 
		 AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }		 
};

//----------------------------------------------------------------------------
template <class A, class B>
inline VDivSparseG<A, B> operator/ ( const ExprSparseG<A>& a, const ExprSparseG<B>& b ) {
	return VDivSparseG<A, B> ( a, b );
}

//----------------------------------------------------------------------------

#endif // EXTEMP_H
