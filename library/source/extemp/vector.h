/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/
 
// ------------------------------------------------------------
//
// vector.h
//
// ------------------------------------------------------------

#ifndef VECTOR_SG_H
#define VECTOR_SG_H
#define  testNow

#include "../abbrevi.h"
#include "../myAssert.h"
#include "../indices/index.h"

#include "../primes/prime.h"
#include "../sgrid/depth.h"
#include "../sgrid/sparseGrid.h"
#include "extempAlg.h"


#include "../sgrid/ListOfDepthOrderedGrids.h"
#include "../sgrid/multilevelSparseGrid.h"



#include "../mympi.h"
#include "../sgrid/multiDepthHashGrid.h"






/////////////////////////////////////////////
//  Restriction 
/////////////////////////////////////////////

/*
template <class A>
class ExprSparseG_Restriction_Index {
   public:
     ExprSparseG_Restriction_Index(const ExprSparseG<A>& a, const IterationDepth& iterDepth_) 
           : a_(a), iterDepth(iterDepth_) { }

    // inline double getValue(double* data, const IndexDimension& I) const { 
	//	    return a_.getValue(data,I);  
     //}
     inline double getValue(int i, const IndexDimension& I)const{
          return a_.getValue(i,I); 
    }
     
     

     ExpressionDescription getDescription() const {
		    return a_.getDescription(); 
     }	
     
     AdaptiveSparseGrid_Base* getSparseGrid() const { return  a_.getSparseGrid(); }
     
     IterationDepth getIterDepth() const { return iterDepth; }
     
   public:
     //const reference to an object of class A  
     const A& a_;
     const IterationDepth iterDepth;
};

template <class A>
inline ExprSparseG_Restriction_Index<A> operator| ( const ExprSparseG<A>& a, const IterationDepth& iterDepth) {
	return ExprSparseG_Restriction_Index<A>(a,iterDepth);
}


template <class A>
inline ExprSparseG_Restriction_Index<A> operator| (const ExprSparseG<A>& a,int i) {
    myAssert(i > -1);
	IterationDepth iterDepth(i);
    return ExprSparseG_Restriction_Index<A>(a,iterDepth);
}
*/

/////////////////////////////////////////////////////
// Restriction on a subgrid of certain depth
/////////////////////////////////////////////////////


template <class A>
class ExprSparseG_Restriction_SubgridDepth {

public:
    ExprSparseG_Restriction_SubgridDepth(const ExprSparseG<A> &a, const SubgridFixedDepth &iterObj_)
            : a_(a), iterObj(iterObj_) {}

    inline double getValue(int i, const IndexDimension &I) const {

        return a_.getValue(i, I);
    }

    ExpressionDescription getDescription() const {
        return a_.getDescription();
    }

    AdaptiveSparseGrid *getSparseGrid() const { return a_.getSparseGrid(); }

    SubgridFixedDepth getSubgridIterator() const { return iterObj; }

public:
    //const reference to an object of class A
    const A &a_;
    const SubgridFixedDepth iterObj;
};

template<class A>
inline ExprSparseG_Restriction_SubgridDepth<A> operator|(const ExprSparseG<A> &a, const SubgridFixedDepth &iterObj) {

    return ExprSparseG_Restriction_SubgridDepth<A>(a, iterObj);
}
//////////////////////////////////////////////////////
//Restriction on fixed Index


template<class A>
class ExprSparseG_Restriction_FixedIndex {

public:
    ExprSparseG_Restriction_FixedIndex(const ExprSparseG<A> &a, const IndexDimension &iterObj_)
            : a_(a), iterObj(iterObj_) {}

    inline double getValue(int i, const IndexDimension &I) const {
        return a_.getValue(i, I);
    }

    ExpressionDescription getDescription() const {
        return a_.getDescription();
    }

    AdaptiveSparseGrid *getSparseGrid() const { return a_.getSparseGrid(); }

    IndexDimension getIndex() const { return iterObj; }

public:
    //const reference to an object of class A
    const A &a_;
    const IndexDimension iterObj;
};

/*
template<class A>
inline ExprSparseG_Restriction_FixedIndex<A> operator|(const ExprSparseG<A> &a, const IndexDimension &iterObj) {
    return ExprSparseG_Restriction_FixedIndex<A>(a, iterObj);
}
*/




class DoubleSubgridExpr {

public:
    DoubleSubgridExpr(double d_, SubgridFixedDepth itObj_) : d(d_), itObj(itObj_) {};

    SubgridFixedDepth getSubgrid() const { return itObj; };

    double getValue() const { return d; };
private:
    double d;
    SubgridFixedDepth itObj;

};


class DoubleDepthExpr {

public:
    DoubleDepthExpr(double d_, Depth itObj_) : d(d_), itObj(itObj_) {};

    Depth getDepth() const { return itObj; };

    double getValue() const { return d; };
private:
    double d;
    Depth itObj;

};


inline DoubleSubgridExpr
operator|(double d, const SubgridFixedDepth &itObj) {
    return DoubleSubgridExpr(d, itObj);
}


inline DoubleDepthExpr
operator|(double d, const Depth &itObj) {
    return DoubleDepthExpr(d, itObj);
}


class DoubleIndexExpr {

public:
    DoubleIndexExpr(double d_, IndexDimension itObj_) : d(d_), itObj(itObj_) {};

    IndexDimension getIndex() const { return itObj; };

    double getValue() const { return d; };
private:
    double d;
    IndexDimension itObj;

};


inline DoubleIndexExpr
operator|(double d, const IndexDimension &itObj) {
    return DoubleIndexExpr(d, itObj);
}







/////////////////////////////////////////////////////
// Restriction with Level Object





///////////////////////////
// 2. vector
///////////////////////////

template<class B>
class DWrapSim;


 class VectorSparseG : public ExprSparseG<VectorSparseG> {
     template<class B>
     friend
     class DWrapSim;


 public:
     int length=0;
     VectorSparseG(){};

     VectorSparseG(AdaptiveSparseGrid &sg);

     VectorSparseG(AdaptiveSparseGrid *sg);

     VectorSparseG(VectorSparseG &u);


     ~VectorSparseG();


     AdaptiveSparseGrid *getSparseGrid() { return sparseGrid; }

     double *getDatatableVector() { return dataTableVector; }


     template<class A>
     void operator=(const ExprSparseG<A> &a);


     template<class A>
     void operator=(const ExprSparseG_Restriction_SubgridDepth<A> &a);

     template<class A>
     void operator=(const ExprSparseG_Restriction_FixedIndex<A> &a);

     void operator=(const VectorSparseG &v);

     void operator=(VectorSparseG &v);

     inline void operator=(double x);

     bool operator==(const VectorSparseG &v);

     void setMultiLevelValues2(MultiLevelVector &a, Depth &T);

     void addMultiLevelValues2(MultiLevelVector &a, Depth &T);
     void setMultiLevelValues(MultiLevelVector &a, Depth &T);

     inline void operator=(const DoubleSubgridExpr &a);

     inline void operator=(const DoubleIndexExpr &a);

     inline bool workonindex(unsigned long i) {

         if ((sparseGrid->WorkOnHangingNodes) || (!(sparseGrid->WorkOnHangingNodes) && sparseGrid->getActiveTable()[i]))
             return true;

         return false;
     }

     ExpressionDescription getDescription() const { return ExpressionDescription(false); }

     AdaptiveSparseGrid_Base *getSparseGrid() const { return sparseGrid; }

     /**
      * Aussgabe Form:
      *
      *      --x--
      *     |
      *     y
      *     |
      *
      *   x: Dimension 0
      *   y: Dimension 1
      *
      */
     void Print(int level);

     void PrintDouble(int level);

     void PrintDoubleTwoD(int level, Depth T);


     void PrintIntTwoD(int level);

     void PrintDoubleTwoD(int level);


     void PrintDoubleOneD(int level);

     void PrintIntOneD(int level);

    void Print_vtk(std::ostream &Datei);

    void Print_gnu(string name);

    void PrintIndexTwoD(int d, int level);


     inline double getValue(int i, const IndexDimension &I) const { return dataTableVector[i]; };

     inline double getValue(IndexDimension &I) const {
         unsigned long k;
         if (sparseGrid->occupied(k, I)) return dataTableVector[k];
         return 0.0;

     };

     inline double getValue(unsigned long k) const {

         return dataTableVector[k];


     };

     inline double getValue(const IndexDimension &I) const {
         unsigned long k;
         if (sparseGrid->occupied(k, I)) return dataTableVector[k];
         return 0.0;

     };

     inline double getValue(unsigned long k) {
         return dataTableVector[k];
     }

     inline double getValue(const IndexDimension &I, const SingleDepthHashGrid *singleDepthHashGrid) const {
         unsigned long k;
         if (singleDepthHashGrid->occupied(k, I)) return dataTableVector[k];
         return 0.0;

     };

     inline bool setValue(IndexDimension &I, double x) const {
         unsigned long k;
         if (sparseGrid->occupied(k, I)) {
             dataTableVector[k] = x;
             return true;
         }
         return false;

     };

     inline bool setValue(unsigned long k, double x) const {

         dataTableVector[k] = x;
         return true;


     };

     inline bool addToValue(unsigned long k, double x) const {
         #pragma omp atomic
         dataTableVector[k] += x;
         return true;
     };

     inline bool addToValueNoOMP(unsigned long k, double x) const {
         dataTableVector[k] += x;
         return true;
     };
     inline bool addToValue(IndexDimension &I, double x) const {
         unsigned long k;
         if (sparseGrid->occupied(k, I)) {
             dataTableVector[k] += x;
             return true;
         }
         return false;

     };



     ///// MPI

     void Broadcast(int rank);

     bool mpi_doit();

     void sendTo(int torank);

     void ReduceSum(int rank);


     Process *getProcess() { return sparseGrid->mpi; }

 protected:
     AdaptiveSparseGrid *sparseGrid;
     double *dataTableVector;
 private:

     bool constructed=false;


 };





///////////////////////////
// 3. Other functions, max, ...
///////////////////////////

template <class A, class B>
double product(const ExprSparseG<A>& a, const ExprSparseG<B>& b );


template <class A>
double Maximum ( const ExprSparseG<A>& a );


template <class A>
double Minimum ( const ExprSparseG<A>& a );



//////////////////////////////////////////
// and it's implementation
template<typename A, typename B> double product(const ExprSparseG<A>& a, const ExprSparseG<B>& b)
{
    const A& ao (a);
    const B& bo (b);
    
    AdaptiveSparseGrid_Base* sparseGridA = ao.getSparseGrid();
    AdaptiveSparseGrid_Base* sparseGridB = bo.getSparseGrid();
    
    // todo check if a,b have same sparse grids
    
     unsigned long endIndex   = sparseGridA->maximalOccupiedSecondTable;
     //double* dataTable        = sparseGrid->dataTable;
     dataInteger* secondTable = sparseGridA->secondTable;
     //unsigned int numberOfData = sparseGrid->numberOfData;

     
    
    if(ao.getDescription().isIndexNeeded()|| bo.getDescription().isIndexNeeded()) {
        double value = 0.0;
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(secondTable[i]!=0) {
               
            //double x = ao.getValue(&dataTable[i * numberOfData],sparseGrid->getIndexOfTable(i));
            double ax = ao.getValue(i,sparseGridA->getIndexOfTable(i));
            double bx = bo.getValue(i,sparseGridB->getIndexOfTable(i));
     
	        value = value + ax*bx;
            }
        }        
        return value;
     }
     else {
        double value = 0.0;
        IndexDimension Idummy;
        for(unsigned long i = 0;i < endIndex; ++i) {
            if(secondTable[i]!=0) {
                double ax = ao.getValue(i,Idummy);
                double bx = bo.getValue(i,Idummy);
     
                value = value + ax*bx;
            }
        }        
        return value;  
     }
}











//////////////////////////////////////////////////////////////////////////////////////////////
  

template <class A>
void VectorSparseG::operator=(const ExprSparseG<A>& a) {

    if (mpi_doit()){
     const A& ao ( a );
     
     unsigned long endIndex    = sparseGrid->maximalOccupiedSecondTable;
     //double* dataTable         = sparseGrid->dataTable;
     dataInteger* secondTable  = sparseGrid->secondTable;
     //unsigned int numberOfData = sparseGrid->numberOfData;
     
     if(ao.getDescription().isIndexNeeded()) {
        
        for(unsigned long i = 0;i < endIndex; ++i) {

            if (workonindex(i))
                dataTableVector[i] = ao.getValue(i, sparseGrid->getIndexOfTable(i));


        }
     }
     else {
       IndexDimension Idummy;
       for(unsigned long i = 0;i < endIndex; ++i) {
         if(secondTable[i]!=0)
            //dataTable[number + i * numberOfData] = ao.getValue(&dataTable[i * numberOfData],Idummy);
             
              dataTableVector[i]=ao.getValue(i,Idummy);
       }
     }
    }
}



//test
//IndexDimension I;
//I = iteratorSG.getCurrentPoint();
//cout << "   IteratorSG::hasNext() " << I.coordinate(0)
//     << " " << I.coordinate(1) << " i: " << i << endl;






template <class A>
void VectorSparseG::operator=(const ExprSparseG_Restriction_SubgridDepth<A>& a) {

    if (mpi_doit()) {
        //unsigned long endIndex = sparseGrid->maximalOccupiedSecondTable;
        //double* dataTable         = sparseGrid->dataTable;
        dataInteger *secondTable = sparseGrid->secondTable;
        //unsigned int numberOfData = sparseGrid->numberOfData;

        //hier wird gesamtes subgrid kopiert. braucht man das?
        SubgridFixedDepth subgrid = a.getSubgridIterator();
        SubgridFixedDepth::iterator itobj(subgrid);

        for (bool weiter = true; weiter; weiter = itobj.next()) {

            unsigned long i = itobj.geti();
            if (a.getDescription().isIndexNeeded()) {
                if (secondTable[i] != 0)
                    if (workonindex(i))
                        dataTableVector[i] = a.getValue(i, itobj.getPoint());
            } else {
                IndexDimension Idummy;
                if (secondTable[i] != 0)
                    if (workonindex(i))
                        dataTableVector[i] = a.getValue(i, Idummy);
            }
        }
    }
}


template<class A>
void VectorSparseG::operator=(const ExprSparseG_Restriction_FixedIndex<A> &a) {

    if (mpi_doit()) {
        //unsigned long endIndex = sparseGrid->maximalOccupiedSecondTable;
        //double* dataTable         = sparseGrid->dataTable;
        //dataInteger *secondTable = sparseGrid->secondTable;
        //unsigned int numberOfData = sparseGrid->numberOfData;


        IndexDimension Index = a.getIndex();
        unsigned long k;
        if (sparseGrid->occupied(k, Index)) {
            if (workonindex(k))
                dataTableVector[k] = a.getValue(k, Index);
        }


        /*   for (int i = 0; i <= endIndex; i++) {

               IndexDimension localIndex = sparseGrid->getIndexOfTable(i);
               if (localIndex == Index)
                   if (a.getDescription().isIndexNeeded()) {
                       if (workonindex(i))
                           dataTableVector[i] = a.getValue(i, Index);
                   } else {
                       IndexDimension Idummy;
                       if (workonindex(i))
                           dataTableVector[i] = a.getValue(i, Idummy);
                   }
           }*/
    }
}


void VectorSparseG::operator=(const DoubleSubgridExpr &a) {

    if (mpi_doit()) {
        SubgridFixedDepth subgrid = a.getSubgrid();
        SubgridFixedDepth::iterator itobj(subgrid);

        double value = a.getValue();

        for (bool weiter = true; weiter; weiter = itobj.next()) {
            unsigned long i = itobj.geti();
            if (workonindex(i))
                dataTableVector[i] = value;
        }
    }
}


void VectorSparseG::operator=(const DoubleIndexExpr &a) {


    if (mpi_doit()) {
        IndexDimension Index = a.getIndex();
        unsigned long k;
        if (sparseGrid->occupied(k, Index)) {
            if (workonindex(k))
                dataTableVector[k] = a.getValue();
        }

        /* double value = a.getValue();

         for (int i = 0; i <= sparseGrid->maximalOccupiedSecondTable; i++) {
             IndexDimension localIndex = sparseGrid->getIndexOfTable(i);
             if (Index == localIndex)
                 if (workonindex(i))
                     dataTableVector[i] = value;
         }*/
    }

}






/////////////////////////////////////////////////////////

inline void VectorSparseG::operator=(double x) {

    if (mpi_doit()) {
        unsigned long endIndex = sparseGrid->maximalOccupiedSecondTable;
        //double* dataTable        = sparseGrid->dataTable;
        dataInteger *secondTable = sparseGrid->secondTable;
        //unsigned int numberOfData = sparseGrid->numberOfData;



        for (unsigned long i = 0; i <= endIndex; ++i) {
            if (secondTable[i] != 0) {
                if (sparseGrid->workonindex(i))
                    dataTableVector[i] = x;
            }
        }
    }

}



#endif //VECTOR_SG_H



