// ------------------------------------------------------------
//
// shift.h
//
// ------------------------------------------------------------
#ifndef SHIFT_SG_H
#define SHIFT_SG_H


#include "../myAssert.h"
#include "extempAlg.h"
#include "../sgrid/sparseGrid.h"
#include "../indices/index.h"
#include "vector.h"




//todo: welche functions inline?



////////////////////////////////////////
// 1. class ShiftExpr
// 2. class ShiftOperator
//
////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// 1.) class ShiftExpr
////////////////////////////
/**
 * 
 *Expression Template for a ShiftOperator 
 *
 *
 * 
 * 
 **/


template <class A>
class ShiftExpr : public ExprSparseG<ShiftExpr<A> >{
friend VectorSparseG;
   
   
    const A& a_;
    int d_; //>>> Dimension
    Direction s_; //>>> Indicates wether to shift left or right
    int level_;
    bool hier;
    
public: 
    ShiftExpr(const A& a, int d, Direction s):a_(a),d_(d), s_(s){
        hier = true;
    }
    ShiftExpr(const A& a, int d, Direction s, int level):a_(a),d_(d), s_(s), level_(level){
        hier = false;
    }
    
    
    inline double getValue(int i, const IndexDimension& K)const;

    
    ExpressionDescription getDescription() const { return ExpressionDescription(true); }
    AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }
};






template<typename A>
inline double ShiftExpr<A>::getValue(int i, const IndexDimension& K)const{
      
        if(hier == true ){
        // do hierarchical shift    
         
            AdaptiveSparseGrid_Base* sparseGrid = a_.getSparseGrid();
            //unsigned long endIndex   = sparseGrid->maximalOccupiedSecondTable;
            
            
            dataInteger* secondTable = sparseGrid->getSecondTable();
            dataInteger* primeTable = sparseGrid->getPrimeTable();
            
            if(secondTable[i]!=0) {
                IndexDimension I = K;
                bool nonext = true;
               
                //if(I.isNotAtBoundary() == false) {
                    //return 0;  
                    //return a_.getDataTableVector(i);  
                //}
                IndexDimension II;
                unsigned long kk;
                while(nonext){
                    if (s_ == Left) {
                        if (sparseGrid->occupied(kk, I.nextLeft(d_)))
                            II = I.nextLeft(d_);
                        else return 0;
                    } else {
                        if (sparseGrid->occupied(kk, I.nextRight(d_)))
                            II = I.nextRight(d_);
                        else return 0;
                    }

                    unsigned long indArray =  sparseGrid->hash(II);
                    unsigned long k = primeTable[indArray];
                           
                    if(k == 0) { 
                      
                            return a_.getValue(i,K);  
                       
               
                    } else { // anaylse what is going on 
                    k = k-1; // shift since data are stored with shift
                        if(II == sparseGrid->getIndexOfTable(k)) {
                            return a_.getValue(k,K);
                        } else { // search in second table;
                            unsigned long iNext = secondTable[k];
                            while(iNext > 1) {
                            k = iNext - 2;
                                if(II == sparseGrid->getIndexOfTable(k)) {
                                return a_.getValue(k,K) ;
                                nonext = false;
                                }
                            iNext = secondTable[k];
                            }
                        }
                    }
                }
            }
        } else {
            AdaptiveSparseGrid_Base* sparseGrid = a_.getSparseGrid();
            //unsigned long endIndex   = sparseGrid->maximalOccupiedSecondTable;
                        
            dataInteger* secondTable = sparseGrid->getSecondTable();
            dataInteger* primeTable = sparseGrid->getPrimeTable();
            
            int maxdepth = sparseGrid->getMaxDepth(d_);
           
            
            //assert( level_ <= maxdepth);
            maxdepth = level_;

            if(secondTable[i]!=0) {
                IndexDimension I = sparseGrid->getIndexOfTable(i);
                bool nonext = true;
                // todo: andere If Bedingung
                if(I.isNotAtBoundary() == false) {
                    return 0;  
                    //return a_.getDataTableVector(i);  
                }
                IndexDimension II;
                while(nonext){
                    if (s_ == Left){ 
                        II = I.nextLeft(d_,maxdepth);
                    }else{
                        II = I.nextRight(d_,maxdepth);
                    }

                    unsigned long indArray =  sparseGrid->hash(II);
                    unsigned long k = primeTable[indArray];
                           
                    if(k == 0) { 
                        maxdepth= maxdepth-1;
                        if (maxdepth == 0){
                            return a_.getValue(i,K);  
                        }
               
                    } else { // anaylse what is going on 
                    k = k-1; // shift since data are stored with shift
                        if(II == sparseGrid->getIndexOfTable(k)) {
                            return a_.getValue(k,K);
                        } else { // search in second table;
                            unsigned long iNext = secondTable[k];
                            while(iNext > 1) {
                            k = iNext - 2;
                                if(II == sparseGrid->getIndexOfTable(k)) {
                                return a_.getValue(k,K) ;
                                nonext = false;
                                }
                            iNext = secondTable[k];
                            }
                        }
                    }
                }
            }
        }
        return 0.0;
}







 
//////////////////////////////////////////////////////////////////////////////
// 2.) class ShiftOperator
//////////////////////
 
        

/* @{ */        
/**
 *Example:\verbatim
   ShiftOperator Left(0,LeftShift);
   VectorSparseG a(AdaptiveSparseGrid);
   VectorSparseG b(AdaptiveSparseGrid);
   ...
   a = L(b);    
   \endverbatim
 *
 *Then a contains the hierarchical to the left shifted values of b. 
 *
 *
 **/
class ShiftOperator{
public:
/**
 * @param int d: dimension 
 * @param LeftRight s: shift left or right here we have a hierarchical shift, i.e. we get the "left" father and "right" father
 * 
 * **/
    ShiftOperator(int d, Direction s):d_(d),s_(s){hier =true;}
        
 /**
 * @param d: dimension 
 * @param s: shift left (= LeftShift) or right (=RightShift)
 * @param level: shift on a fixed level
 * 
 * 
 * **/       
    ShiftOperator(int d, Direction s, int level):d_(d),s_(s), level_(level){hier =false;}

    template <class A>
    ShiftExpr<A> operator()(ExprSparseG<A>& a){
       if ( hier == true ){
           return ShiftExpr<A> (a,d_,s_);
        } else {
            return ShiftExpr<A> (a,d_,s_,level_);
        } 
    }
    
private:
    int d_;
    Direction s_;
    
    bool hier;
    int level_;
};
//@}

/////////////////////////////////////////////////////////
// 3.) ShiftExprDimension
////////////////////////////////////////////////////////
enum Compass {
    North, NorthEast, East, SouthEast, South, SouthWest, West, NorthWest, Center
};


template <class A>
class ShiftDimensionExpr : public ExprSparseG<ShiftDimensionExpr<A> >{
friend VectorSparseG;
   
   
    const A& a_;
    Compass s_; 
  
    int level;
    bool hier;

public: 
    ShiftDimensionExpr(const A& a, Compass s, bool h):a_(a), s_(s),hier(h){}
    ShiftDimensionExpr(const A& a, Compass s, int l, bool h):a_(a),s_(s),level(l), hier(h){}
 
 
    inline double getValue(int i, const IndexDimension& K) const;
 
    
    ExpressionDescription getDescription() const { return ExpressionDescription(true); }
    AdaptiveSparseGrid* getSparseGrid() const { return  a_.getSparseGrid(); }
};


template<typename A> double ShiftDimensionExpr<A>::getValue(int i, const IndexDimension& K) const
{
    
    
   
    if (hier == true){
      
    AdaptiveSparseGrid* sparseGrid = a_.getSparseGrid();
    //unsigned long endIndex   = sparseGrid->maximalOccupiedSecondTable;
            
            
    dataInteger* secondTable = sparseGrid->getSecondTable();
    dataInteger* primeTable = sparseGrid->getPrimeTable();
            
    if(secondTable[i]!=0) {
        IndexDimension I = K;
        bool nonext = true;
       //if(I.isNotAtBoundary() == false) return a_.getValue(i,K);  
        IndexDimension II;
        while(nonext){
            if (s_ == North)II = I.nextLeft(1);
            if (s_== NorthEast) II = (I.nextLeft(1)).nextRight(0);
            if (s_==East) II = I.nextRight(0);
            if (s_ == SouthEast) II = (I.nextRight(0)).nextRight(1);
            if (s_ == South) II = I.nextRight(1);
            if (s_ == SouthWest) II = (I.nextRight(1)).nextLeft(0);
            if (s_ == West) II = I.nextLeft(0);
            if (s_ == NorthWest) II= (I.nextLeft(0)).nextLeft(1);


            unsigned long indArray =  sparseGrid->hash(II);
        unsigned long k = primeTable[indArray];
        
        if(k == 0) {
            return 0.0;
            //return a_.getValue(i,K);
            
        } else {
            // anaylse what is going on 
            k = k-1; // shift since data are stored with shift
            if(II == sparseGrid->getIndexOfTable(k)) {
                return a_.getValue(k,II);
                
            } else {
                // search in second table;
                unsigned long iNext = secondTable[k];
                    while(iNext > 1) {
                        k = iNext - 2;
                        if(II == sparseGrid->getIndexOfTable(k)) {
                            return a_.getValue(k,II) ;
                            nonext = false;
                        }
                    iNext = secondTable[k];
                    }
                }
            }
        }
    }
      } else {
      
          AdaptiveSparseGrid* sparseGrid = a_.getSparseGrid();
            //unsigned long endIndex   = sparseGrid->maximalOccupiedSecondTable;
                        
            dataInteger* secondTable = sparseGrid->getSecondTable();
            dataInteger* primeTable = sparseGrid->getPrimeTable();
            

       
            
            //assert( level_ <= maxdepth);
            int maxdepth = level;

            if(secondTable[i]!=0) {
                IndexDimension I =K;
                bool nonext = true;
                // todo: andere If Bedingung
                /*if(I.isNotAtBoundary() == false) {
                    return 0;  
                    //return a_.getDataTableVector(i);  
                }*/
                IndexDimension II;
                while(nonext){
                    
                    if (s_ == North)II = I.nextLeft(1,maxdepth);
                    if (s_== NorthEast) II = (I.nextLeft(1,maxdepth)).nextRight(0,maxdepth);
                    if (s_==East) II = I.nextRight(0,maxdepth);
                    if (s_ == SouthEast) II = (I.nextRight(0,maxdepth)).nextRight(1,maxdepth);
                    if (s_ == South) II = I.nextRight(1,maxdepth);
                    if (s_ == SouthWest) II = (I.nextRight(1,maxdepth)).nextLeft(0, maxdepth);
                    if (s_ == West){
                    
                        II = I.nextLeft(0, maxdepth);
             
                    }
                    if (s_ == NorthWest) II= (I.nextLeft(0, maxdepth)).nextLeft(1, maxdepth);


                    unsigned long indArray =  sparseGrid->hash(II);
                    unsigned long k = primeTable[indArray];
                           
                    if(k == 0) {
                       
                      return  0;
                        /*maxdepth= maxdepth-1;
                        if (maxdepth == 0){
                            return a_.getValue(i,K);  
                        }*/
               
                    } else {
                    // anaylse what is going on 
                    k = k-1; // shift since data are stored with shift
                        if(II == sparseGrid->getIndexOfTable(k)) {
                               
                            return a_.getValue(k,II);
                        } else { // search in second table;
                            unsigned long iNext = secondTable[k];
                            while(iNext > 1) {
                            k = iNext - 2;
                                if(II == sparseGrid->getIndexOfTable(k)) {
                                return a_.getValue(k,II) ;
                                nonext = false;
                                }
                            iNext = secondTable[k];
                            }
                        }
                    }
                }
            }
    }
    return 0.0;
}

    
    


////////////////////////////////////////////////////////////////////////
// 4.) class ShiftOperatorDimension
/////////////////////////////////////////



/* @{ */  
/*class ShiftOperatorDimension{
   
public:
*//**
 * 
 * @param s: shift hierarchical in "compass"- direction
 * 
 * **//*
    ShiftOperatorDimension(Compass s):s_(s){
        hier = true;
        myAssert(DimensionSparseGrid == 2);
    }
  
    
    ShiftOperatorDimension(Compass s, int l):s_(s), level(l){
        hier = false;
        myAssert(DimensionSparseGrid == 2);
    }
       
      


    template <class A>
    ShiftDimensionExpr<A> operator()(ExprSparseG<A>& a){
      if (hier == true){
           return ShiftDimensionExpr<A> (a,s_, hier);
      } else {
          return ShiftDimensionExpr<A>(a,s_,level,hier);
      }


        
    }
    
private:
  
    Compass s_;
    int level;
    bool hier;
    
  };*/
//@}


#endif
