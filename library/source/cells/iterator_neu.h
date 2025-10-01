#ifndef ITERATORNEU_H
#define ITERATORNEU_H

#include "celldimension.h"
#include "../myMath/myMath.h"

template<unsigned int N>
struct PowerOfTwo {
    static const unsigned int value = 2 * PowerOfTwo<N-1>::value;
};

template<>
struct PowerOfTwo<0> {
    static const unsigned int value = 1;
};


template <unsigned int N>
struct TriangularNumber {
    static const unsigned int value = N * (N + 1) / 2;
};


class RectangularIteratorCells {
public:
    RectangularIteratorCells(CellDimension& A_, CellDimension &B_);




    CellDimension getCell() { return P; }

    bool goon() { return weiter; }

    void operator++();



protected:
    CellDimension A;
    CellDimension B;
    int level[DimensionSparseGrid];
    CellDimension P;
    bool weiter;

};


class CellIndexDirection{
public:

    CellIndexDirection(){
        for(int d=0; d<DimensionSparseGrid; d++) directions[d]=Left;
    };
    void setDirection(int d, Direction dir){directions[d]=dir;}
    Direction getDirection(int d){return directions[d];}
private:
    Direction directions[DimensionSparseGrid];
};



class CellIndexIterator{
public:
    CellIndexIterator(CellDimension* parent_){

        IndexDimension centerOfCell;
        for(int d=0; d< DimensionSparseGrid; d++){
            level[d]= parent_->getDepth(d);
            centerOfCell.replace(d, parent_->getIndex(d));
        }

        A = centerOfCell;
        B = centerOfCell;
        for(int d=0; d<DimensionSparseGrid; d++){
            A = A.nextLeft(d,level[d]);
            B = B.nextRight(d,level[d]);
            level[d]--;
            cellIndexDirection.setDirection(d,Left);
        }

        //CellIndexDirection cellIndexDirection_cp = cellIndexDirection;

        //P=A;
        //weiter = true;

        //cellIndexDirection1 = new CellIndexDirection[number_indices];
        //Indices = new IndexDimension[number_indices];
        /*for(int i=0; i< number_indices;i++){
            Indices[i]=P;
            cellIndexDirection1[i] = cellIndexDirection;
            operator++();
        }*/

        //cellIndexDirection = cellIndexDirection_cp;
        P=A;
        weiter=true;


    };


    CellIndexIterator(CellIndexIterator const &iterator_parent){
        A = iterator_parent.A;
        B = iterator_parent.B;
        for(int d=0; d<DimensionSparseGrid; d++){
            level[d]=iterator_parent.level[d];
            cellIndexDirection.setDirection(d,iterator_parent.getCellIndexDirection().getDirection(d));
        }

        P=iterator_parent.P;
        weiter = iterator_parent.weiter;


    }

/*    CellIndexIterator(CellIndexIterator  &iterator_parent, int iter_index){
        A = iterator_parent.A;
        B = iterator_parent.B;
        for(int d=0; d<DimensionSparseGrid; d++){
            level[d]=iterator_parent.level[d];
        }


        cellIndexDirection = iterator_parent.getCellIndexDirection(iter_index);
        P=iterator_parent.getIndex(iter_index);

        weiter = true;


    }*/
    CellIndexIterator(CellDimension* parent_, IndexDimension& indexDimension, Direction* dir){

        IndexDimension centerOfCell;
        for(int d=0; d< DimensionSparseGrid; d++){
            level[d]= parent_->getDepth(d);
            centerOfCell.replace(d, parent_->getIndex(d));
        }

        A = centerOfCell;
        B = centerOfCell;
        for(int d=0; d<DimensionSparseGrid; d++){
            A = A.nextLeft(d,level[d]);
            B = B.nextRight(d,level[d]);
            level[d]--;


        }
        P=indexDimension;
        weiter = true;



    };

    IndexDimension getIndex() {
        //return Indices[current_i];
        return P;
    }

    //IndexDimension getIndex(int i){return Indices[i];}

    IndexDimension getIndexByBinary(int i){


        bool binary[DimensionSparseGrid];

        for (int d = DimensionSparseGrid - 1; d >= 0; --d) {
            binary[d] = i & 1;
            i >>= 1;
        }

        IndexDimension Index=A;
        for(int d=0; d<DimensionSparseGrid; d++){
            if(binary[d]){
                Index = Index.nextRight(d,level[d]);
                cellIndexDirection.setDirection(d,Right);
            }else{
                cellIndexDirection.setDirection(d,Left);
            }
        }
        return Index;
    }






    CellIndexDirection getCellIndexDirection() const {return cellIndexDirection;}
    //CellIndexDirection getCellIndexDirection(int i) const {return cellIndexDirection1[i];}

    bool goon() {
               return weiter; }

    void operator++();




private:


    IndexDimension A;
    IndexDimension B;
    IndexDimension P;

    int number_indices = PowerOfTwo<DimensionSparseGrid>::value;

    bool weiter;
    int level[DimensionSparseGrid];
    //Direction directions[DimensionSparseGrid];
    CellIndexDirection cellIndexDirection;
    //CellIndexDirection cellIndexDirection1[PowerOfTwo<DimensionSparseGrid>::value];
    //IndexDimension Indices[PowerOfTwo<DimensionSparseGrid>::value];

};



#endif