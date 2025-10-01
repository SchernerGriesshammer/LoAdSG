//
// Created by to35jepo on 4/11/23.
//

#include "RectangularIterator.h"
#include "../sgrid/depth.h"

RectangularIterator::RectangularIterator(IndexDimension &A_, IndexDimension &B_)
        : A(A_), B(B_) {
    // mit P wird iteriert
    P = A;
    weiter = true;

    Depth T = Depth(A);

    for (int d = 0; d < DimensionSparseGrid; ++d) {
        level[d] = T.at(d);
    }

}

RectangularIterator::RectangularIterator(IndexDimension &A_, IndexDimension &B_, Depth T) : A(A_), B(B_) {

    // mit P wird iteriert
    P = A;
    weiter = true;
    for (int d = 0; d < DimensionSparseGrid; ++d) {

        level[d] = T.at(d);

    }

}
void RectangularIterator::operator++() {

    for(int d=0;d<DimensionSparseGrid;++d) {
        if (IndexOneD(P.getIndex(d)) == IndexOneD(B.getIndex(d))) {
            P.replace(d, A.getIndex(d));
        } else {
            IndexOneD indexOneD = IndexOneD(P.getIndex(d)).nextRight(level[d]);
            P.replace(d, indexOneD);

            return;
        }
    }
    weiter = false;
}

void RectangularIterator::nextnext() {
    for(int d=0;d<DimensionSparseGrid;++d) {
        if (IndexOneD(P.getIndex(d)) == IndexOneD(B.getIndex(d))) {
            P.replace(d, A.getIndex(d));
        } else {
            IndexOneD indexOneD = IndexOneD(P.getIndex(d)).nextRight(level[d]);
            //P.replace(d, indexOneD);
            if(level[d]>1)
            P = P.nextRight(d,level[d]).nextRight(d,level[d]);
            else
            P = P.nextRight(d,level[d]);
            return;
        }
    }
    weiter = false;
}

void RectangularIterator::operator2(int dim) {

    int dd = dim;
    int d;


    do {
        d = dd % DimensionSparseGrid;
        if (IndexOneD(P.getIndex(d)) == IndexOneD(B.getIndex(d))) {
            P.replace(d, A.getIndex(d));
        } else {
            IndexOneD indexOneD = IndexOneD(P.getIndex(d)).nextRight(level[d]);
            P.replace(d, indexOneD);
            return;
        }
        dd++;

    } while (dd % (DimensionSparseGrid) != (dim % DimensionSparseGrid));
    weiter = false;
}


RectangularIteratorExact::RectangularIteratorExact(IndexDimension &A_, IndexDimension &B_, Depth T) : RectangularIterator(A_,
                                                                                                                 B_) {
    for(int d=0; d<DimensionSparseGrid; d++)
        level[d]=T.at(d);

}

void RectangularIteratorExact::operator++() {

    for(int d=0;d<DimensionSparseGrid;++d) {
        if (IndexOneD(P.getIndex(d)) == IndexOneD(B.getIndex(d))) {
            P.replace(d, A.getIndex(d));
        } else {
            IndexOneD indexOneD = IndexOneD(P.getIndex(d)).nextRight(level[d]);
            P.replace(d, indexOneD);
            return;
        }
    }
    weiter = false;

}

