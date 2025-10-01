#include "iterator_neu.h"

RectangularIteratorCells::RectangularIteratorCells(CellDimension &A_, CellDimension &B_)
        : A(A_), B(B_) {
    // mit P wird iteriert
    P = A;
    weiter = true;



    for (int d = 0; d < DimensionSparseGrid; ++d) {
        level[d] = A.getDepth(d);
    }

}


void RectangularIteratorCells::operator++() {

    for(int d=0;d<DimensionSparseGrid;++d) {
        if (P.getIndex(d) == B.getIndex(d)) {
            P.replace(d, A.getIndex(d));
        } else {
            CellOneD cellOneD = CellOneD(P.getIndex(d)).nextRightLevel(level[d]);
            P.replace(d, cellOneD);

            return;
        }
    }
    weiter = false;
}


void CellIndexIterator::operator++() {
    for(int d=0;d<DimensionSparseGrid;++d) {
        if (IndexOneD(P.getIndex(d)) == IndexOneD(B.getIndex(d))) {
            P.replace(d, A.getIndex(d));
            cellIndexDirection.setDirection(d,Left);
        } else {
            IndexOneD indexOneD = IndexOneD(P.getIndex(d)).nextRight(level[d]);
            P.replace(d, indexOneD);
            cellIndexDirection.setDirection(d,Right);

            return;
        }
    }
    weiter = false;

}


