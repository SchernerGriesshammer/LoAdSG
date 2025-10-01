//
// Created by to35jepo on 4/11/23.
//
#include "../indices/index.h"
#ifndef SGRUN_RECTANGULARITERATOR_H
#define SGRUN_RECTANGULARITERATOR_H

/**
 * iteriert alle Gitterpunkte von A bis B
 * Annahme A und B haben genau die gleiche Tiefe
 * Die Iteration erfolgt dann  2*T(A) oder 2*T(B)
 **/
class RectangularIterator {
public:
    RectangularIterator(IndexDimension &A_, IndexDimension &B_);

    /**
     *
     * @param A_
     * @param B_
     * @param T  BSP: Rectangle soll von A nach B, diese haben aber unterschiedliche Tiefen. Mit T kann man angeben,
     * bezüglich welcher Tiefen das Rechteck durchlaufen werden soll.
     * wird eventuell noch gebraucht?
     */
    RectangularIterator(IndexDimension &A_, IndexDimension &B_, Depth T);

    IndexDimension getIndex() { return P; }

    bool goon() { return weiter; }

    void begin(){
        P=A;
    weiter = true;}

    void operator++();

    void nextnext();

    /**
     * Neuer Operator muss noch getestet und umbenannt werden. Verändert die Reihenfolge.
     * @param dim
     */
    void operator2(int dim);

protected:
    IndexDimension A;
    IndexDimension B;
    int level[DimensionSparseGrid];
    IndexDimension P;
    bool weiter;

};

class RectangularIteratorExact: public RectangularIterator{


public:
    RectangularIteratorExact(IndexDimension &A_, IndexDimension &B_, Depth T);
    void operator++();



};
#endif //SGRUN_RECTANGULARITERATOR_H
