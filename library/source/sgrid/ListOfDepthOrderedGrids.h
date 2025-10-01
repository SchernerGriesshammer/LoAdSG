/****************************
 * Test for ordered grids
 ***************************/

#ifndef LISTOFORDEREDGRIDS_H
#define LISTOFORDEREDGRIDS_H

#include <list>
#include <vector>
#include <iostream>


/**
 * Objekt, das alle Gitter Punkte mit Tiefe T enthält.
 *
 * Beispiel Anwendung:
 * @code
 *
    ListOfDepthOrderedSubgrids orderedSubgrids(grid);
    ListOfDepthOrderedSubgrids::iterator outerIterationsSubgrids(orderedSubgrids);
    outerIterationsSubgrids.gotoBegin();

        //Gehe alle Tiefen von Coarse to Fine durch
    do {
        Depth T = outerIterationsSubgrids.getSubgrid()->getT();
        SubgridFixedDepth* subgrid = outerIterationsSubgrids.getSubgrid();

    }while(outerIterationsSubgrids.next());
    @endcode
 **/
class SubgridFixedDepth {
private:
    class IndexPlusi {
    public:
        IndexPlusi(IndexDimension I_, unsigned long i_) : I(I_), i(i_) {}

        IndexDimension I;
        unsigned long i;

        bool operator==(IndexPlusi vergleich) {
            if (I == vergleich.I) return true;
            if (i == vergleich.i) return true;
            return false;

        }

        void print() {
            I.Print();
        }
    };

public:
    SubgridFixedDepth(Depth &T_);


    /**
     * @param p muss hier Tiefe 
     * @param T 
     * haben
     **/
    SubgridFixedDepth(Depth &T_, IndexDimension p, unsigned long i);


    Depth getT() { return T; }

    /**
     * @param p ein neuer Punkt, der noch nicht im SubgridFixedDepth enthalten ist. Das wird geprüft.
     **/
    bool addPoint(IndexDimension &p, unsigned long i);

    bool isOccupied(IndexDimension &p, unsigned long i);


    /**
     * @code
     *  SubgridFixedDepth::iterator iter(subgrid);
        do{
            IndexDimension Index = iter.getPoint();
        }while(iter.next());
        @endcode
    */
    class iterator {
    public:
        iterator(SubgridFixedDepth &dataStructure_);

        void gotobegin();

        bool next();

        bool empty();

        IndexDimension getPoint() {

            return (*it).I;
        }

        unsigned long geti() {

            return it->i;
        }

    private:
        SubgridFixedDepth *dataStructure;
        std::list<SubgridFixedDepth::IndexPlusi>::iterator it;
    };

private:
    Depth T;
    std::list<SubgridFixedDepth::IndexPlusi> listPoints;
};

/**
 * Liste von Gittern mit unterschiedlicher Tiefe T
 **/
class ListOfDepthOrderedSubgrids {
public:
    ListOfDepthOrderedSubgrids();

    ListOfDepthOrderedSubgrids(AdaptiveSparseGrid_Base &grid);

    ListOfDepthOrderedSubgrids(AdaptiveSparseGrid_Base &grid, bool *restrictions);


    void addPoint(IndexDimension &p, unsigned long i);


    class iterator {
    public:
        iterator(ListOfDepthOrderedSubgrids &dataStructure_);

        void gotoBegin();

        void gotoEnd();

        bool previous();

        bool next();

        bool empty();


        int get_t() { return t; };

        Depth getDepth() {
            return (*it).getT();
        }

        SubgridFixedDepth *getSubgrid() { return &(*it); }

    private:
        ListOfDepthOrderedSubgrids *dataStructure;
        std::list<SubgridFixedDepth>::iterator it;
        int t;
    };


    bool isIncluded(Depth &T);

    inline void sortDepths(bool *Restrictions);

    inline void SortiereTiefenBoundary(bool *Restrictions);


    SubgridFixedDepth *getSubgrid(Depth T) {
        iterator iter(*this);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getSubgrid()->getT();
            if (Tlocal == T) return iter.getSubgrid();
        } while (iter.next());

        cout << "error no subgrid of required depth" << endl;
        return iter.getSubgrid();

    };

    std::list<Depth> *getSortierteTiefen() { return &SortierteTiefen; }

private:
    std::vector<std::list<SubgridFixedDepth>> TiefeVector;
    std::list<Depth> SortierteTiefen;
    int MaximalDepth[DimensionSparseGrid];

    inline void recursiveDepth(bool *Restrictions, bool checkProlongation, int d, Depth *T);

    inline void recursiveDepthBoundary(bool *Restrictions, bool checkProlongation, int d, Depth *T);

    void recursiveDepthCoarse(int d, Depth *T);
};


/////////////////////////////////// inline
void ListOfDepthOrderedSubgrids::sortDepths(bool *Restrictions) {

    SortierteTiefen.clear();
    Depth T(0);
    recursiveDepth(Restrictions, true, DimensionSparseGrid - 1, &T);

}


inline void ListOfDepthOrderedSubgrids::SortiereTiefenBoundary(bool *Restrictions) {

    SortierteTiefen.clear();
    Depth T(0);
    recursiveDepthBoundary(Restrictions, true, DimensionSparseGrid - 1, &T);

}


inline void ListOfDepthOrderedSubgrids::recursiveDepth(bool *Restrictions, bool checkProlongation, int d, Depth *T) {
    if (checkProlongation && d > 0) {
        if (!Restrictions[d]) {
            for (int t = 1; t <= MaximalDepth[d]; t++) {
                T->set(t, d);
                recursiveDepth(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (checkProlongation && d == 0) {
        if (!Restrictions[d]) {
            for (int t = 1; t <= MaximalDepth[d]; t++) {
                T->set(t, d);
                recursiveDepth(Restrictions, false, DimensionSparseGrid - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, false, DimensionSparseGrid - 1, T);
        }
    }

    if ((!checkProlongation) && d >= 0) {
        if (Restrictions[d]) {
            for (int t = MaximalDepth[d]; t > 1; t--) {
                T->set(t, d);
                recursiveDepth(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepth(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (d < 0 && isIncluded(*T)) {
        SortierteTiefen.push_back(*T);
    }

}


inline void
ListOfDepthOrderedSubgrids::recursiveDepthBoundary(bool *Restrictions, bool checkProlongation, int d, Depth *T) {

    if (checkProlongation && d > 0) {
        if (!Restrictions[d]) {
            for (int t = 0; t <= MaximalDepth[d]; t++) {
                T->set(t, d);
                recursiveDepthBoundary(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepthBoundary(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (checkProlongation && d == 0) {
        if (!Restrictions[d]) {
            for (int t = 0; t <= MaximalDepth[d]; t++) {
                T->set(t, d);
                recursiveDepthBoundary(Restrictions, false, DimensionSparseGrid - 1, T);
            }
        } else {
            recursiveDepthBoundary(Restrictions, false, DimensionSparseGrid - 1, T);
        }
    }

    if ((!checkProlongation) && d >= 0) {
        if (Restrictions[d]) {
            for (int t = MaximalDepth[d]; t >= 1; t--) {
                T->set(t, d);
                recursiveDepthBoundary(Restrictions, checkProlongation, d - 1, T);
            }
        } else {
            recursiveDepthBoundary(Restrictions, checkProlongation, d - 1, T);
        }
    }

    if (d < 0 && isIncluded(*T)) {
        SortierteTiefen.push_back(*T);
    }

}

#endif
