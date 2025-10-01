//
// Created by scherner on 01.03.21.
//

#ifndef GRUN_KOMPONENTE_H
#define GRUN_KOMPONENTE_H
#include "sparseGrid.h"




class ZusammenhangsKomponente {
public:
    ZusammenhangsKomponente(AdaptiveSparseGrid_Base *grid_, Depth T_);

    ~ZusammenhangsKomponente() { delete[] marker; };
    enum Marker {
        nicht, arbeite
    };

    /**
     * @return ein index mit Markierung
     **/
    unsigned long findAnicht(unsigned long iNow, unsigned long& iNow_mapping);

    /**
     * setzt Maximum und Minimum der Komponente
     * @param iNow mit einer nicht Markierung
     * @return zugehöriger Index
     **/
    IndexDimension StartSearchComponent(unsigned long iNow) {

        maxIndex = grid->getIndexOfTable(iNow);
        minIndex = maxIndex;
        T = Depth(maxIndex);
        return minIndex;

    };

    IndexDimension getMinIndex() { return minIndex; };

    IndexDimension getMaxIndex() { return maxIndex; };

    Depth getDepth() { return T; };


    virtual void recursiveMarkArbeite(IndexDimension indexNow, unsigned long iNow);

protected:
    AdaptiveSparseGrid_Base *grid;
    Marker *marker;
    unsigned long secondTableLength;

    IndexDimension maxIndex;
    IndexDimension minIndex;

    Depth T;
};


class ZusammenhangsKomponente_Neumann : public ZusammenhangsKomponente {
public:
    ZusammenhangsKomponente_Neumann(AdaptiveSparseGrid_Base *grid_, Depth T_) : ZusammenhangsKomponente(grid_, T_) {};

    inline void recursiveMarkArbeite(IndexDimension indexNow, unsigned long iNow) {
        marker[iNow] = arbeite;


        //Rekursion
        for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
            double val;
            IndexDimension J = indexNow.nextFive_Neumann(&mc, T, &val);

            unsigned long indexOfData;
            if (grid->occupied(indexOfData, J)) {
                if (grid->workonindex(indexOfData)) {

                    Marker m = marker[indexOfData];
                    if (m != arbeite) {
                        if (Depth(J) == T) {
                            maxIndex = IndexDimension::Maximum(maxIndex, J);
                            minIndex = IndexDimension::Minimum(minIndex, J);
                            recursiveMarkArbeite(J, indexOfData);
                        } else {

                            maxIndex = IndexDimension::Maximum(maxIndex, J);
                            minIndex = IndexDimension::Minimum(minIndex, J);
                        }
                    }
                }
            }
        }
        //max min erneuern -> wichtig: Tiefe ändert sich dadurch nicht
        maxIndex = IndexDimension::Maximum(maxIndex, indexNow);
        minIndex = IndexDimension::Minimum(minIndex, indexNow);
    };

    inline void markArbeite(unsigned long iNow) {
        marker[iNow] = arbeite;
    }
};


#endif //GRUN_KOMPONENTE_H
