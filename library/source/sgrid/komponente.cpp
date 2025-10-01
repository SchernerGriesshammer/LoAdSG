//
// Created by scherner on 29.04.21.
//

#include "komponente.h"



unsigned long ZusammenhangsKomponente::findAnicht(unsigned long iNow, unsigned long &iNow_mapping) {
    SingleDepthHashGrid &depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(T);
    const auto &mapping = depthGrid._mapPosToGridPos;


    size_t i=iNow_mapping+1;


    for (; i < mapping.size(); i++) {
        unsigned long k=mapping[i];
        if (marker[k] == nicht) {
            if (grid->workonindex(k)) {
                return k;
                iNow_mapping=i;
            }
        }
    }



    unsigned long iMax = grid->getMaximalOccupiedSecondTable();
    /* for (unsigned long i = iNow + 1; i < iMax; ++i) {
         IndexDimension Index = grid->getIndexOfTable(i);
         if (marker[i] == nicht && Depth(Index) == T) {
             if (grid->workonindex(i))
                 return i;
         }
     }*/

    return iMax; // also nichts gefunden
}

void ZusammenhangsKomponente::recursiveMarkArbeite(IndexDimension indexNow, unsigned long iNow) {
    marker[iNow] = arbeite;


    //Rekursion
    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
        IndexDimension J = indexNow.nextFive(&mc, T);

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
}

ZusammenhangsKomponente::ZusammenhangsKomponente(AdaptiveSparseGrid_Base *grid_, Depth T_) : T(T_) {

    grid = grid_;
    secondTableLength = grid->getMaximalOccupiedSecondTable();
    marker = new Marker[secondTableLength];
    for (unsigned long i = 0; i < secondTableLength; ++i)
        marker[i] = nicht;
};