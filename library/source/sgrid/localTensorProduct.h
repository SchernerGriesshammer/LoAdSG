//
// Created by scherner on 21.09.21.
//
#include "../iterator/RectangularIterator.h"

#ifndef GRUN_LOCALTENSORPRODUCT_H
#define GRUN_LOCALTENSORPRODUCT_H


/**
 * local class for tensor product construction
 *
 *      unsigned long iNow = constructorTensor.findAnichtWithDepth(-1, T);
        while (iNow < maximalOccupiedSecondTable) {


            IndexDimension indexNow = constructorTensor.StartSearchComponent(iNow);


            // Zusammenhangskomponente wird untersucht und bereits vorhandene Punkte auf arbeite gesetzt
            constructorTensor.recursiveMarkArbeite2(indexNow, iNow);

            constructorTensor.checkRechteck(indexNow);


            constructorTensor.checkStorage();
            iNow = constructorTensor.findAnichtWithDepth(iNow, T);

        }
 **/
class LocalForTensor {
public:
    explicit LocalForTensor(AdaptiveSparseGrid *grid_);

    ~LocalForTensor();

    enum Marker {
        nicht, fertig, arbeite
    };

    /**
     * @return ein index mit Markierung
     **/
    unsigned long findAnicht(unsigned long iNow);

    unsigned long findAnichtWithDepth(unsigned long iNow, Depth &T_);

    /**
     * setzt Maximum und Minimum der Komponente
     * @param iNow mit einer nicht Markierung
     * @return zugehöriger Index
     **/
    IndexDimension StartSearchComponent(unsigned long iNow);

    /**
     * untersucht rekursive die 3^d-1 Nachbarpunkte
     * @param touch eine andere Zusammenhangskomponente wurde vorher berührt
     * @param indexNow index zu iNow, Startpunkt
     * @return true falls ein Punkt Markierung von fertig zu arbeite ändert; das heisst eine andere Zusammenhangskomponente wurde berührt
     **/
    void recursiveMarkArbeite(IndexDimension indexNow, unsigned long iNow);

    void recursiveMarkArbeite2(IndexDimension indexNow, unsigned long iNow);

    bool checkRechteck(IndexDimension indexNow);

    /**
     * markiert alle Punkte zu arbeite und fügt fehlende Punkte hinzu
     * @return true falls eine andere Zusammenhangskomponente wurde berührt
     **/
    bool markArbeiteAndAddRechteck();

    void checkStorage() {} //  muss man noch implementieren !!!!!!
    void markRechteck(Marker mm);

    void PrintIMaxIMin() {
        cout << "print max and min" << endl;
        maxIndex.PrintCoord();
        cout << endl;
        minIndex.PrintCoord();
    }

    IndexDimension getMin(){return minIndex;}
    IndexDimension getMax(){return maxIndex;}

    void PrintDepth() {
        for (int i = 0; i < DimensionSparseGrid; i++) cout << T.at(i) << endl;
    }
/*void AdaptiveSparseGrid::CompleteToLocalTensorProductGrid() {
    // Grid hat immer zuerst hierarchische Struktur und dann noch als TP erweitert
    completeGrid();

   //1. Start
   LocalForTensor constructorTensor(this);

   // Iteration
    for (unsigned long iNow = constructorTensor.findAnicht(-1); iNow < maximalOccupiedSecondTable;
         iNow = constructorTensor.findAnicht(iNow)) {

        IndexDimension indexNow = constructorTensor.StartSearchComponent(iNow);

        // Zusammenhangskomponente wird untersucht und bereits vorhandene Punkte auf arbeite gesetzt
        constructorTensor.recursiveMarkArbeite(indexNow, iNow);


        bool touch = false;
        do {
           // fehlende Rechteckpunkte hinzufügen, markieren arbeite und touch feststellen
           touch = constructorTensor.markArbeiteAndAddRechteck();
           constructorTensor.checkStorage();

           if(touch){
               constructorTensor.markRechteck(LocalForTensor::nicht);
               constructorTensor.recursiveMarkArbeite(indexNow,iNow);
           }

       }while(touch);

       constructorTensor.markRechteck(LocalForTensor::fertig);
       constructorTensor.checkStorage();

   if(TEST)  PrintActiveHanging(3);
   }
   // completeGrid();

}*/
protected:
    AdaptiveSparseGrid *grid;
    Marker *marker;
    unsigned long secondTableLength;

    IndexDimension maxIndex;
    IndexDimension minIndex;

    Depth T;
};

bool checkRechteck(IndexDimension minIndex,IndexDimension maxIndex, Depth T, AdaptiveSparseGrid* grid);
#endif //GRUN_LOCALTENSORPRODUCT_H
