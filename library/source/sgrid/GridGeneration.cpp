//
// Created by to35jepo on 4/12/23.
//

#include <omp.h>
#include "GridGeneration.h"
#include "../iterator/RectangularIterator.h"
#include "../iterator/depthIterator.h"
#include "../mympi.h"


void support(const IndexDimension& P, IndexDimension& Imin, IndexDimension& Imax){
    Depth TP(P);
    Imin = P;
    Imax = P;
    for(int d=0; d<DimensionSparseGrid; d++){

        int t= TP.at(d);

        if(t<=2){
            Imin.replace(d,0);
            Imax.replace(d,1);
        }else {
            if (P.isLinksRandNah(d)){
                Imin.replace(d,0);
            }else{
                Imin = Imin.nextLeft(d).nextLeft(d,t).nextLeft(d,t);
            }
            if (P.isRechtsRandNah(d)){
                Imax.replace(d,1);
            }else{
                Imax = Imax.nextRight(d).nextRight(d,t).nextRight(d,t);
            }
        }
    }
}
bool overlap2(const IndexDimension& P, const IndexDimension& Q,IndexDimension& Imin, IndexDimension& Imax){

    Depth TP(P);
    Depth TQ(Q);

    IndexDimension IminP = P;
    IndexDimension ImaxP = P;

    IndexDimension IminQ = Q;
    IndexDimension ImaxQ = Q;

    for(int d=0; d<DimensionSparseGrid; d++){

        int t= TP.at(d);

        if(t<=2){
            IminP.replace(d,0);
            ImaxP.replace(d,1);

        }else {
            if (P.isLinksRandNah(d)){
                IminP.replace(d,0);
            }else{
                IminP = IminP.nextLeft(d).nextLeft(d,t).nextLeft(d,t);
            }
            if (P.isRechtsRandNah(d)){
                ImaxP.replace(d,1);
            }else{
                ImaxP = ImaxP.nextRight(d).nextRight(d,t).nextRight(d,t);
            }
        }



        t= TQ.at(d);

        if(t<=2){
            IminQ.replace(d,0);
            ImaxQ.replace(d,1);

        }else {
            if (Q.isLinksRandNah(d)){
                IminQ.replace(d,0);
            }else{
                IminQ = IminQ.nextLeft(d).nextLeft(d,t).nextLeft(d,t);
            }
            if (Q.isRechtsRandNah(d)){
                ImaxQ.replace(d,1);
            }else{
                ImaxQ = ImaxQ.nextRight(d).nextRight(d,t).nextRight(d,t);
            }
        }

        if(ImaxQ.coordinate(d)<IminP.coordinate(d) || ImaxP.coordinate(d)<IminQ.coordinate(d)) {
            return false;
        }


        if(IminP.coordinate(d)>IminQ.coordinate(d))
            Imin.replace(d,IminP.getIndex(d));
        else
            Imin.replace(d,IminQ.getIndex(d));

        if(ImaxP.coordinate(d)<ImaxQ.coordinate(d))
            Imax.replace(d,ImaxP.getIndex(d));
        else
            Imax.replace(d,ImaxQ.getIndex(d));
    }

    return true;



}

bool overlap(unsigned long p,unsigned long q, IndexDimension& Imin, IndexDimension& Imax, AdaptiveSparseGrid& grid){
    IndexDimension IminP,IminQ,ImaxP,ImaxQ;

    // ToDo: Speichere IminP und ImaxP für jeden Punkt einmal!

    IminP = grid.getSupportIMin(p);
    ImaxP = grid.getSupportIMax(p);

    IminQ = grid.getSupportIMin(q);
    ImaxQ = grid.getSupportIMax(q);






    for(int d=0; d<DimensionSparseGrid; d++){
        if(ImaxQ.coordinate(d)<IminP.coordinate(d) || ImaxP.coordinate(d)<IminQ.coordinate(d)) {
            return false;
        }


        if(IminP.coordinate(d)>IminQ.coordinate(d))
            Imin.replace(d,IminP.getIndex(d));
        else
            Imin.replace(d,IminQ.getIndex(d));

        if(ImaxP.coordinate(d)<ImaxQ.coordinate(d))
            Imax.replace(d,ImaxP.getIndex(d));
        else
            Imax.replace(d,ImaxQ.getIndex(d));
    }

    return true;

}

Depth max(Depth& TP, Depth& TQ){
    Depth T=TQ;
    for(int d=0;d<DimensionSparseGrid;d++){
        if(TP.at(d)>TQ.at(d))
            T.set(TP.at(d),d);
    }
    return T;
}

bool refill(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q){

    bool returnbool=false;

        IndexDimension Imin, Imax;



        if (overlap(p, q, Imin, Imax, grid)) {


            bool addNodes = false;


            for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
                IndexDimension I = iter.getIndex();
                Depth Tneu(I);
                unsigned long k;
                if (Tneu == T) {
                    if (grid.occupied(k, I)) {
                            addNodes = true;
                            break;
                    }
                }
            }

            if (addNodes) {
                for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
                    IndexDimension I = iter.getIndex();
                    if(I.isNotAtBoundary()){
                        Depth Tneu(I);
                        unsigned long k;
                        if (Tneu == T) {
                            if (!grid.occupied(k, I)) {
                                grid.AddPoint(I, false);
                                returnbool = true;
                            }
                        }
                    }
                }
            }
        }


    return returnbool;
}

bool refill_NEU(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q, std::vector<IndexDimension>& pointsToAdd){

    bool returnbool=false;

    IndexDimension Imin, Imax;



    if (overlap(p, q, Imin, Imax, grid)) {


        bool addNodes = false;



        for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
            IndexDimension I = iter.getIndex();
            Depth Tneu(I);
            unsigned long k;
            if (Tneu == T) {
                if (grid.occupied(k, I)) {
                    addNodes = true;
                    break;
                }
            }
        }

        if (addNodes) {

            for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
                IndexDimension I = iter.getIndex();
                if(I.isNotAtBoundary()){
                    Depth Tneu(I);
                    unsigned long k;
                    if (Tneu == T) {
                        if (!grid.occupied(k, I)) {

                            //grid.AddPoint(I, false);

                            pointsToAdd.push_back(I);
                            returnbool = true;
                        }
                    }
                }
            }
        }
    }


    return returnbool;
}



bool refill2(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T){

    bool returnbool = false;


    IndexDimension Imin, Imax;





   if (overlap2(P, Q, Imin, Imax)) {


            bool addNodes = false;
            bool FullNode = false;


/*
           for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
                IndexDimension I = iter.getIndex();
                Depth Tneu(I);
                unsigned long k;
                if (Tneu == T) {
                    bool occu = grid.occupied(k,I);
                    if(!occu && FullNode){
                        addNodes = true;
                        break;
                    }
                    if(occu && !FullNode){
                        FullNode = true;
                    }
                }
            }*/

            if (addNodes) {
                for(RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter){
                    IndexDimension I = iter.getIndex();
                    if(I.isNotAtBoundary())
                        if(grid.AddPoint(I, false)){
                            returnbool = true;
                        }
                }
            }
        }




    return returnbool;
}

void addPoints(AdaptiveSparseGrid& dGrid){


    for (bool changed = true; changed;) {
        changed = false;
        for (unsigned long k = 0; k < dGrid.getMaximalOccupiedSecondTable(); k++) {
            if (dGrid.getActiveTable()[k]) {

                for (unsigned long j = k + 1; j < dGrid.getMaximalOccupiedSecondTable(); j++) {
                    if (dGrid.getActiveTable()[j]) {
                        IndexDimension P = dGrid.getIndexOfTable(k);
                        IndexDimension Q = dGrid.getIndexOfTable(j);
/*                        if (refill(dGrid, P, Q)) {
                            changed = true;
                        }*/
                    }
                }
            }
        }
    }
}

void addPoints2(AdaptiveSparseGrid& dGrid){

    DepthList list(dGrid);
    for (bool changed = true; changed;) {
        changed = false;

        for (auto it = list.begin_all(); it != list.end_all(); ++it) {
            Depth T1 = *it;
            auto it2 = it;
            ++it2;
            for (; it2 != list.end_all(); ++it2) {
                Depth T2 = *it2;


                int lone = 0;
                for (int d = 0; d < DimensionSparseGrid; d++) {
                    if (T1.at(d) != T2.at(d)) {
                        lone++;
                    }
                }

                if (lone >= 2) {

                    Depth T = max(T1,T2);




                    SingleDepthHashGrid &depthGrid1 = dGrid.getMultiDepthHashGrid()->getGridForDepth(T1);
                    SingleDepthHashGrid &depthGrid2 = dGrid.getMultiDepthHashGrid()->getGridForDepth(T2);

                    const auto &mapping1 = depthGrid1._mapPosToGridPos;
                    const auto &mapping2 = depthGrid2._mapPosToGridPos;
                    for (size_t i = 0; i < mapping1.size(); i++) {
                        IndexDimension I1 = dGrid.getIndexOfTable(mapping1[i]);

                        if (dGrid.getActiveTable()[mapping1[i]]) {
                            for (size_t j = 0; j < mapping2.size(); j++) {
                                if (dGrid.getActiveTable()[mapping2[j]]) {

                                    IndexDimension I2 = dGrid.getIndexOfTable(mapping2[j]);


                                   if (refill(dGrid, I1, I2,T,mapping1[i],mapping2[j])) {

                                        changed = true;

                                    }
                                }
                            }
                        }
                    }
                }

            }

        }
    }



/*
    for (bool changed = true; changed;) {
        changed = false;
        for (unsigned long k = 0; k < dGrid.getMaximalOccupiedSecondTable(); k++) {
            if (dGrid.getActiveTable()[k]) {

                for (unsigned long j = k + 1; j < dGrid.getMaximalOccupiedSecondTable(); j++) {
                    if (dGrid.getActiveTable()[j]) {
                        IndexDimension P = dGrid.getIndexOfTable(k);
                        IndexDimension Q = dGrid.getIndexOfTable(j);
                        if (refill(dGrid, P, Q)) {
                            changed = true;
                        }
                    }
                }
            }
        }
    }*/
}

void addPoints2_NEU(AdaptiveSparseGrid& dGrid, DepthList list) {

    int rank=0; int size=1;
#ifdef MY_MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

#endif
    //DepthList list(dGrid);

    std::vector<Depth> dlist;
    for (auto it = list.begin_all(); it != list.end_all(); ++it) {
        Depth T1 = *it;
        dlist.push_back(T1);
    }

    for (bool changed = true; changed;) {
        changed = false;

        std::vector<IndexDimension> pointsToAdd;

        for (int i1 = 0; i1 < dlist.size(); i1++) {
                if (i1 % dlist.size() == rank) {
                    Depth T1 = dlist.at(i1);
                    auto i2 = i1;
                    ++i2;
                    for (; i2 < dlist.size(); i2++) {
                        Depth T2 = dlist.at(i2);


                        int lone = 0;
                        for (int d = 0; d < DimensionSparseGrid; d++) {
                            if (T1.at(d) != T2.at(d)) {
                                lone++;
                            }
                        }

                        if (lone >= 2) {

                            Depth T = max(T1, T2);


                            SingleDepthHashGrid &depthGrid1 = dGrid.getMultiDepthHashGrid()->getGridForDepth(T1);
                            SingleDepthHashGrid &depthGrid2 = dGrid.getMultiDepthHashGrid()->getGridForDepth(T2);

                            const auto &mapping1 = depthGrid1._mapPosToGridPos;
                            const auto &mapping2 = depthGrid2._mapPosToGridPos;

                            #pragma omp parallel for schedule(dynamic)
                            for (size_t i = 0; i < mapping1.size(); i++) {
                                IndexDimension I1 = dGrid.getIndexOfTable(mapping1[i]);

                                if (dGrid.getActiveTable()[mapping1[i]]) {
                                    for (size_t j = 0; j < mapping2.size(); j++) {
                                        if (dGrid.getActiveTable()[mapping2[j]]) {


                                            IndexDimension I2 = dGrid.getIndexOfTable(mapping2[j]);


                                            if (refill_NEU(dGrid, I1, I2, T, mapping1[i], mapping2[j], pointsToAdd)) {
                                                changed = true;

                                            }


                                        }
                                    }
                                }
                            }

                        }

                    }


                    for (IndexDimension I: pointsToAdd) {
                        dGrid.AddPoint(I, false);
                    }

#ifdef MY_MPI_ON
                    int numberAddPoints[size];
                    int numberAddPoints_send = pointsToAdd.size();
                    int numberAddPoints_rec = 0;
                    MPI_Request request;


                    for (int i = 0; i < size; i++) {
                        if (i != rank) {

                            MPI_Isend(&numberAddPoints_send, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
                            MPI_Irecv(&numberAddPoints_rec, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
                            numberAddPoints[i] = numberAddPoints_rec;
                            if (numberAddPoints_rec > 0) changed = true;
                        }
                    }
                    for (int i = 0; i < size; i++) {
                        if (i != rank) {
                            int index[DimensionSparseGrid];
                            if(pointsToAdd.size()>0){
                            for (IndexDimension I: pointsToAdd) {
                                for (int d = 0; d < DimensionSparseGrid; d++) {
                                    index[d] = I.getIndex(d);
                                }
                                MPI_Isend(&index, DimensionSparseGrid, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
                            }
                            }
                            if(numberAddPoints[i]>0){
                            IndexDimension I;
                            for (int k = 0; k < numberAddPoints[i]; k++) {

                                MPI_Irecv(&index, DimensionSparseGrid, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
                                for (int d = 0; d < DimensionSparseGrid; d++) {
                                    I.replace(d, index[d]);
                                }
                                dGrid.AddPoint(I, false);
                            }
                            }
                        }
                    }
                   MPI_Barrier(MPI_COMM_WORLD);
#endif
                }
        }

    }

}



void addPoints(AdaptiveSparseGrid& dGrid, IndexDimension P){


                for (unsigned long j = 0; j < dGrid.getMaximalOccupiedSecondTable(); j++) {
                    if (dGrid.getActiveTable()[j]) {
                        IndexDimension Q = dGrid.getIndexOfTable(j);
            /*            if(refill(dGrid, P, Q)) {

                        }*/
                    }
                }

}