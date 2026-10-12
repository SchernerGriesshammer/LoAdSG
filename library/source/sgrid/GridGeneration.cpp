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

// refill_v2 / refill_NEU_v2: functionally identical to refill()/refill_NEU() (same
// D_(p,q)^full completion criterion), but faster in two ways:
//  1. The "does this box already contain an occupied depth-T point" check and the
//     "add every missing depth-T point" fill are done in a single RectangularIteratorExact
//     pass instead of two. Cells seen before the first occupied hit are buffered in
//     `pending` and only committed if a hit is found later in the same pass; cells seen
//     after the hit are added immediately. This preserves the original semantics (add
//     nothing unless some point in the box already exists) while walking the box once.
//  2. occupied() checks go through the depth-T SingleDepthHashGrid (fetched once, outside
//     the loop) instead of the whole-grid hash table. Every candidate cell here already
//     satisfies Tneu == T, and every point added to the grid is unconditionally registered
//     in its own depth's SingleDepthHashGrid too (see setIndexInTable()), so this is the
//     same lookup restricted to a much smaller, depth-local table.
bool refill_v2(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q){

    bool returnbool = false;

    IndexDimension Imin, Imax;

    if (overlap(p, q, Imin, Imax, grid)) {

        SingleDepthHashGrid& depthGrid = grid.getMultiDepthHashGrid()->getGridForDepth(T);

        bool addNodes = false;
        std::vector<IndexDimension> pending;

        for (RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
            IndexDimension I = iter.getIndex();
            Depth Tneu(I);
            if (!(Tneu == T)) continue;

            unsigned long k;
            if (depthGrid.occupied(k, I)) {
                addNodes = true;
            } else if (I.isNotAtBoundary()) {
                if (addNodes) {
                    grid.AddPoint(I, false);
                    returnbool = true;
                } else {
                    pending.push_back(I);
                }
            }
        }

        if (addNodes) {
            for (const IndexDimension& I : pending) {
                grid.AddPoint(I, false);
                returnbool = true;
            }
        }
    }

    return returnbool;
}

// Parallel-safe counterpart to refill_v2(), matching refill_NEU()'s contract: never calls
// AddPoint() (not thread-safe), instead appends missing points to the caller-owned
// `pointsToAdd` buffer. Callers running this inside an OpenMP region must ensure every
// depth T passed in already has its SingleDepthHashGrid created (see addPoints2_omp_v2()) -
// getGridForDepth() lazily inserts into a shared unordered_multimap on a miss, which is not
// thread-safe.
bool refill_NEU_v2(AdaptiveSparseGrid& grid, IndexDimension P, IndexDimension Q, Depth T, unsigned long p, unsigned long q, std::vector<IndexDimension>& pointsToAdd){

    bool returnbool = false;

    IndexDimension Imin, Imax;

    if (overlap(p, q, Imin, Imax, grid)) {

        SingleDepthHashGrid& depthGrid = grid.getMultiDepthHashGrid()->getGridForDepth(T);

        bool addNodes = false;
        std::vector<IndexDimension> pending;

        for (RectangularIteratorExact iter(Imin, Imax, T); iter.goon(); ++iter) {
            IndexDimension I = iter.getIndex();
            Depth Tneu(I);
            if (!(Tneu == T)) continue;

            unsigned long k;
            if (depthGrid.occupied(k, I)) {
                addNodes = true;
            } else if (I.isNotAtBoundary()) {
                if (addNodes) {
                    pointsToAdd.push_back(I);
                    returnbool = true;
                } else {
                    pending.push_back(I);
                }
            }
        }

        if (addNodes) {
            for (const IndexDimension& I : pending) {
                pointsToAdd.push_back(I);
                returnbool = true;
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

void addPoints2_omp(AdaptiveSparseGrid& dGrid){

    struct RefillCandidate {
        IndexDimension I1, I2;
        Depth T;
        unsigned long p, q;
    };

    // The set of active points (and which depth each belongs to) never changes during
    // addPoints2: refill/refill_NEU only ever add new points as hanging (AddPoint(I,false)),
    // and never flip an existing point's active flag. So this depth-pair / active-point-pair
    // enumeration is invariant across "changed" passes - build it once.
    std::vector<RefillCandidate> candidates;

    DepthList list(dGrid);
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
                    if (dGrid.getActiveTable()[mapping1[i]]) {
                        IndexDimension I1 = dGrid.getIndexOfTable(mapping1[i]);

                        for (size_t j = 0; j < mapping2.size(); j++) {
                            if (dGrid.getActiveTable()[mapping2[j]]) {
                                IndexDimension I2 = dGrid.getIndexOfTable(mapping2[j]);

                                candidates.push_back({I1, I2, T, mapping1[i], mapping2[j]});
                            }
                        }
                    }
                }
            }
        }
    }

    long numCandidates = (long) candidates.size();
    int maxThreads = omp_get_max_threads();

    for (bool changed = true; changed;) {
        changed = false;

        // Each thread accumulates into its own slot - no shared vector, no critical section.
        // AddPoint (not thread-safe: mutates shared hash tables) is never called inside the
        // parallel region; all points found this pass are applied afterward, serially.
        std::vector<std::vector<IndexDimension>> perThreadPoints(maxThreads);
        std::vector<char> perThreadChanged(maxThreads, 0);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<IndexDimension> &localPoints = perThreadPoints[tid];

            #pragma omp for schedule(dynamic)
            for (long idx = 0; idx < numCandidates; idx++) {
                const RefillCandidate &c = candidates[idx];
                if (refill_NEU(dGrid, c.I1, c.I2, c.T, c.p, c.q, localPoints)) {
                    perThreadChanged[tid] = 1;
                }
            }
        }

        for (int t = 0; t < maxThreads; t++) {
            if (perThreadChanged[t]) changed = true;
            for (const IndexDimension &I : perThreadPoints[t]) {
                dGrid.AddPoint(I, false);
            }
        }
    }
}

// Same depth-pair / active-point-pair enumeration as addPoints2(), but calling refill_v2()
// (single-pass box scan, depth-local occupied() lookups) instead of refill().

void addPoints2_v2(AdaptiveSparseGrid& dGrid){

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

                                   if (refill_v2(dGrid, I1, I2,T,mapping1[i],mapping2[j])) {

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
}

// Same candidate-list / per-thread-buffer structure as addPoints2_omp(), but calling
// refill_NEU_v2() instead of refill_NEU(). refill_NEU_v2() reads (and, on a miss, would
// insert into) MultiDepthHashGrid's shared depth map via getGridForDepth() - not safe to
// call concurrently on an unseen depth. Every candidate's T is therefore pre-created here,
// serially, before the parallel region; every T value refill_NEU_v2() will ever see across
// every "changed" pass comes from this same fixed candidate list (see addPoints2_omp()'s
// comment on why the candidate list itself is pass-invariant), and any point later added by
// AddPoint() lands in one of these same pre-created depth buckets (it's added at exactly
// depth T, see refill_NEU_v2()) - so no further creation is ever needed inside the loop.
void addPoints2_omp_v2(AdaptiveSparseGrid& dGrid){

    struct RefillCandidate {
        IndexDimension I1, I2;
        Depth T;
        unsigned long p, q;
    };

    std::vector<RefillCandidate> candidates;

    DepthList list(dGrid);
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
                    if (dGrid.getActiveTable()[mapping1[i]]) {
                        IndexDimension I1 = dGrid.getIndexOfTable(mapping1[i]);

                        for (size_t j = 0; j < mapping2.size(); j++) {
                            if (dGrid.getActiveTable()[mapping2[j]]) {
                                IndexDimension I2 = dGrid.getIndexOfTable(mapping2[j]);

                                candidates.push_back({I1, I2, T, mapping1[i], mapping2[j]});
                            }
                        }
                    }
                }
            }
        }
    }

    // Return value intentionally discarded: this call is here purely for its side effect.
    // getGridForDepth() inserts a new SingleDepthHashGrid into the shared depth map on a
    // miss - not safe to do concurrently. Forcing every needed depth to exist now, serially,
    // means every getGridForDepth() call made later from inside the parallel region below
    // (via refill_NEU_v2) is guaranteed to only ever hit the safe, read-only "found it" path.
    for (const RefillCandidate &c : candidates) {
        (void) dGrid.getMultiDepthHashGrid()->getGridForDepth(c.T);
    }

    long numCandidates = (long) candidates.size();
    int maxThreads = omp_get_max_threads();

    for (bool changed = true; changed;) {
        changed = false;

        std::vector<std::vector<IndexDimension>> perThreadPoints(maxThreads);
        std::vector<char> perThreadChanged(maxThreads, 0);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            std::vector<IndexDimension> &localPoints = perThreadPoints[tid];

            #pragma omp for schedule(dynamic)
            for (long idx = 0; idx < numCandidates; idx++) {
                const RefillCandidate &c = candidates[idx];
                if (refill_NEU_v2(dGrid, c.I1, c.I2, c.T, c.p, c.q, localPoints)) {
                    perThreadChanged[tid] = 1;
                }
            }
        }

        for (int t = 0; t < maxThreads; t++) {
            if (perThreadChanged[t]) changed = true;
            for (const IndexDimension &I : perThreadPoints[t]) {
                dGrid.AddPoint(I, false);
            }
        }
    }
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