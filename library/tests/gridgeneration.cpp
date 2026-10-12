#include <filesystem>
#include <iomanip>
#include <omp.h>

#include "../source/sparseEXPDE.h"


void printIndex(IndexDimension I) {
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        cout << IndexOneD(I.getIndex(d)).coordinate() << (d + 1 < DimensionSparseGrid ? ", " : "");
    }
}

/**
 * Compares two adaptive sparse grids and reports whether they are similar,
 * i.e. whether they contain the same set of grid points and whether each
 * point has the same status (active vs. hanging) in both grids.
 *
 * @param gridA first grid
 * @param gridB second grid
 * @param verbose if true, differences are printed to stdout
 * @return true if both grids contain the same points with the same
 *         active/hanging status, false otherwise
 **/
bool gridsAreSimilar(AdaptiveSparseGrid &gridA, AdaptiveSparseGrid &gridB, bool verbose = true) {
    bool similar = true;

    if (gridA.getDOFS() != gridB.getDOFS()) {
        similar = false;
        if (verbose)
            cout << "Number of active nodes differs: " << gridA.getDOFS()
                 << " vs. " << gridB.getDOFS() << endl;
    }

    if (gridA.getHangingNodes() != gridB.getHangingNodes()) {
        similar = false;
        if (verbose)
            cout << "Number of hanging nodes differs: " << gridA.getHangingNodes()
                 << " vs. " << gridB.getHangingNodes() << endl;
    }

    // every point of gridA must be in gridB with the same active/hanging status
    for (unsigned long k = 0; k < gridA.getMaximalOccupiedSecondTable(); ++k) {
        if (gridA.getSecondTable()[k] == 0) continue;

        IndexDimension I = gridA.getIndexOfTable(k);
        bool activeInA = gridA.getActiveTable()[k];

        unsigned long indexOfData;
        if (!gridB.occupied(indexOfData, I)) {
            similar = false;
            if (verbose) {
                cout << "Point (";
                printIndex(I);
                cout << ") is present in first grid but missing in second grid." << endl;
            }
            continue;
        }

        bool activeInB = gridB.getActiveTable()[indexOfData];
        if (activeInA != activeInB) {
            similar = false;
            if (verbose) {
                cout << "Point (";
                printIndex(I);
                cout << ") is " << (activeInA ? "active" : "hanging") << " in first grid but "
                     << (activeInB ? "active" : "hanging") << " in second grid." << endl;
            }
        }
    }

    // every point of gridB must also be in gridA (points only in gridB)
    for (unsigned long k = 0; k < gridB.getMaximalOccupiedSecondTable(); ++k) {
        if (gridB.getSecondTable()[k] == 0) continue;

        IndexDimension I = gridB.getIndexOfTable(k);
        unsigned long indexOfData;
        if (!gridA.occupied(indexOfData, I)) {
            similar = false;
            if (verbose) {
                cout << "Point (";
                printIndex(I);
                cout << ") is present in second grid but missing in first grid." << endl;
            }
        }
    }

    if (verbose) {
        if (similar)
            cout << "Grids are similar: same active and hanging nodes." << endl;
        else
            cout << "Grids are NOT similar." << endl;
    }

    return similar;
}


int main(int argc, char **argv) {





    AdaptiveSparseGrid grid;


    IndexDimension centerPoint;
    IndexDimension leftCenter=centerPoint;
    for (int d=0; d<DimensionSparseGrid; ++d)
    {
        leftCenter = leftCenter.leftSon(d).leftSon(d).leftSon(d).leftSon(d);
    }

    IndexDimension rightCenter=centerPoint;
    for (int d=0; d<DimensionSparseGrid; ++d)
    {
    rightCenter = rightCenter.rightSon(d).rightSon(d).rightSon(d).rightSon(d);
    }


    grid.AddPoint(leftCenter);
    grid.AddPoint(rightCenter);



    grid.completeDirichletGrid();


    cout << "active nodes " <<endl;
    for (unsigned long k=0; k< grid.getMaximalOccupiedSecondTable(); k++)
    {
        IndexDimension I = grid.getIndexOfTable(k);
        if (grid.getActiveTable()[k]) I.PrintCoord();
    }

    cout << endl;
    cout << " hanging nodes " << endl;
    for (unsigned long k=0; k < grid.getMaximalOccupiedSecondTable(); k++)
    {
        IndexDimension I = grid.getIndexOfTable(k);
        if (!grid.getActiveTable()[k]) I.PrintCoord();
    }

    exit(1);


    AdaptiveSparseGrid grid_test;

    grid_test.AddPoint(leftCenter);
    grid_test.AddPoint(rightCenter);

    //grid_test.completeDirichletGrid_timed();


    bool similar = gridsAreSimilar(grid, grid_test, false);
    cout << "grid and grid_test are " << (similar ? "similar" : "not similar") << endl;
    cout << "number of grid points" << grid.getDOFS() << " " << grid.getHangingNodes() << endl;

    AdaptiveSparseGrid grid_fast;

    grid_fast.AddPoint(leftCenter);
    grid_fast.AddPoint(rightCenter);

    //grid_fast.completeDirichletGrid_fast_timed();

    bool fastSimilar = gridsAreSimilar(grid, grid_fast, false);
    cout << "grid and grid_fast are " << (fastSimilar ? "similar" : "not similar") << endl;

    AdaptiveSparseGrid grid_omp;

    grid_omp.AddPoint(leftCenter);
    grid_omp.AddPoint(rightCenter);

    cout << "OMP max threads: " << omp_get_max_threads() << endl;

    grid_omp.completeDirichletGrid_omp_timed();

    bool ompSimilar = gridsAreSimilar(grid, grid_omp, false);
    cout << "grid and grid_omp are " << (ompSimilar ? "similar" : "not similar") << endl;
    cout << "grid_omp number of grid points " << grid_omp.getDOFS() << " " << grid_omp.getHangingNodes() << endl;

    AdaptiveSparseGrid grid_v2;

    grid_v2.AddPoint(leftCenter);
    grid_v2.AddPoint(rightCenter);

    //grid_v2.completeDirichletGrid_v2_timed();

    bool v2Similar = gridsAreSimilar(grid_omp, grid_v2, false);
    cout << "grid and grid_v2 are " << (v2Similar ? "similar" : "not similar") << endl;
    cout << "grid_v2 number of grid points " << grid_v2.getDOFS() << " " << grid_v2.getHangingNodes() << endl;

    AdaptiveSparseGrid grid_omp_v2;

    grid_omp_v2.AddPoint(leftCenter);
    grid_omp_v2.AddPoint(rightCenter);

    grid_omp_v2.completeDirichletGrid_omp_v2_timed();

    bool ompV2Similar = gridsAreSimilar(grid_omp, grid_omp_v2, false);
    cout << "grid and grid_omp_v2 are " << (ompV2Similar ? "similar" : "not similar") << endl;
    cout << "grid_omp_v2 number of grid points " << grid_omp_v2.getDOFS() << " " << grid_omp_v2.getHangingNodes() << endl;

    return 0;
}