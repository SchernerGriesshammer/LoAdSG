#include <filesystem>
#include <iomanip>

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
        leftCenter = leftCenter.leftSon(d);
    }

    IndexDimension rightCenter=centerPoint;
    for (int d=0; d<DimensionSparseGrid; ++d)
    {
        rightCenter = rightCenter.rightSon(d);
    }


    grid.AddPoint(leftCenter);
    grid.AddPoint(rightCenter);

    grid.completeDirichletGrid();


    AdaptiveSparseGrid grid_test;

    grid_test.AddPoint(leftCenter);
    grid_test.AddPoint(rightCenter);

    grid_test.completeDirichletGrid_timed();


    bool similar = gridsAreSimilar(grid, grid_test, false);
    cout << "grid and grid_test are " << (similar ? "similar" : "not similar") << endl;
    cout << "number of grid points" << grid.getDOFS() << " " << grid.getHangingNodes() << endl;

    AdaptiveSparseGrid grid_fast;

    grid_fast.AddPoint(leftCenter);
    grid_fast.AddPoint(rightCenter);

    grid_fast.completeDirichletGrid_fast_timed();

    bool fastSimilar = gridsAreSimilar(grid, grid_fast, false);
    cout << "grid and grid_fast are " << (fastSimilar ? "similar" : "not similar") << endl;



    return 0;
}