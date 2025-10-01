#include "source/sparseEXPDE.h"

int main(int argc, char **argv) {

    AdaptiveSparseGrid regular_grid;

    int level=3;
    IndexDimension centerPoint;

    regular_grid.AddRecursiveSonsOfPoint(centerPoint,level);

    regular_grid.PrintActiveHanging(level);
    regular_grid.PrintActiveHanging(2);
    regular_grid.Print_gnu("grid.gnu");

    return 0;
}