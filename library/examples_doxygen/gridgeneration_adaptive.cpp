//! [start]
#include "source/sparseEXPDE.h"

int main(int argc, char **argv) {

    AdaptiveSparseGrid adaptiveSparseGrid;
    int level = 3;


    IndexDimension centerPoint;
    IndexDimension bottomLeft = centerPoint.leftSon(0).leftSon(1).leftSon(0).leftSon(1);
    IndexDimension topRight = centerPoint.rightSon(0).rightSon(1).rightSon(0).rightSon(1);

    adaptiveSparseGrid.AddPoint(bottomLeft);
    adaptiveSparseGrid.AddPoint(topRight);
    //! [start]

    //! [completegrid]
    adaptiveSparseGrid.completeDirichletGrid();
    //! [completegrid]

    //! [end]


    adaptiveSparseGrid.PrintActiveHanging(level);
    adaptiveSparseGrid.Print_gnu("grid.gnu");
    adaptiveSparseGrid.Print_vtk("grid.vtk");

    cout << adaptiveSparseGrid.getDOFS() << "  " << adaptiveSparseGrid.getHangingNodes() << endl;

    return 0;
}
//! [end]