
#include "myMath/myMath.h"
#include "abbrevi.h"
#include "myAssert.h"
#include "mympi.h"
#include "mycomm.h"
#include "indices/index.h"
#include "sgrid/sparseGrid.h"
#include "sgrid/depth.h"




#include "extemp/extempAlg.h"
#include "extemp/vector.h"
#include "extemp/multilevelvector.h"
#include "extemp/coordinate.h"
#include "extemp/shift.h"

#include "extemp/funktor.h"
#include "extemp/funktorExample.h"

#include "extemp/system.h"



#include "primes/prime.h"
#include "primes/prime.cc"

#include "sgrid/ListOfDepthOrderedGrids.h"
#include "sgrid/komponente.h"
#include "sgrid/multilevelSparseGrid.h"


#include "sgoperation/matrix_operations.h"
#include "sgoperation/matrix.h"


#include "fullgrid/FullGrid.h"

#include "tests/testing.h"

#include "cells/celldimension.h"
#include "cells/CellStructure.h"




#include "BasisTransformations/BasisTransformations.h"
#include "BasisTransformations/BasisTransformations_inhomog.h"

#include "MatrixVectorMultiplication/MatrixVectorInhomogen.h"
#include "MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "MatrixVectorMultiplication/RHS.h"





#include "stencils/PoissonStencil.h"
#include "stencils/MassStencil.h"
#include "stencils/Stencil.h"
//#include "stencils/StencilSGpp.h"
#include "stencils/StencilMC.h"

#include "localStiffnessMatrices/LocalStiffnessMatrices.h"
#include "localStiffnessMatrices/LocalStiffnessMatricesMemoryDistribution.h"
#include "localStiffnessMatrices/LocalStiffnessMatricesDynamicDistribution.h"
#include "localStiffnessMatrices/ADDLocalStiffnessMatrices.h"
#include "localStiffnessMatrices/LocalStiffnessMatricesFixedDistribution.h"



#include "applications/cg_method.h"
#include "applications/norms.h"
#include "applications/condition.h"
#include "applications/PDE.h"
#include "applications/H1_Norm.h"

#include "tests/testing.h"


#include "sgrid/refinement.h"
#include "sgrid/GridGeneration.h"





