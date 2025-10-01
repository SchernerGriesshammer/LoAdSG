//
// Created by to35jepo on 5/3/24.
//

#ifndef RUN_STENCILSGPP_H
#define RUN_STENCILSGPP_H
#include "Stencil.h"

#include "../../../SGpp/base/src/sgpp_base.hpp"


enum QuadratureType {
    SGadaptive,
    SGregular,
    MonteCarlo,
    MonteCarlo_Hierarchichal,
    MonteCarloFixed
};


template <typename F>
class StencilSGpp: public StencilTemplate {
public:


    int level;
    bool adaptive;
    int numberofgridpoints=0;
    CellDimension cellDimension;
    QuadratureType quadrature_type;
    int gridpoints=2147483647;
    std::unique_ptr<sgpp::base::Grid> grid_MC;
    int maxlevel=0;
    int getMinimumGridPointsOfACell(){
        int localvalue = gridpoints;
        int global_min;
        MPI_Reduce(&localvalue, &global_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        return global_min;
    };

    int mcPoints=1;

    StencilSGpp(AdaptiveSparseGrid& grid_, const F &varCoeff_, int level_=1,QuadratureType quadrature_type_=SGadaptive) : StencilTemplate(grid_), varCoeff(varCoeff_), level(level_), quadrature_type(quadrature_type_){

        if(quadrature_type==MonteCarloFixed) {
            grid_MC = std::unique_ptr<sgpp::base::Grid>(sgpp::base::Grid::createLinearGrid(DimensionSparseGrid));
            grid_MC->getGenerator().regular(level);


        }
        //gridpoints = grid_MC->getSize();
        mcPoints=grid.getDOFS();

    };




    inline double returnValue(const IndexDimension &Index, const MultiDimCompass &mc){cout << "StencilSGPP returnValue not implmented" << endl; exit(1);}

    inline double integration(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {
        double value = 0.0;


        if(quadrature_type == SGadaptive) {
            value = integration_adaptive(cell,dirP,dirQ);
        } else if(quadrature_type == SGregular){
            value = integration_regular(cell,dirP,dirQ);
        } else if(quadrature_type == MonteCarlo){
            value = integration_MC(cell,dirP,dirQ);
        } else if(quadrature_type == MonteCarlo_Hierarchichal){
            value = integration_MC_AD(cell,dirP,dirQ);
        }else if(quadrature_type == MonteCarloFixed) {
            value = integration_MC_fixedCell(cell, dirP, dirQ);
        }


        return value;
    }
    inline double integration_regular(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);

    inline double integrationLocallyAdaptive(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);


    inline double integration_adaptive(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);

    inline double integration_gausslegendre(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);



    inline double integration_MC(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);


    inline double integration_MC_fixedCell(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);

    inline double integration_MC_fixedCell2(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);



    inline double integration_MC_fixed(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);



    inline double integration_MC_AD(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ);
protected:

    F varCoeff;



    //in dieser Klasse soll die Funktion erstellt werden, die dann über das Einheitsquadrat integriert werden soll. Stichwort: Substitution
    class FunctionClass{


    public:
        FunctionClass (CellDimension& cell_, CellIndexDirection& dirP_, CellIndexDirection &dirQ_,F varCoeff_ ): cell(cell_),dirP(dirP_),dirQ(dirQ_),varCoeff(varCoeff_){
            // transformation of variable xi gives
            // J : xi -> (bi-ai)xi+ai
            // det(J')=(b1-a1)*...*(bd-ad)
            det_J_prime=1.0;
            for(int i=0; i<DimensionSparseGrid; i++){
                double h = cell.getRightBoundary(i)-cell.getLeftBoundary(i);
                det_J_prime*=h;
                meshwidth[i]=h;

            }

        }
        double evaluation(int dim, double* x, void* clientdata){
            // transform variable
            // J : xi -> (bi-ai)xi+ai
            // det(J')=(b1-a1)*...*(bd-ad)

            for(int i=0; i<DimensionSparseGrid; i++){
                x[i]=meshwidth[i]*x[i]+cell.getLeftBoundary(i);
            }






            double value = 1.0;

            for(int i = 0; i < DimensionSparseGrid; ++i )
            {

                // === Muliply by the first base function ===
                if(dirP.getDirection(i)==Left){
                    // p(x) = (x-b)/(a-b)=(b-x)/h
                    value *= (cell.getRightBoundary(i)-x[i]) / (meshwidth[i]);

                }else{
                    // p(x) = (x-a)/(b-a)
                    value *= (x[i] - cell.getLeftBoundary(i)) / (meshwidth[i]);
                }

                // === Muliply by the second base function ===
                if(dirQ.getDirection(i)==Left){
                    // p(x) = (x-b)/(a-b)=(b-x)/h
                    value *= (cell.getRightBoundary(i)-x[i]) / (meshwidth[i]);


                }else{
                    // p(x) = (x-a)/(b-a)
                    value *= (x[i] - cell.getLeftBoundary(i)) / (meshwidth[i]);
                }


            }


            value *= varCoeff(x);


            return value;
        };




        double detJPrime(){return det_J_prime;};

        static double evaluationWrapper(int dim, double* x, void* clientdata) {
            return static_cast<FunctionClass*>(clientdata)->evaluation(dim, x, nullptr);
        }
    private:
        double det_J_prime = 1.0;
        double meshwidth[DimensionSparseGrid];
        CellDimension cell;
        CellIndexDirection dirP;
        CellIndexDirection dirQ;
        F varCoeff;
    };



};

template<typename F>
double StencilSGpp<F>::integration_MC_AD(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int levels[DimensionSparseGrid];
    int maxlevel=0;
    for(int d=0; d<DimensionSparseGrid;d++){
        levels[d]=grid.getMaxDepth(d)-(cell.getDepth(d));
        maxlevel<levels[d]?maxlevel=levels[d]:maxlevel;
    }





    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level
    grid_sgpp->getGenerator().regular(1);


    sgpp::base::DataVector alpha(gridStorage.getSize());
    sgpp::base::DataVector refined(gridStorage.getSize());


    std::vector<size_t> addedPoints;


    for (int step = 0; step < maxlevel; step++) {

        int alpha_c=0;
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            alpha[i]=1.0;


            // do not allow refinement if max level in one direction is reached
            for(size_t k=0; k< DimensionSparseGrid; k++){
                if(gp.getLevel(k)> levels[k]){
                    alpha[i]=0.0;
                    break;
                }

            }



            // ensure that grid point can only be refined once
            if(refined[i]>1.0) alpha[i]=0.0;


            if(alpha[i]>0.0) {
                alpha_c++;

                // this means that this grid point can't be refined in the next refinement step
                refined[i] = 2.0;
            }


        }

        if (alpha_c<1)break;

        sgpp::base::SurplusRefinementFunctor functor(alpha,alpha_c, 0.5);
        grid_sgpp->getGenerator().refine(functor, &addedPoints);


        alpha.resize(gridStorage.getSize());
        refined.resize(gridStorage.getSize());

        if(addedPoints.size()<1)break;
        addedPoints.clear();
    }

    double p[DimensionSparseGrid];
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
        sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
        for(int j=0; j<DimensionSparseGrid;j++){
            p[j] = gp.getStandardCoordinate(j);
        }
        alpha[i] = function.evaluation(DimensionSparseGrid, p, NULL);
    }

    if(numberofgridpoints<gridStorage.getSize()){
        numberofgridpoints=gridStorage.getSize();
        cellDimension = cell;
    }



    std::unique_ptr<sgpp::base::OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid_sgpp))
            ->doHierarchisation(alpha);

    sgpp::base::OperationQuadratureMC opMC(*grid_sgpp,gridStorage.getSize());




    //cout << " do mc" << endl;
    double res = opMC.doQuadrature(alpha);
    //double res = opMC.doQuadratureFunc(FunctionClass::evaluationWrapper, &function);



    res *= function.detJPrime();

    return res;
}

template<typename F>
double StencilSGpp<F>::integration_MC_fixed(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {

    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int dim = DimensionSparseGrid;



    std::unique_ptr<sgpp::base::Grid> grid2(sgpp::base::Grid::createLinearBoundaryGrid(DimensionSparseGrid));




    sgpp::base::OperationQuadratureMC opMC(*grid2, grid.getDOFS());

    double res = opMC.doQuadratureFunc(FunctionClass::evaluationWrapper, &function);
    res *= function.detJPrime();

    return res;
}

template<typename F>
double StencilSGpp<F>::integration_MC_fixedCell(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int l=0;
    for(int d=0; d<DimensionSparseGrid;d++){
        l+=cell.getDepth(d)-1;
    }

    int vol_cell= POW2(unsigned (l));

    if(vol_cell>grid_MC->getSize()) mcPoints=1;
    else mcPoints = grid_MC->getSize()/vol_cell;



    sgpp::base::OperationQuadratureMC opMC(*grid_MC,mcPoints);

    double res = opMC.doQuadratureFunc(FunctionClass::evaluationWrapper, &function);



    res *= function.detJPrime();

    return res;
}

template<typename F>
double StencilSGpp<F>::integration_MC_fixedCell2(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);



    mcPoints = grid_MC->getSize();
    sgpp::base::OperationQuadratureMC opMC(*grid_MC,mcPoints);

    double res = opMC.doQuadratureFunc(FunctionClass::evaluationWrapper, &function);



    res *= function.detJPrime();

    return res;
}



template<typename F>
double StencilSGpp<F>::integration_MC(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int levels[DimensionSparseGrid];
    int maxlevel=0;
    for(int d=0; d<DimensionSparseGrid;d++){
        levels[d]=grid.getMaxDepth(d)-(cell.getDepth(d));
        maxlevel<levels[d]?maxlevel=levels[d]:maxlevel;
    }





    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level
    grid_sgpp->getGenerator().regular(1);


    sgpp::base::DataVector alpha(gridStorage.getSize());
    sgpp::base::DataVector refined(gridStorage.getSize());


    std::vector<size_t> addedPoints;


    for (int step = 0; step < maxlevel; step++) {

        int alpha_c=0;
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            alpha[i]=1.0;


            // do not allow refinement if max level in one direction is reached
            for(size_t k=0; k< DimensionSparseGrid; k++){
                if(gp.getLevel(k)> levels[k]){
                    alpha[i]=0.0;
                    break;
                }

            }



            // ensure that grid point can only be refined once
            if(refined[i]>1.0) alpha[i]=0.0;


            if(alpha[i]>0.0) {
                alpha_c++;

                // this means that this grid point can't be refined in the next refinement step
                refined[i] = 2.0;
            }


        }

        if (alpha_c<1)break;

        sgpp::base::SurplusRefinementFunctor functor(alpha,alpha_c, 0.5);
        grid_sgpp->getGenerator().refine(functor, &addedPoints);


        alpha.resize(gridStorage.getSize());
        refined.resize(gridStorage.getSize());

        if(addedPoints.size()<1)break;
        addedPoints.clear();
    }

/*        double p[DimensionSparseGrid];
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            for(int j=0; j<DimensionSparseGrid;j++){
                p[j] = gp.getStandardCoordinate(j);
            }
            alpha[i] = function.evaluation(DimensionSparseGrid, p, NULL);
        }

        if(numberofgridpoints<gridStorage.getSize()){
            numberofgridpoints=gridStorage.getSize();
            cellDimension = cell;
        }



        std::unique_ptr<sgpp::base::OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid_sgpp))
                ->doHierarchisation(alpha);*/

    sgpp::base::OperationQuadratureMC opMC(*grid_sgpp,gridStorage.getSize());
/*        std::unique_ptr<sgpp::base::OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid_sgpp))
                ->doHierarchisation(alpha);*/



    //cout << " do mc" << endl;
    //double res = opMC.doQuadrature(alpha);
    double res = opMC.doQuadratureFunc(FunctionClass::evaluationWrapper, &function);



    res *= function.detJPrime();

    return res;
}

template<typename F>
double
StencilSGpp<F>::integration_gausslegendre(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int levels[DimensionSparseGrid];
    int maxlevel=0;
    for(int d=0; d<DimensionSparseGrid;d++){
        levels[d]=grid.getMaxDepth(d)-(cell.getDepth(d));
        maxlevel<levels[d]?maxlevel=levels[d]:maxlevel;
    }





    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level
    grid_sgpp->getGenerator().regular(1);


    sgpp::base::DataVector alpha(gridStorage.getSize());
    sgpp::base::DataVector refined(gridStorage.getSize());


    std::vector<size_t> addedPoints;


    for (int step = 0; step < maxlevel; step++) {

        int alpha_c=0;
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            alpha[i]=1.0;


            // do not allow refinement if max level in one direction is reached
            for(size_t k=0; k< DimensionSparseGrid; k++){
                if(gp.getLevel(k)> levels[k]){
                    alpha[i]=0.0;
                    break;
                }

            }



            // ensure that grid point can only be refined once
            if(refined[i]>1.0) alpha[i]=0.0;


            if(alpha[i]>0.0) {
                alpha_c++;

                // this means that this grid point can't be refined in the next refinement step
                refined[i] = 2.0;
            }


        }

        if (alpha_c<1)break;

        sgpp::base::SurplusRefinementFunctor functor(alpha,alpha_c, 0.5);
        grid_sgpp->getGenerator().refine(functor, &addedPoints);


        alpha.resize(gridStorage.getSize());
        refined.resize(gridStorage.getSize());

        if(addedPoints.size()<1)break;
        addedPoints.clear();
    }


    double p[DimensionSparseGrid];
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
        sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
        for(int j=0; j<DimensionSparseGrid;j++){
            p[j] = gp.getStandardCoordinate(j);
        }
        alpha[i] = function.evaluation(DimensionSparseGrid, p, NULL);
    }


    double res;



    res*= function.detJPrime();

    return res;
}

template<typename F>
double StencilSGpp<F>::integration_adaptive(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int levels[DimensionSparseGrid];
    int maxlevel=0;
    for(int d=0; d<DimensionSparseGrid;d++){
        levels[d]=grid.getMaxDepth(d)-(cell.getDepth(d));
        maxlevel<levels[d]?maxlevel=levels[d]:maxlevel;
    }





    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level
    grid_sgpp->getGenerator().regular(1);


    sgpp::base::DataVector alpha(gridStorage.getSize());
    sgpp::base::DataVector refined(gridStorage.getSize());


    std::vector<size_t> addedPoints;


    for (int step = 0; step < maxlevel; step++) {

        int alpha_c=0;
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            alpha[i]=1.0;


            // do not allow refinement if max level in one direction is reached
            for(size_t k=0; k< DimensionSparseGrid; k++){
                if(gp.getLevel(k)> levels[k]){
                    alpha[i]=0.0;
                    break;
                }

            }



            // ensure that grid point can only be refined once
            if(refined[i]>1.0) alpha[i]=0.0;


            if(alpha[i]>0.0) {
                alpha_c++;

                // this means that this grid point can't be refined in the next refinement step
                refined[i] = 2.0;
            }


        }

        if (alpha_c<1)break;

        sgpp::base::SurplusRefinementFunctor functor(alpha,alpha_c, 0.5);
        grid_sgpp->getGenerator().refine(functor, &addedPoints);


        alpha.resize(gridStorage.getSize());
        refined.resize(gridStorage.getSize());

        if(addedPoints.size()<1)break;
        addedPoints.clear();
    }


    double p[DimensionSparseGrid];
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
        sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
        for(int j=0; j<DimensionSparseGrid;j++){
            p[j] = gp.getStandardCoordinate(j);
        }
        alpha[i] = function.evaluation(DimensionSparseGrid, p, NULL);
    }

    if(numberofgridpoints<gridStorage.getSize()){
        numberofgridpoints=gridStorage.getSize();
        cellDimension = cell;
    }


    // Prints Grid
/*
        sgpp::base::GridPrinter gridPrinter(*grid_sgpp);
        string filename = "sgpp_grid.gnu";
        gridPrinter.printSparseGrid(alpha,filename,false);
*/

    std::unique_ptr<sgpp::base::OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid_sgpp))
            ->doHierarchisation(alpha);

    // direct quadrature
    std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
            sgpp::op_factory::createOperationQuadrature(*grid_sgpp));
    double res = opQ->doQuadrature(alpha);

    res *= function.detJPrime();

    return res;
}

template<typename F>
double
StencilSGpp<F>::integrationLocallyAdaptive(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);


    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearBoundaryGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level
    grid_sgpp->getGenerator().regular(level);


    sgpp::base::DataVector alpha(gridStorage.getSize());
    sgpp::base::DataVector refined(gridStorage.getSize());

    double p[DimensionSparseGrid];
    sgpp::base::DataVector funEvals(gridStorage.getSize());
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
        sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
        for(size_t k=0; k<DimensionSparseGrid; k++){
            p[k]=gp.getStandardCoordinate(k);
        }
        funEvals[i] = function.evaluation(DimensionSparseGrid,p,NULL);
        refined[i]=0.0;

    }

    std::vector<size_t> addedPoints;

    alpha.copyFrom(funEvals);
    sgpp::op_factory::createOperationHierarchisation(*grid_sgpp)->doHierarchisation(alpha);





    int number_refinements=10;
    for (int step = 0; step < number_refinements-1; step++) {


        int alpha_c=0;
        for (size_t i = 0; i < gridStorage.getSize(); i++) {
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);

            for(size_t k=0; k< DimensionSparseGrid; k++){
                if(gp.getLevel(k)>= grid.getMaxDepth(int(k))-cell.getDepth(int(k))){

                    alpha[i]=0.0;
                    break;
                }

            }

            if(refined[i]>1.0) alpha[i]=0.0;


            if(alpha[i]>0.0) {
                alpha_c++;
                refined[i] = 2.0;
            }





        }

        cout << " alpha > 0 " << alpha_c << endl;



        sgpp::base::SurplusRefinementFunctor functor(alpha,alpha_c, 0.0);
        grid_sgpp->getGenerator().refine(functor, &addedPoints);

        cout << " added points " << addedPoints.size() << endl;







        alpha.resize(gridStorage.getSize());
        funEvals.resize(gridStorage.getSize());
        refined.resize(gridStorage.getSize());


        for (size_t i = 0; i < addedPoints.size(); i++) {
            size_t seq = addedPoints[i];
            sgpp::base::GridPoint& gp = gridStorage.getPoint(seq);
            for(size_t k=0; k<DimensionSparseGrid; k++){
                p[k]=gp.getStandardCoordinate(k);
            }
            funEvals[seq] = function.evaluation(DimensionSparseGrid,p,NULL);
        }



        alpha.copyFrom(funEvals);
        sgpp::op_factory::createOperationHierarchisation(*grid_sgpp)->doHierarchisation(alpha);

        int max=0;
        int sum=0;
        for(size_t i=0; i<gridStorage.getSize(); i++){
            sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
            for(int d=0; d<DimensionSparseGrid;d++){
                sum+=gp.getLevel(d);
            }
            if(max<sum)max=sum;
            sum=0;
        }

        cout << "max " << max << endl;





        if(addedPoints.size()<1)break;
        addedPoints.clear();
    }


    cout << " final points "  << gridStorage.getSize() << endl;




    // direct quadrature
    std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
            sgpp::op_factory::createOperationQuadrature(*grid_sgpp));
    double res = opQ->doQuadrature(alpha);

    res *= function.detJPrime();

    return res;
}

template<typename F>
double StencilSGpp<F>::integration_regular(CellDimension &cell, CellIndexDirection &dirP, CellIndexDirection &dirQ) {


    FunctionClass function(cell,dirP,dirQ,varCoeff);

    int dim = DimensionSparseGrid;
    std::unique_ptr<sgpp::base::Grid> grid_sgpp(sgpp::base::Grid::createLinearBoundaryGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid_sgpp->getStorage();

    // create regular grid, level 3

    grid_sgpp->getGenerator().regular(level);


    //ab hier soll SGpp anfangen, function.evaluation() ist dann die Funktion die von SGpp integriert werden soll

    sgpp::base::DataVector alpha(gridStorage.getSize());
    double p[DimensionSparseGrid];
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
        sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
        for(int j=0; j<DimensionSparseGrid;j++){
            p[j] = gp.getStandardCoordinate(j);
        }
        alpha[i] = function.evaluation(DimensionSparseGrid, p, NULL);
    }

    std::unique_ptr<sgpp::base::OperationHierarchisation>(sgpp::op_factory::createOperationHierarchisation(*grid_sgpp))
            ->doHierarchisation(alpha);

    // direct quadrature
    std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
            sgpp::op_factory::createOperationQuadrature(*grid_sgpp));
    double res = opQ->doQuadrature(alpha);

    res *= function.detJPrime();




    return res;
}


#endif //RUN_STENCILSGPP_H
