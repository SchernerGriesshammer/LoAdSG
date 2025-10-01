//
// Created by to35jepo on 9/19/24.
//

#ifndef RUN_STENCILMC_H
#define RUN_STENCILMC_H


#include <random>
#include "Stencil.h"


template <typename F>
class StencilMC: public StencilTemplate{

public:

    StencilMC(AdaptiveSparseGrid& grid_, const F &varCoeff_, int c_) : StencilTemplate(grid_), varCoeff(varCoeff_), c(c_){
        dofs=grid.getDOFS();
        gen =std::mt19937(rd());
    };


    int c=1;


    inline double integration(CellDimension& cell, CellIndexDirection& dirP, CellIndexDirection& dirQ){
        FunctionClass functionClass(cell,dirP,dirQ);
        double value = 1.0;

        int numberMC = 1;
        int l=0;

        for(int d=0; d<DimensionSparseGrid;d++){
            l+=cell.getDepth(d)-1;

        }


        int vol_cell= POW2(unsigned (l));


        if(vol_cell>dofs) numberMC=1;
        else numberMC= dofs/vol_cell;

        numberMC*=c;







        double sum=0.0;


        for(int j=0; j<numberMC;j++){
            double x[DimensionSparseGrid];
            for(int d=0;d<DimensionSparseGrid;d++){
                x[d]=std::uniform_real_distribution<>(cell.getLeftBoundary(d),cell.getRightBoundary(d))(gen);
            }
            double vC = varCoeff(x)*functionClass.evaluation(x);

            sum+=vC;
        }
        value = (1.0/vol_cell)*(sum/double(numberMC));

        return value;



    }
    class FunctionClass{


    public:
        FunctionClass (CellDimension& cell_, CellIndexDirection& dirP_, CellIndexDirection &dirQ_): cell(cell_),dirP(dirP_),dirQ(dirQ_){

            for(int i=0; i<DimensionSparseGrid; i++){
                double h = cell.getRightBoundary(i)-cell.getLeftBoundary(i);

                meshwidth[i]=h;

            }

        }
        double evaluation(double* x){






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





            return value;
        };






    private:

        double meshwidth[DimensionSparseGrid];
        CellDimension cell;
        CellIndexDirection dirP;
        CellIndexDirection dirQ;


    };
    F varCoeff;
    int dofs;
    std::random_device rd;
    std::uniform_real_distribution<> dist[DimensionSparseGrid];
    std::mt19937 gen;
};


template <typename F>
class RHS_MC: public StencilTemplate{

public:

    RHS_MC(AdaptiveSparseGrid& grid_, const F &varCoeff_, int c_) : StencilTemplate(grid_), varCoeff(varCoeff_), c(c_){
        dofs=grid.getDOFS();
        gen =std::mt19937(rd());
    };


    int c=1;


    inline double integration(CellDimension& cell, CellIndexDirection& dirP, CellIndexDirection& dirQ){
        FunctionClass functionClass(cell,dirP,dirQ);
        double value = 1.0;

        int numberMC = 1;
        int l=0;

        for(int d=0; d<DimensionSparseGrid;d++){
            l+=cell.getDepth(d)-1;

        }


        int vol_cell= POW2(unsigned (l));


        if(vol_cell>dofs) numberMC=1;
        else numberMC= dofs/vol_cell;

        numberMC*=c;







        double sum=0.0;


        for(int j=0; j<numberMC;j++){
            double x[DimensionSparseGrid];
            for(int d=0;d<DimensionSparseGrid;d++){
                x[d]=std::uniform_real_distribution<>(cell.getLeftBoundary(d),cell.getRightBoundary(d))(gen);
            }
            double vC = varCoeff(x)*functionClass.evaluation(x);

            sum+=vC;
        }
        value = (1.0/vol_cell)*(sum/double(numberMC));

        return value;



    }
    class FunctionClass{


    public:
        FunctionClass (CellDimension& cell_, CellIndexDirection& dirP_, CellIndexDirection &dirQ_): cell(cell_),dirP(dirP_),dirQ(dirQ_){

            for(int i=0; i<DimensionSparseGrid; i++){
                double h = cell.getRightBoundary(i)-cell.getLeftBoundary(i);

                meshwidth[i]=h;

            }

        }
        double evaluation(double* x){

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





            return value;
        };






    private:

        double meshwidth[DimensionSparseGrid];
        CellDimension cell;
        CellIndexDirection dirP;
        CellIndexDirection dirQ;


    };
    F varCoeff;
    int dofs;
    std::random_device rd;
    std::uniform_real_distribution<> dist[DimensionSparseGrid];
    std::mt19937 gen;
};


#endif //RUN_STENCILMC_H
