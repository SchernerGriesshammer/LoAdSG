#ifndef SGINTEGRATOR_H
#define SGINTEGRATOR_H

#include <iostream>
#include <algorithm>
#include <vector>
#include "curvedPoissonIntegrator.h"
#include "../math_lib.h"
#include <queue>
#include <cmath>
#include <string>

typedef double (IntegratorPoissonCurved::*Curvedfunc)(D3vector& x);

class SGPoint{
private:
    D3vector coords;//coordinates of the point
    D3vector levels;//indicates the level of refinement
    D3vector locs;//indicates the location
    D3vector supp_sizes;//size of the support of the basis function along each dimension

    D3vector limit1,limit2;

public:
    SGPoint(){

    }

    SGPoint(D3vector levels_,D3vector locs_,\
            D3vector &limit1, D3vector &limit2):levels(levels_),locs(locs_){
        this->limit1 = limit1;
        this->limit2 = limit2;

        //Assuming level starts from 0. Level 0 ->  Boundary basis functions
        for(int i=0;i<3;i++){
            if(levels[i] == 0){
                supp_sizes[i] = limit2[i]-limit1[i];
                if(locs[i] == 0){
                    coords[i] = limit1[i];
                }
                else if(locs[i] == 2){
                    coords[i] = limit2[i];
                }
            }
            else{
                supp_sizes[i] = (limit2[i]-limit1[i])/pow(2,levels[i]-1);
                coords[i] = limit1[i]+locs[i]*supp_sizes[i]*0.5;
            }
        }

    }

    SGPoint& operator=(SGPoint& point){
        this->coords = point.coords;
        this->levels = point.levels;
        this->locs = point.locs;
        this->supp_sizes = point.supp_sizes;

        return *this;
    }

    void Print(){
        std::cout<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
    }

    double volume(int dim);

    D3vector getCoords(){
        return coords;
    }

    D3vector getSizes(){
        return supp_sizes;
    }

    D3vector getLevels(){
        return levels;
    }

    D3vector getLocs(){
        return locs;
    }

    void reset(D3vector levels_,D3vector locs_,D3vector &limit1, D3vector &limit2);
};


class SGIntegrator{
private:
    std::queue<SGPoint> refine_queue;
    std::vector<SGPoint> refined_set;
    double tol;//Refinement criteria
    int dim;//dimension

    //Limits of integration
    D3vector limit1;
    D3vector limit2;
public:
    SGIntegrator(double tol_, int dim_):tol(tol_),dim(dim_){

    }

    //Set the tolerance
    void  setTol(double &tol_){
        this->tol =  tol_;
    }

    //Check whether a given point lies on the boundary
    bool isBoundaryPoint(SGPoint &p);

    //Refine the grid by adding neighbours
    void refine(SGPoint &p);

    void decideMinMax(D3vector &min,  D3vector &max, D3vector &levels);

    //Find the contribution to the total integral from the given subspace of basis functions
    double subspaceIntegral(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate, int level_x, int level_y, \
                            int level_z);

    //Find the hierarchical surplus
    double findHierarchicalSurplus(IntegratorPoissonCurved& ipc, SGPoint &p, Curvedfunc func);

    //Perform the adaptive integration with the given integral limits
    double integrate(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate,D3vector &limit1, D3vector &limit2);

    //Perform the integration uising a classical smolyak grid
    double smolyakIntegrate(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate, D3vector &limit1, D3vector &limit2, int max_level);

};

inline bool operator==(SGPoint point1, SGPoint point2){
    return(point1.getCoords() == point2.getCoords() &&\
           point1.getLevels() == point2.getLevels() &&\
           point1.getLocs() == point2.getLocs());
}

//a utility function used in findHierarchicalSurplus
int corn(int index){
    if(index == 0){
        return -1;
    }
    else{
        return 1;
    }
}

void SGPoint::reset(D3vector levels_,D3vector locs_,D3vector &limit1, D3vector &limit2){
    this->levels = levels_;
    this->locs = locs_;

    this->limit1 = limit1;
    this->limit2 = limit2;

    //Assuming level starts from 0. Level 0 =>  Boundary basis functions
    for(int i=0;i<3;i++){
        if(levels[i] == 0){
            supp_sizes[i] = limit2[i]-limit1[i];
            if(locs[i] == 0){
                coords[i] = limit1[i];
            }
            else if(locs[i] == 2){
                coords[i] = limit2[i];
            }
        }
        else{
            supp_sizes[i] = (limit2[i]-limit1[i])/pow(2,levels[i]-1);
            coords[i] = limit1[i]+locs[i]*supp_sizes[i]*0.5;
        }
    }

}

double SGPoint::volume(int dim){
    double vol = 1.0;

    for(int i=0;i<dim;i++){
        vol *= 0.5*supp_sizes[i];
    }
    return vol;
}

bool SGIntegrator::isBoundaryPoint(SGPoint &p){
    bool test = false;
    D3vector levels = p.getLevels();
    for(int i=0;i<this->dim;i++){
        test = test || (levels[i] == 0);
    }
    return test;
}

void SGIntegrator::refine(SGPoint& current){
    D3vector levels  = current.getLevels();
    D3vector  locs = current.getLocs();

    D3vector levels_temp = levels,locs_temp = locs;

    /****Add all the neighbours(2*dim) to the refinement queue****/
    for(int i=0;i<this->dim;i++){
        //If it is a boundary point, add just one neighbour
        if(levels[i] == 0){
            levels_temp[i] = 1;
            locs_temp[i] = 1;
            current.reset(levels_temp,locs_temp,limit1,limit2);
            refine_queue.push(current);

            //Resetting the point to the its original location to prepare for refinement along next direction
            levels_temp[i] = levels[i];
            locs_temp[i] = locs[i];
        }
        //else add both neighbours
        else{
            //Refine
            levels_temp[i]++;

            //Right refine
            locs_temp[i] = 2*locs_temp[i]+1;
            current.reset(levels_temp,locs_temp,limit1,limit2);
            refine_queue.push(current);

            //Left refine
            locs_temp[i] -= 2;
            current.reset(levels_temp,locs_temp,limit1,limit2);
            refine_queue.push(current);

            //Resetting the point to the its original location to prepare for refinement along next direction
            levels_temp[i] = levels[i];
            locs_temp[i] = locs[i];
        }
    }


}

double SGIntegrator::findHierarchicalSurplus(IntegratorPoissonCurved& ipc, SGPoint &p, Curvedfunc func){

    D3vector h = p.getSizes();
    D3vector coords = p.getCoords();
    D3vector levels = p.getLevels();

    //surplus = f(coords)
    double surplus = (ipc.*func)(coords);

    //surplus -= 0.5*(..........)
    for(int i=0;i<this->dim;i++){
        if(levels[i] != 0){
            coords[i] += 0.5*h[i];
            surplus -= 0.5*(ipc.*func)(coords);

            coords[i] -= h[i];
            surplus -= 0.5*(ipc.*func)(coords);
        }
        coords = p.getCoords();
    }

    //surpplus += 0.25*(............)
    int i1=0,i2=0;
    int num_planes = ((this->dim)*(this->dim-1))/2;
    for(int i=0;i<num_planes;i++){
        i1 = i;
        i2 = (i+1)%(this->dim);

        if(levels[i1] != 0 && levels[i2] != 0){
            coords[i1] += 0.5*h[i1];
            coords[i2] += 0.5*h[i2];
            surplus += 0.25*(ipc.*func)(coords);

            coords[i2] -= h[i2];
            surplus += 0.25*(ipc.*func)(coords);

            coords[i1] -= h[i1];
            surplus += 0.25*(ipc.*func)(coords);

            coords[i2] += h[i2];
            surplus += 0.25*(ipc.*func)(coords);
        }

        coords = p.getCoords();
    }


    //surplus -= 0.125*(.................)
    if(this->dim == 3  && levels[0] != 0 && levels[1] != 0 && levels[2] != 0){
           for(int i=0;i<2;i++){
             for(int j=0;j<2;j++){
                for(int k=0;k<2;k++){
                    coords = p.getCoords();
                    coords[0] += 0.5*corn(i)*h[0];
                    coords[1] += 0.5*corn(j)*h[1];
                    coords[2] += 0.5*corn(k)*h[2];

                    surplus -= 0.125*(ipc.*func)(coords);
                }
            }
        }
    }
    //std::cout<<"Point is: "<<std::endl;
    //p.Print();
    //std::cout<<"Surplus is "<<surplus<<std::endl;

    return surplus;

}

//Check levels and dimension of the problem and decide the locations of points to traverse for the subspace integral
void SGIntegrator::decideMinMax(D3vector &min,  D3vector &max, D3vector &levels){
    for(int i=0;i<3;i++){
        if(this->dim < i+1){
            min[i] = max[i] =  0;
        }
        else{
            if(levels[i] == 0){
                min[i] = 0;
                max[i] = 2;
            }
            else{
                min[i] = 1;
                max[i] = pow(2,levels[i])-1;
            }

        }
    }
}

double SGIntegrator::subspaceIntegral(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate,int level_x,int level_y,\
                                      int level_z){
    SGPoint grid_point;
    D3vector levels(level_x,level_y,level_z);
    D3vector locs,boundary_coords;

    double subspace_integral = 0.0;

    D3vector max,min;
    decideMinMax(min,max,levels);

    for(int i = min[0]; i <= max[0] ;i += 2){
        for(int j = min[1]; j <= max[1];j += 2){
            for(int k = min[2]; k <= max[2];k += 2){
                locs[0] = i; locs[1] = j; locs[2] = k;

                grid_point.reset(levels,locs,limit1,limit2);

                subspace_integral += findHierarchicalSurplus(ipc,grid_point,funcToIntegrate)*grid_point.volume(this->dim);
            }
        }
    }
    return subspace_integral;
}

double SGIntegrator::smolyakIntegrate(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate, D3vector &limit1, D3vector &limit2, int max_level){
    double integral = 0.0;

    //Initialize limits of integration
    this->limit1 = limit1;
    this->limit2 = limit2;

    //check the dimension of the problem
    bool dim_2 = (this->dim >= 2);
    bool dim_3 = (this->dim >= 3);
    int max_y,max_z;
    for(int i=0;i<=max_level;i++){
        if(dim_2)
            max_y = max_level-i+this->dim-1;
        else
            max_y = 0;
        for(int j=0;j<=max_y;j++){
            if(dim_3)
                max_z = max_level-i-j+this->dim-1;
            else
                max_z = 0;
            for(int k=0;k<=max_z;k++){
                //Increment the integral using basis functions of the subspace W(i,j,k)
                integral += subspaceIntegral(ipc,funcToIntegrate,i,j,k);
            }
        }
    }

    return integral;
}
double SGIntegrator::integrate(IntegratorPoissonCurved& ipc, Curvedfunc funcToIntegrate,D3vector &limit1, D3vector &limit2){
    double integral = 0.0;
    double loc_integral = 1.0;
    double surplus = 0.0;

    //Initialize the limits of integration
    this->limit1 = limit1;
    this->limit2 = limit2;

    //Create the first points along the corners
    D3vector levels_0(0.0,0.0,0.0);
    D3vector locs;
    SGPoint current;

    bool dim_2 = (this->dim >= 2);
    bool dim_3 = (this->dim >= 3);
    int max_y,max_z;
    //Iterate through all the corners
    for(int i = 0;i <= 2;i += 2){
        max_y = dim_2 ? 2:0;

        for(int j = 0;j <= max_y;j += 2){
            max_z = dim_3 ? 2:0;

            for(int k = 0;k <= max_z; k += 2){
                locs[0] = i; locs[1] = j; locs[2] = k;
                current.reset(levels_0,locs,limit1,limit2);
                refine_queue.push(current);
            }
        }
    }


    while(!refine_queue.empty()){
        current = refine_queue.front();
        //current.Print();
        //If the point hasn't been refined yet, find the local integral
        if(std::find(refined_set.begin(),refined_set.end(),current) == refined_set.end()){
            //local_integral = surplus*(h_x/2)*(h_y/2)*(h_z/2)
            //current.Print();
            surplus = findHierarchicalSurplus(ipc,current,funcToIntegrate);

            loc_integral = surplus;
            loc_integral = surplus*current.volume(this->dim);
            integral += loc_integral;

            //Pop the current point from the refinement queue & push it to the refined set
            refine_queue.pop();
            refined_set.push_back(current);

            //Refine if hierarchical surplus is greater than tolerance
            if(fabs(surplus) >= this->tol  || current.getLevels() == levels_0){
                refine(current);
            }
        }

        //if it has been refined already, just pop
        else{
            refine_queue.pop();
        }
    }

    return integral;
}


#endif // SGINTEGRATOR_H
