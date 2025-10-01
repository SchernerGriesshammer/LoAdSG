#ifndef STENCIL_H
#define STENCIL_H

#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include "matrix.h"

class Stencil{
private:
    int rows;

    std::vector<std::map<int,double> > data;

public:
    Stencil(int rows_){
        this->rows = rows_;
        std::map<int,double> init;
        init.insert(std::pair<int,double>(0,0.0));
        for(int i=0;i<rows;i++){
            data.push_back(init);
        }
    }

    ~Stencil(){

    }

    int getRows(){
        return rows;
    }
    bool exists(int row, int col){
        std::map<int,double> &rowmap = data[row];
        if(rowmap.find(col) != rowmap.end()){
            return true;
        }
        else return false;
    }

    void create(int row, int col){
        std::map<int,double> &rowmap = data[row];
        rowmap.insert(std::pair<int,double>(col,0.0));
    }

    double& value(int row, int col){
        std::map<int,double> &rowmap = data[row];
        return rowmap[col];
    }

    Matrix multiply(Matrix &b){
        Matrix c(b);
        if(b.getCols() != 1 || b.getRows() != rows){
            std::cerr<<"Stencil vector multiplication not possible";
        }
        for(int i=0;i<rows;i++){
            c.value(i,0) = 0;
            for(std::map<int,double>::iterator it = data[i].begin();it != data[i].end(); it++){
                c.value(i,0) +=  b.value(it->first,0)*(it->second);
            }
        }
        return c;
    }
};
#endif // STENCIL_H

