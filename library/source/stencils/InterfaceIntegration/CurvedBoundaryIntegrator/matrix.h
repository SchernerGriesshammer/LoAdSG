#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>

class Matrix{
private:
    int rows,cols;

    std::vector<double> data;

public:
    Matrix(int dim1, int dim2);

    Matrix(const Matrix& M);

    ~Matrix(){
    }
    Matrix& operator=(const Matrix&  M);

    double operator[](int index) const;

    double& operator[](int index);

    Matrix operator*(const Matrix& M);

    Matrix subMatrix(int i, int j);

    double det();

    Matrix invert();

    double value(int i1, int i2) const{
        return data[i1*cols+i2];
    }

    double& value(int i1, int i2){
        return data[i1*cols+i2];
    }

    int getRows() const{
        return rows;
    }

    int getCols() const{
        return cols;
    }

    void scale(double s){
        for(int i=0;i<rows*cols;i++){
            data[i] = s*data[i];
        }
    }

    void Print(){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                std::cout<<value(i,j)<<"\t";
            }
            std::cout<<std::endl;
        }
    }

    Matrix transpose();

    std::vector<double> getData() const{
        return data;
    }

    void exportToFile(const char *filename);

    double norm2(){
        double vecnorm = 0.0;
        for(int i=0;i<rows;i++){
            for(int  j=0;j<cols;j++){
                vecnorm +=  value(i,j)*value(i,j);
            }
        }
        vecnorm  = vecnorm/(rows*cols);
        return(std::sqrt(vecnorm));
    }

    Matrix operator+(const Matrix& m){
        Matrix newmat(*this);
        for(int i=0;i<rows;i++){
            for(int  j=0;j<cols;j++){
                newmat.value(i,j) += m.value(i,j);
            }
        }

        return newmat;
    }

    Matrix operator-(const Matrix& m){
        Matrix newmat(*this);
        for(int i=0;i<rows;i++){
            for(int  j=0;j<cols;j++){
                newmat.value(i,j) -= m.value(i,j);
            }
        }
        return newmat;
    }

    double dp(const Matrix& m){
        double res =  0.0;
        for(int i=0;i<rows;i++){
            res += value(i,0)*m.value(i,0);
        }
        return res;
    }

    Matrix sc(double val){
        Matrix newmat(*this);
        for(int i=0;i<rows;i++){
            for(int  j=0;j<cols;j++){
                newmat.value(i,j) = val*newmat.value(i,j);
            }
        }
        return newmat;
    }

};


Matrix::Matrix(int dim1_, int dim2_){
   rows = dim1_;
   cols = dim2_;

   data.resize(rows*cols);

   std::fill(data.begin(),data.end(),0.0);
}

Matrix::Matrix(const Matrix &M){

    rows = M.getRows();
    cols = M.getCols();

    data = M.getData();
}

Matrix& Matrix::operator=(const Matrix &M){

    rows = M.getRows();
    cols = M.getCols();

    data = M.getData();

    return *this;
}

double& Matrix::operator[](int index){
    return data[index];
}

double Matrix::operator[](int index) const{
    return data[index];
}


Matrix Matrix::operator*(const Matrix& M){
    if(this->cols != M.getRows()){
        std::cerr<<"Dimensions incompatible for matrix multiplication"<<std::endl;
        exit(-1);
    }

    int r = this->rows;
    int c = M.getCols();

    Matrix result(r,c);


    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            for(int k=0;k<this->cols;k++){
                result.value(i,j) += value(i,k)*M.value(k,j);
            }
        }
    }

    return result;
}

Matrix Matrix::subMatrix(int r, int c){
    Matrix sub(rows-1,cols-1);
    int temp_i,temp_j;

    for(int i=0;i<rows-1;i++){
        for(int j=0;j<cols-1;j++){
            temp_i = i;
            temp_j = j;

            if(i >= r){
                temp_i = i+1;
            }
            if(j >= c){
                temp_j = j+1;
            }
            sub.value(i,j) = value(temp_i,temp_j);
        }
    }

    return sub;
}

Matrix Matrix::transpose(){
    Matrix  trans(cols,rows);

    for(int i=0;i<cols;i++){
        for(int j=0;j<rows;j++){
            trans.value(i,j) = this->value(j,i);
        }
    }

    return trans;
}

Matrix Matrix::invert(){
    if(det() == 0){
        std::cerr<<"Matrix not invertible"<<std::endl;
        exit(-1);
    }
    Matrix inv(rows,cols);
    for(int i=0;i<rows;i++){
        for(int j=0; j<cols;j++){
            Matrix temp = subMatrix(i,j);
            inv.value(i,j) = temp.det()*pow(-1,i+j);
        }
    }

    inv = inv.transpose();
    inv.scale(1.0/(this->det()));

    return inv;
}
double Matrix::det(){
    double det = 0.0;
    if(rows != cols){
        std::cerr<<"Not a square matrix. Hence can't compute the determinant."<<std::endl;
        exit(-1);
    }

    if(rows == 1){
        return data[0];
    }

    if(rows == 2){
        return (value(0,0)*value(1,1)-value(0,1)*value(1,0));
    }

    if(rows == 3){
        for(int i=0;i<3;i++){
            Matrix temp(this->subMatrix(0,i));

            det += pow(-1,i)*value(0,i)*temp.det();
        }
        return det;
    }
    return 0.0;
}

void Matrix::exportToFile(const char* filename){
    std::ofstream outfile;
    outfile.open(filename);

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            outfile<<value(i,j)<<std::endl;
        }
    }

    outfile.close();
}
#endif // MATRIX_H
