#include "matrix_operations.h"



/////////////////////////////////////////////////////////////////
// 1) mapping of matrix indices
/////////////////////////////////////////////////////////////////

// maps matrix indices (i,j) of a bandmatrix of bandwith 5 to an index k of a one dimensional array of size (I*5)
int convert(int i, int j, int I){
    int savei = i;
    if (I <6){return i*I+j;}
    if (j>2&&j < I-2){
        i = i-((j-3)+1);}
    if( j == I-2) i = i-j+3;
    if (j == I-1) i = i-j+4;
    
   if(I>5){ if ( i > 4 || i <0){  cout << "error" << endl;
       cout << " i " << i << " j  " << j << endl;
       cout <<"save i  " << savei << "   I " << I << endl;}}
       
  
        return i*I+j;
    
    
}



///////////////////////////////////////////////////////////////////
// 2) LR Decomposition, forward and backward substitution
///////////////////////////////////////////////////////////////////


void LR(double * L , double* R,int I){
    
    // L[i][i] = 1.0;i.e L is UnitMatrix
    for( int i=0; i< I ; i++){
    
        L[convert(i,i,I)]=1.0;
    }
    
    
   for (int k=0; k< I-1; k++){
        
       for (int i=k+1; i < min(k+3,I); i++){
            
            // L[i][k]=R[i][k]/R[k][k];
            L[convert(i,k,I)]=R[convert(i,k,I)]/R[convert(k,k,I)];
            
            for (int j=k; j < min(k+3,I) ; j++){
                
                //R[i][j]=R[i][j]-L[i][k]*R[k][j];
                R[convert(i,j,I)]=R[convert(i,j,I)]-L[convert(i,k,I)]*R[convert(k,j,I)];
            }
        }
    }
}


    
 void forward(double *y, double * lambda, double *L,int I){

     for (int k=0; k < I-1; k++){
         for (int i=k+1; i < min(k+3,I); i++){

             y[i]=y[i]-L[convert(i,k,I)]*y[k];
         }
     }
 }
 

 
 void backward(double* c, double* y, double* R, int I){
     for(int k = I-1; k>-1;k--){
         c[k]=y[k];
         for (int j=k+1;j < min(k+3,I); j++){
             c[k]=c[k]-R[convert(k,j,I)]*c[j];
         }
         c[k]=c[k]/R[convert(k,k,I)];
     }
 }




