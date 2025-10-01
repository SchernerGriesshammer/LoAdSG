// ------------------------------------------------------------
// main.cc
//
// ------------------------------------------------------------

#include <iostream>
#include <fstream>

//test NOW
using namespace std;

#include "../interfaceMatrices.h"
#include "../constantIntegrators.h"
#include "../interatorBasisFunction.h"
#include "curvedPoissonIntegrator.h"
#include "SGIntegrator.h"
#include "stencil.h"


//test NOW
//using namespace std;

#define Dim 3

D3vector transform_cyl_NE_inner ( double t) {
  return D3vector ( cos ( t*0.5*M_PI )- ( 1-t ), sin ( t*0.5*M_PI )-t, 0.0 );
}

D3vector transform_cyl_NE ( double t) {
  return 2.0 * D3vector ( cos ( t*0.5*M_PI )- ( 1-t ), sin ( t*0.5*M_PI )-t, 0.0 );
}

/*D3vector transform_cyl_NE_inner ( double t) {
  return D3vector ( cos ( t*0.5*M_PI ), sin ( t*0.5*M_PI ), 0.0 );
}

D3vector transform_cyl_NE ( double t) {
  return 2.0 * D3vector ( cos ( t*0.5*M_PI ), sin ( t*0.5*M_PI ), 0.0 );
}*/

D3vector transform_const ( double t) {
  return D3vector(0.0,0.0,0.0);
}


//Transformation functions and their derivatives
D3vector transform_SD(double t){
    return D3vector(1.0+t,0.0,0.0);
}

D3vector transDerivs_SD(double t){
    return D3vector(1.0,0.0,0.0);
}

D3vector transform_ND(double t){
    return D3vector(0.0,1.0+t,0.0);
}

D3vector transDerivs_ND(double t){
    return D3vector(0.0,1.0,0.0);
}

D3vector transform_ST(double t){
    return D3vector(1.0+t,0.0,0.5);
}

D3vector transDerivs_ST(double t){
    return D3vector(1.0,0.0,0.0);
}

D3vector transform_NT(double t){
    return D3vector(0.0,1.0+t,0.5);
}

D3vector transDerivs_NT(double t){
    return D3vector(0.0,1.0,0.0);
}


D3vector transform_WD(double t){
    return D3vector(cos(t*0.5*M_PI), sin(t*0.5*M_PI), 0.0 );
}

D3vector transDerivs_WD(double t){
    return D3vector(-0.5*M_PI*sin(t*0.5*M_PI),0.5*M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_ED(double t){
    return D3vector(2.0*cos(t*0.5*M_PI), 2.0*sin(t*0.5*M_PI), 0.0 );
}

D3vector transDerivs_ED(double t){
    return D3vector(-1.0*M_PI*sin(t*0.5*M_PI),M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_WT(double t){
    return D3vector(cos(t*0.5*M_PI), sin(t*0.5*M_PI), 0.5 );
}

D3vector transDerivs_WT(double t){
    return D3vector(-0.5*M_PI*sin(t*0.5*M_PI),0.5*M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_ET(double t){
    return D3vector(2.0*cos(t*0.5*M_PI), 2.0*sin(t*0.5*M_PI), 0.5);
}

D3vector transDerivs_ET(double t){
    return D3vector(-1.0*M_PI*sin(t*0.5*M_PI),M_PI*cos(t*0.5*M_PI),0.0);
}


D3vector transform_SW(double t){
    return D3vector(1.0,0.0,0.5*t);
}

D3vector transDerivs_SW(double t){
    return D3vector(0.0,0.0,0.5);
}


D3vector transform_SE(double t){
    return D3vector(2.0,0.0,0.5*t);
}

D3vector transDerivs_SE(double t){
    return D3vector(0.0,0.0,0.5);
}


D3vector transform_NW(double t){
    return D3vector(0.0,1.0,0.5*t);
}

D3vector transDerivs_NW(double t){
    return D3vector(0.0,0.0,0.5);
}

D3vector transform_NE(double t){
    return D3vector(0.0,2.0,0.5*t);
}

D3vector transDerivs_NE(double t){
    return D3vector(0.0,0.0,0.5);
}

bool isBoundaryNode(int i, int j, int k, int N){
    if(i == 0 || i == N || j == 0 || j == N || k == 0 || k == N){
        return true;
    }
    else return false;
}

//return the global index of a node, given the indices of the element and the type of basis function
int decideGlobalIndex(int i, int j, int k, int N, IteratorBasisFunction<3> U){

    int ind[3] = {0,0,0};
    for(int d = 0; d<3; d++){
     //   std::cout<<U.getBasisTypeNum(d)<<"  ";
        if(U.getBasisTypeCoord(d) == rightBasis){
            ind[d]++;
        }
    }

    //Check if it is a boundary node. If yes return -1
    if(isBoundaryNode(i+ind[0],j+ind[1],k+ind[2],N)){
        return -1;
    }

    //Transform to the corresponding indices in the (N-1)*(N-1)*(N-1) grid of inner points
    i  = i+ind[0]-1;j = j+ind[1]-1;k = k+ind[2]-1;
    return(k*(N-1)*(N-1)+j*(N-1)+i);
}

//Perform assembly for the element (i,j,k)
void addComponents(Matrix& A, Matrix& b, int i, int j, int k, int N, IntegratorPoissonCurved& ipc){
    BasisFunctionType u[3];
    BasisFunctionType v[3];

    double dh = 1.0/(double)N;
    Curvedfunc uAv = &IntegratorPoissonCurved::uAv;
    Curvedfunc fv = &IntegratorPoissonCurved::fv;

    SGIntegrator IG(20,0.0001);

    //Limits of integration for the element
    D3vector lim1(i*dh,j*dh,k*dh),lim2((i+1)*dh,(j+1)*dh,(k+1)*dh);
    ipc.setLimits(lim1,lim2);

    int row,col;//indices of the global matrix to which to add to
    for(IteratorBasisFunction<Dim> iterV;iterV.hasNext();iterV.next()) {
          row = decideGlobalIndex(i,j,k,N,iterV);
          if(row == -1){
              continue;
          }

          //Determine type of v
          for(int d=0;d<Dim;++d) {
            v[d] =iterV.getBasisTypeCoord(d);
          }

          //Set v
          ipc.setState(v,false);

          b.value(row,0) += IG.smolyakIntegrate(ipc,fv,lim1,lim2,10);//Add the contribution to RHS
          for(IteratorBasisFunction<Dim> iterU;iterU.hasNext();iterU.next()){
              col =  decideGlobalIndex(i,j,k,N,iterU);
              if(col == -1){
                  continue;
              }

              //Determine type of u
              for(int d=0;d<Dim;++d){
                u[d] =iterU.getBasisTypeCoord(d);
              }

              //Set u
              ipc.setState(u,true);

              //Add contribution to the global matrix

              A.value(row,col) +=  IG.smolyakIntegrate(ipc, uAv, lim1, lim2, 10);

          }
    }
}

void addComponentsSparse(Stencil& A, Matrix& b, int i, int j, int k, int N, IntegratorPoissonCurved& ipc){
    BasisFunctionType u[3];
    BasisFunctionType v[3];

    double dh = 1.0/(double)N;
    Curvedfunc uAv = &IntegratorPoissonCurved::uAv;
    Curvedfunc fv = &IntegratorPoissonCurved::fv;

    SGIntegrator IG(20,0.0001);

    //Limits of integration for the element
    D3vector lim1(i*dh,j*dh,k*dh),lim2((i+1)*dh,(j+1)*dh,(k+1)*dh);
    ipc.setLimits(lim1,lim2);

    int row,col;//indices of the global matrix to which to add to
    double entry = 0.0;
    for(IteratorBasisFunction<Dim> iterV;iterV.hasNext();iterV.next()) {
          row = decideGlobalIndex(i,j,k,N,iterV);
          if(row == -1){
              continue;
          }

          //Determine type of v
          for(int d=0;d<Dim;++d) {
            v[d] =iterV.getBasisTypeCoord(d);
          }

          //Set v
          ipc.setState(v,false);

          b.value(row,0) += IG.smolyakIntegrate(ipc,fv,lim1,lim2,14);//Add the contribution to RHS
          for(IteratorBasisFunction<Dim> iterU;iterU.hasNext();iterU.next()){
              col =  decideGlobalIndex(i,j,k,N,iterU);
              if(col == -1){
                  continue;
              }

              //Determine type of u
              for(int d=0;d<Dim;++d){
                u[d] =iterU.getBasisTypeCoord(d);
              }

              //Set u
              ipc.setState(u,true);

              //Add contribution to the global matrix

              entry =  IG.smolyakIntegrate(ipc, uAv, lim1, lim2, 10);
              if(entry != 0.0){
                  if(A.exists(row,col)){
                      A.value(row,col) += entry;
                  }
                  else{
                      A.create(row,col);
                      A.value(row,col) +=  entry;
                  }
              }
          }
    }
}

int bicgStab(Stencil &A, Matrix &b, Matrix &x, double tol, int max_itr){
    int n = b.getRows();
    double rho,rho_old,alpha,beta,w;
    Matrix xi(n,1),ri(n,1),r_cap(n,1),p(n,1),v(n,1),s(n,1),t(n,1);

    ri = b-A.multiply(xi);
    r_cap = ri;
    int i = 0;

    while(ri.norm2() >= tol*b.norm2() && i <= max_itr){
        rho = r_cap.dp(ri);

        if(rho == 0.0){
            std::cerr<<"BiCGSTAB  breakdown"<<std::endl;
            return 0;
        }

        if(i == 0){
            p = ri;
        }
        else{
            beta = (rho*alpha)/(rho_old*w);
            p =  ri+(p-v.sc(w)).sc(beta);
        }

        v = A.multiply(p);

        alpha = rho/(r_cap.dp(v));

        s = ri-v.sc(alpha);

        if(s.norm2() <= tol*b.norm2()){
            xi = xi+p.sc(alpha);
            x = xi;
            ri = s;

            std::cout<<"BICGSTAB:Iterations run = "<<i<<std::endl;
            return 1;
        }

        t = A.multiply(s);
        w = (t.dp(s))/(t.dp(t));

        xi  = xi+p.sc(alpha)+s.sc(w);

        ri =  s-t.sc(w);

        rho_old = rho;
        i++;

        std::cout<<i<<std::endl;
    }

    std::cout<<"BICGSTAB:Maximum iterations reached"<<std::endl;
    x = xi;

    return 1;
}
int main(int argc, char** argv) {
  cout.precision(10);
  cout.setf(std::ios::fixed,std::ios::floatfield);

  
  // Example: POissin Curved
  //////////////////////////////////
  D3vector corners[8];
  D3vector (*transformEdge[12])(double);
  D3vector (*transDerivs[12])(double);

  // this describes a small ring segment of radii r, R and length L:
  double r = 1.0;
  double R = 2.0;
  double L = 0.5;
    
  corners[WSDdir3D] = D3vector ( r,0.0,0.0 );
  corners[ESDdir3D] = D3vector ( R,0.0,0.0 );
  corners[WNDdir3D] = D3vector ( 0.0,r,0.0 );
  corners[ENDdir3D] = D3vector ( 0.0,R,0.0 );
  corners[WSTdir3D] = D3vector ( r,0.0,L );
  corners[ESTdir3D] = D3vector ( R,0.0,L );
  corners[WNTdir3D] = D3vector ( 0.0,r,L );
  corners[ENTdir3D] = D3vector ( 0.0,R,L );

  /*transformEdge[SDed] = transform_const;
  transformEdge[NDed] = transform_const;
  transformEdge[STed] = transform_const;
  transformEdge[NTed] = transform_const;
  transformEdge[WDed] = transform_cyl_NE_inner;
  transformEdge[EDed] = transform_cyl_NE;
  transformEdge[WTed] = transform_cyl_NE_inner;
  transformEdge[ETed] = transform_cyl_NE;
  transformEdge[SWed] = transform_const;
  transformEdge[SEed] = transform_const;
  transformEdge[NWed] = transform_const;
  transformEdge[NEed] = transform_const;*/

  transformEdge[SDed] = transform_SD;
  transformEdge[NDed] = transform_ND;
  transformEdge[STed] = transform_ST;
  transformEdge[NTed] = transform_NT;
  transformEdge[WDed] = transform_WD;
  transformEdge[EDed] = transform_ED;
  transformEdge[WTed] = transform_WT;
  transformEdge[ETed] = transform_ET;
  transformEdge[SWed] = transform_SW;
  transformEdge[SEed] = transform_SE;
  transformEdge[NWed] = transform_NW;
  transformEdge[NEed] = transform_NE;

  transDerivs[SDed] = transDerivs_SD;
  transDerivs[NDed] = transDerivs_ND;
  transDerivs[STed] = transDerivs_ST;
  transDerivs[NTed] = transDerivs_NT;
  transDerivs[WDed] = transDerivs_WD;
  transDerivs[EDed] = transDerivs_ED;
  transDerivs[WTed] = transDerivs_WT;
  transDerivs[ETed] = transDerivs_ET;
  transDerivs[SWed] = transDerivs_SW;
  transDerivs[SEed] = transDerivs_SE;
  transDerivs[NWed] = transDerivs_NW;
  transDerivs[NEed] = transDerivs_NE;

   
  IntegratorPoissonCurved ipc(transformEdge,transDerivs,corners);

  //This code snippet is to check the correctness of the mapping function
  /*D3vector coords(0,0,0);
  for(int i=0;i<26;i++){
      for(int j=0;j<26;j++){
          for(int k=0;k<26;k++){
              coords = ipc.Map2(i/25.0,j/25.0,k/25.0);
              std::cout<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<std::endl;
          }
      }
  }*/

  //Number of elements = N*N*N
  int N = atoi(argv[1]);

  double dh = 1.0/(double)N;

  //Tolerance for the Krylov solver
  double tol =  0.000001;
  int max_itr = 1000;

  //Number of unknowns
  int n = (N-1)*(N-1)*(N-1);

  Matrix b(n,1), u(n,1), uactual(n,1), error(n,1);//Global stiffness matrix and the RHS vector

  Stencil As(n);

  //Iterate through all the elements and perform assembly
   for(int i=0;i<N;i++){
       for(int j=0;j<N;j++){
           for(int k=0;k<N;k++){
           //   std::cout<<i<<" "<<j<<" "<<k<<std::endl;
              addComponentsSparse(As,b,i,j,k,N,ipc);

           }
       }
   }

   D3vector x_trans(0,0,0),x(0,0,0);
   for(int k = 0;k<N-1;k++){
       for(int j=0;j<N-1;j++){
           for(int i=0;i<N-1;i++){
            x[0] = (i+1)*dh;
            x[1] = (j+1)*dh;
            x[2] = (k+1)*dh;

            x_trans = ipc.Map2(x[0],x[1],x[2]);

            uactual.value(k*(N-1)*(N-1)+j*(N-1)+i,0) = x_trans[0]*x_trans[1]*x_trans[2]*(x_trans[2]-0.5)*cos(0.5*M_PI*(x_trans[0]*x_trans[0]+x_trans[1]*x_trans[1]))*sin(0.5*M_PI*(x_trans[0]*x_trans[0]+x_trans[1]*x_trans[1]));
           }
       }
   }

  //Solve the linear system using bicgstab
  bicgStab(As,b,u,tol,max_itr);

  //Compute average error
  double avg_err = 0.0;
  for(int i=0;i<n;i++){
      error.value(i,0) = fabs(u.value(i,0)-uactual.value(i,0));
      avg_err += error.value(i,0);
  }
  avg_err = avg_err/n;

  std::cout<<std::endl;
  std::cout<<"Average error in the solution is "<<avg_err<<std::endl;
  //Export the global stiffness matrix and RHS vector to separate files
  //A.exportToFile("globalStiffnessMatrix.txt");
  //b.exportToFile("rhs.txt");
}


