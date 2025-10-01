/**********************************************************************************
 * Copyright 2016 Christoph Pflaum, Parikshit Upadhyaya
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 **********************************************************************************/

#ifndef INTERATORPOISSONCURVED_H
#define INTERATORPOISSONCURVED_H

#include "../math_lib.h"
#include "../abbrevi_common.h"
#include "matrix.h"

/**
 * integrator of Poisson equation in 3D in curved domain
 * Parikshit Upadhyaya
 */
class IntegratorPoissonCurved : public InterfaceLocalStiffnessMatrices<3> {
   public:
     /**
      * @param transformEdge[k] is a mapping from [0,1] -> D3vector     for every k=0,...,11
      *                         it must be transformEdge[k](0) = (0,0,0) = transformEdge[k](1)
      * @param corner contains the coordinates of the eight corners of the deformed cube
      */
      IntegratorPoissonCurved(D3vector (*transformEdge_[12])(double), D3vector (*transDerivs_[12])(double),D3vector corner_[8]){
          //Transformation functions for the edges
          for(int i=0;i<12;++i){
              transformEdge[i] = transformEdge_[i];
              transDerivs[i] = transDerivs_[i];
          }

          //Corners of the curvilinear block
          for(int i=0;i<8;++i){
              corner[i] = corner_[i];
          }
      }

      //The mapping function that maps a unit cube to the curvilinear block.....(Old version)
      D3vector Map1(double eta, double xi, double phi) const;

      //The mapping function that maps a unit cube to the curvilinear block.....(New version)
      D3vector Map2(double s, double t, double u) const;

      //(Bit interpolation )Return x when bit = 1, else return 1-x
      double bitInterp(int bit,double  x) const;

      //Return the index of the edge, given the dimension index and 2 bits
      int edIndex(int dim_ind, int b1, int b2) const;

      //Return the index of the corner given the 3 bits
      int coIndex(int b1, int b2, int b3) const;

      //If dim_ref == dim, then call bitSign, else call bitInterp
      double interpOrSign(int b,int dim, int dim_ref,D3vector& x);

      //Return -1.0 or 1.0 depending on sign
      double bitSign(int b);

      //Needed inside MapDeriv() for computation of edge indices
      int edIndexDeriv(int edge_ref, int b1, int b2, int dim_ref);

      //Functions that describe the problem
      double A(D3vector& x_,int i, int j);

      double b(D3vector& x_,int i);

      double c(D3vector& x_);

      //RHS for the solution u(x,y,z) = xyz(z-0.5)*cos(0.5*pi*(x²+y²))*sin(0.5*(x²+y²))
      double f1(D3vector& x_);

      //evaluate A at x and store it in mat
      void evalA(Matrix& mat, D3vector& x);

      //Evaluate the transpose of the jacobian at x and store in J
      void evalJ(Matrix& J, D3vector& x);

      //Evaluate the basis function u at x
      double evalU(D3vector& x, BasisFunctionType *state);

      //Evaluate the gradient of u at x
      void evaldel(Matrix& delu, D3vector& x, BasisFunctionType *state);

      //Used to compute one column of the transpose of jacobian. eg: dim = 0 means  d_phi/d_s at x
      D3vector MapDeriv(int dim, D3vector& x);

      //Compute the LHS of the integration to fill in the entries of the local stiffness matrix
      double stencil_integration(double p_left[3], double p_right[3], BasisFunctionType u[3], BasisFunctionType v[3]);

      //Functions that describe the the variational formulation. Will go as input to the sparse grid integrator
      double uAv(D3vector& x);

      double ubv(D3vector& x_);

      double cuv(D3vector& x_);

      double fv(D3vector& x_);

      //if U = true, modify state_u, else modify state_v
      void setState(BasisFunctionType *state, bool U);

      //set the limits of integration
      void  setLimits(D3vector &lim1_, D3vector& lim2_){
          lim1 = lim1_;
          lim2 = lim2_;
      }

   private:   
      D3vector (*transformEdge[12])(double);
      D3vector (*transDerivs[12])(double);
      D3vector corner[8];

      D3vector lim1;
      D3vector lim2;

      BasisFunctionType state_u[3];
      BasisFunctionType state_v[3];

};

//Functions that describe the problem
double IntegratorPoissonCurved::A(D3vector& x_,int i, int j){
    if(i == j){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

double IntegratorPoissonCurved::b(D3vector& x_,int i){
    return 0.0;
}

double IntegratorPoissonCurved::c(D3vector& x_){
    return 0.0;
}

//RHS for the solution u(x,y,z) = xyz(z-0.5)*cos(0.5*pi*(x²+y²))*sin(0.5*(x²+y²))
double IntegratorPoissonCurved::f1(D3vector& x_){
    double x = x_[0], y = x_[1], z = x_[2];

    double pi = M_PI;

    double sin_b = sin(0.5*pi*(x*x+y*y));
    double cos_b = cos(0.5*pi*(x*x+y*y));
    double z_b = (z-0.5);

    double result = 0.0;

    result += 6.0*pi*x*y*z*z_b*sin_b*sin_b - \
              6.0*pi*x*y*z*z_b*cos_b*cos_b - \
              2.0*x*y*sin_b*cos_b + \
              4.0*pi*pi*x*y*y*y*z*z_b*sin_b*cos_b +\
              4.0*pi*pi*x*x*x*y*z*z_b*sin_b*cos_b;

    return result;
}

//Functions that describe the the variational formulation. Will go as input to the sparse grid integrator
double IntegratorPoissonCurved::uAv(D3vector& x){
    D3vector x_new = Map2(x[0],x[1],x[2]);

    Matrix AoPhi(3,3),A_star(3,3),Jt(3,3),J_in(3,3),Jt_in(3,3),delU(1,3),delV(1,3),delVt(3,1);
    double dJ;

    evalA(AoPhi,x_new);
    evalJ(Jt,x);

    //dJ = fabs(Jt.det());
    dJ = Jt.det();
    if(dJ == 0.0 ){
        std::cerr<<"Non invertible Jacobian encountered"<<std::endl;
        return 0.0;
    }

    Jt_in = Jt.invert();
    J_in = Jt_in.transpose();

    A_star = J_in*AoPhi*Jt_in;
    A_star.scale(dJ);

    evaldel(delU,x,state_u);
    evaldel(delV,x,state_v);
    delVt = delV.transpose();

    double result = (delU*A_star*delVt).value(0,0);
    return result;
    //return 2.0;
}

double IntegratorPoissonCurved::ubv(D3vector& x_){
    return 0.0;
}

double IntegratorPoissonCurved::cuv(D3vector& x_){
    return 0.0;
}

double IntegratorPoissonCurved::fv(D3vector& x_){
    D3vector x_new = Map2(x_[0],x_[1],x_[2]);

    Matrix Jt(3,3);
    evalJ(Jt,x_);

    //double result = f1(x_new)*fabs(Jt.det())*evalU(x_,state_v);
    double result = f1(x_new)*Jt.det()*evalU(x_,state_v);

    return result;
}

double IntegratorPoissonCurved::bitInterp(int bit,  double x) const{
    if(bit == 0){
        return 1-x;
    }
    else if(bit == 1){
        return x;
    }
    else{
        std::cerr<<"Incorrect bit value. Can be either 0 or 1"<<std::endl;
        return 0;
    }
}

int IntegratorPoissonCurved::edIndex(int dim_ind, int b1, int b2) const{
  int offset = b2*2+b1;

  return (dim_ind*4+offset);

}

int IntegratorPoissonCurved::coIndex(int b1, int b2, int b3) const{
   return (4*b3+2*b2+b1);
}


double IntegratorPoissonCurved::interpOrSign(int b,int dim, int dim_ref,D3vector& x){
    if(dim == dim_ref){
        return bitSign(b);
    }
    else{
        return bitInterp(b,x[dim]);
    }
}

double IntegratorPoissonCurved::bitSign(int b){
    double ret = (b == 0)?-1.0:1.0;
    return ret;
}


int IntegratorPoissonCurved::edIndexDeriv(int edge_ref, int b1, int b2, int dim_ref){
    if(dim_ref == 0){
        return edIndex(edge_ref/4,b1,b2);
    }
    else if(dim_ref == 1 && edge_ref  == 0){
        return edIndex(edge_ref/4,b1,b2);
    }
    else if(dim_ref == 2 && edge_ref == 8){
        return edIndex(edge_ref/4,b1,b2);
    }
    else{
        return edIndex(edge_ref/4,b2,b1);
    }
}


D3vector IntegratorPoissonCurved::Map2(double s, double t, double u) const{
    D3vector ret(0.0,0.0,0.0);

    //Interpolation along the edges of the domain
    for(int b1 = 0; b1 <=1; b1 ++){
        for(int b2 = 0; b2 <= 1; b2 ++){
            ret  = ret + bitInterp(b1,t)*bitInterp(b2,u)*transformEdge[edIndex(0,b1,b2)](s);
            ret  = ret + bitInterp(b1,s)*bitInterp(b2,u)*transformEdge[edIndex(1,b1,b2)](t);
            ret  = ret + bitInterp(b1,s)*bitInterp(b2,t)*transformEdge[edIndex(2,b1,b2)](u);
        }
    }

    //Interpolation along the corners of the domain
    for(int b1 = 0; b1 <=1; b1 ++){
        for(int b2 = 0; b2 <= 1; b2 ++){
            for(int b3 = 0; b3 <= 1; b3 ++){
                ret = ret - 2*bitInterp(b1,s)*bitInterp(b2,t)*bitInterp(b3,u)*corner[coIndex(b1,b2,b3)];
            }
        }
    }

    return ret;
}

void IntegratorPoissonCurved::evalA(Matrix& mat, D3vector& x){
      for(int i = 0; i<3; i++){
          for(int j = 0; j<3; j++){
              mat.value(i,j) = A(x,i,j);
          }
      }
}

D3vector IntegratorPoissonCurved::MapDeriv(int dim,D3vector& x){

    D3vector ret(0.0,0.0,0.0);

    double s1 = x[dim], s2 = x[(dim+1)%3], s3 = x[(dim+2)%3];
    int ed1 = dim*4, ed2 = ((dim+1)%3)*4, ed3 = ((dim+2)%3)*4;

    //Terms related to edge interpolation
    for(int b1=0; b1<=1; b1++){
        for(int b2=0; b2<=1; b2++){
            ret = ret + bitInterp(b1,s2)*bitInterp(b2,s3)*transDerivs[edIndexDeriv(ed1,b1,b2,dim)](s1);

            ret = ret + bitSign(b1)*bitInterp(b2,s3)*transformEdge[edIndexDeriv(ed2,b1,b2,dim)](s2);

            ret = ret + bitSign(b1)*bitInterp(b2,s2)*transformEdge[edIndexDeriv(ed3,b1,b2,dim)](s3);
        }
    }

    //Terms related to corner interpolation
    for(int b1 = 0; b1 <=1; b1 ++){
        for(int b2 = 0; b2 <= 1; b2 ++){
            for(int b3 = 0; b3 <= 1; b3 ++){
                ret = ret - 2*interpOrSign(b1,0,dim,x)*interpOrSign(b2,1,dim,x)*interpOrSign(b3,2,dim,x)*corner[coIndex(b1,b2,b3)];
            }
        }
    }

    return ret;
}

void IntegratorPoissonCurved::evalJ(Matrix& J, D3vector& x){
      D3vector cols[3];
      cols[0] = MapDeriv(0,x);
      cols[1] = MapDeriv(1,x);
      cols[2] = MapDeriv(2,x);

      for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
              J.value(i,j) = (cols[i])[j];
          }
      }
}

double IntegratorPoissonCurved::evalU(D3vector& x, BasisFunctionType *state){
    double  result = 1.0;
    for(int i=0;i<3;i++){
        if(state[i] == leftBasis){
            result *= (lim2[i]-x[i])/(lim2[i]-lim1[i]);
        }
        else{
            result *= (x[i]-lim1[i])/(lim2[i]-lim1[i]);
        }
    }

    return result;
}

//Evaluate the gradient of u at x
void IntegratorPoissonCurved::evaldel(Matrix& delu, D3vector& x, BasisFunctionType *state){
    double prod[3];
    int ind1,ind2;

    for(int i=0;i<3;i++){
        if(state[i] == leftBasis){
            delu.value(0,i) = -1.0/(lim2[i]-lim1[i]);
            prod[i] = (lim2[i]-x[i])/(lim2[i]-lim1[i]);
        }
        else{
            delu.value(0,i) = 1.0/(lim2[i]-lim1[i]);
            prod[i] = (x[i]-lim1[i])/(lim2[i]-lim1[i]);
        }
    }

    for(int i=0;i<3;i++){
        ind1 = (i+1)%3;
        ind2 = (i+2)%3;
        delu.value(0,i) *= prod[ind1]*prod[ind2];
    }
}

double IntegratorPoissonCurved::stencil_integration(double p_left[3], double p_right[3],
                  BasisFunctionType u[3], BasisFunctionType v[3]) {

     lim1 = D3vector(p_left[0],p_left[1],p_left[2]);
     lim2 = D3vector(p_right[0],p_right[1],p_right[2]);

     //Modifying state_u and state_v lets the functions uAv(), ubv() and cuv() know which types of u and v we are dealing with
     for(int i=0;i<3;i++){
         state_u[i] = u[i];
         state_v[i] = v[i];
     }

     //Divide the variational formulation into 3 parts
     //double integral1 = IG.smolyakIntegrate(*this,this->uAv,lim1,lim2,10);//Integrate grad(u_tilda)*A*grad(v_tilda)
    // double integral2 = IG.smolyakIntegrate(&this->ubv,lim1,lim2,10);//Integrate <grad(v_tilda)*b*>v_tilda
    // double integral3 = IG.smolyakIntegrate(&this->cuv,lim1,lim2,10);//Integrate c*u_tilda*v_tilda

     return 0;
}

void IntegratorPoissonCurved::setState(BasisFunctionType *state, bool U){
    if(U){
        for(int i=0;i<3;i++){
            state_u[i] = state[i];
        }
    }
    else{
        for(int i=0;i<3;i++){
            state_v[i] = state[i];
        }
    }
}

D3vector IntegratorPoissonCurved::Map1(double eta, double xi, double phi) const {
  D3vector Psi[12];
  D3vector PL, PR;
  D3vector PSW, PSE, PNW, PNE;
  D3vector Pres, PT, PD;
  double   p_EW, p_NS;
  double   p;

  for ( int ed = 0; ed < 12; ++ed )  {
    p = find_p ( ( Edges_cell ) ed,eta,xi,phi );  //p is either eta, phi or xi
    Psi[ed] = transformEdge[ed]( p );
    PL = corner[ Transform ( ( Edges_cell ) ed,Ldir1D )  ];  // has to be implemeted to give corner[...]
    PR = corner[ Transform ( ( Edges_cell ) ed,Rdir1D )  ];
    Psi[ed] = Psi[ed] + PL + ( PR-PL ) * p;
  }

  Pres = D3vector ( 0.0,0.0,0.0 );
  for ( int md = 0; md < 3; ++md )
  {
    if ( md==0 ) // EW
    {
      PSW = Psi[SDed];
      PSE = Psi[NDed];
      PNW = Psi[STed];
      PNE = Psi[NTed];

      p_EW = xi;
      p_NS = phi;
    }
    if ( md==1 ) // NS
    {
      PSW = Psi[WDed];
      PSE = Psi[EDed];
      PNW = Psi[WTed];
      PNE = Psi[ETed];

      p_EW = eta;
      p_NS = phi;
    }
    if ( md==2 ) // TD
    {
      PSW = Psi[SWed];
      PSE = Psi[SEed];
      PNW = Psi[NWed];
      PNE = Psi[NEed];

      p_EW = eta;
      p_NS = xi;
    }
    Pres = Pres + PSW +
           ( PSE - PSW ) * p_EW +
           ( PNW - PSW ) * p_NS +
           ( PNE - PSE - PNW + PSW ) * p_EW * p_NS;
  }

  PD = corner[WSDdir3D]  +
       ( corner[ESDdir3D] - corner[WSDdir3D] ) * eta +
       ( corner[WNDdir3D] - corner[WSDdir3D]) * xi +
       ( corner[ENDdir3D] - corner[ESDdir3D] - corner[WNDdir3D] + corner[WSDdir3D]) * eta * xi;//remove all Hs and add corner[...dir3D] to compile

  PT = corner[WSTdir3D] +
       ( corner[ESTdir3D] - corner[WSTdir3D] ) * eta +
       ( corner[WNTdir3D] - corner[WSTdir3D] ) * xi +
       ( corner[ENTdir3D] - corner[ESTdir3D] - corner[WNTdir3D] + corner[WSTdir3D] ) * eta * xi;

  Pres = Pres -
         2.0 * ( PD + ( PT-PD ) * phi );

  return Pres;
}









#endif 
