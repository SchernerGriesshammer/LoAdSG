/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer, Rainer Hartmann
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/




/**
 * integrator of Poisson equation with constant coefficient
 */
template <size_t DIM>
class IntegratorPoisson : public InterfaceLocalStiffnessMatrices<DIM> {
   public:
      IntegratorPoisson() {}
      double stencil_integration(double p_left[DIM],   double p_right[DIM],
				                         BasisFunctionType u[DIM], BasisFunctionType v[DIM]) const;
};


template <size_t DIM>
double IntegratorPoisson<DIM>::stencil_integration(double p_left[DIM],   double p_right[DIM],
				                   BasisFunctionType u[DIM], BasisFunctionType v[DIM]) const {
  static double h[DIM];
  for(int i=0;i<DIM;++i) h[i] = p_right[i] - p_left[i];
  
  double sum=0.0;
  for(int i=0;i<DIM;++i) {  // int du/dxi * dv/dxi  d(x1...xd)
      double prod=1.0;
      for(int j=0;j<DIM;++j) { // integral factor in direction j 
          if(i==j) {           // factor: int du/dxi * dv/dxi dxi
	     if(u[j]==v[j]) prod =   prod / h[j];
	     else           prod = - prod / h[j];
	  }
	  else {               // factor: int du/dxi * dv/dxi dxj  == int u*v dxj
	     if(u[j]==v[j]) prod =  prod * (1.0/3.0) * h[j];
	     else           prod =  prod * (1.0/6.0) * h[j];	    
	  }
      }
      sum = sum + prod;
  }
  return sum;   // sum = sum_i int du/dxi * dv/dxi  d(x1...xd)

}



/**
 * integrator of Helmholtz with constant coefficient
 */
template <size_t DIM>
class IntegratorHelmConstant : public InterfaceLocalStiffnessMatrices<DIM> {
   public:
      IntegratorHelmConstant(double valueC_) : valueC(valueC_) {}
      double stencil_integration(double p_left[DIM],   double p_right[DIM],
				                         BasisFunctionType u[DIM], BasisFunctionType v[DIM]) const;
   private:
     double valueC;
};


template <size_t DIM>
double IntegratorHelmConstant<DIM>::stencil_integration(double p_left[DIM],   double p_right[DIM],
				                   BasisFunctionType u[DIM], BasisFunctionType v[DIM]) const {
  static double h[DIM];
  for(int i=0;i<DIM;++i) h[i] = p_right[i] - p_left[i];
  
  double prod=1.0;
  for(int j=0;j<DIM;++j) { // integral factor in direction j 
    if(u[j]==v[j]) prod =  prod * (1.0/3.0) * h[j];
    else           prod =  prod * (1.0/6.0) * h[j];	    
  }
  return prod * valueC;
}



