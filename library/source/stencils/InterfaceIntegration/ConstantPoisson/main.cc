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

//test NOW
//using namespace std;

#define Dim 2

int main(int argc, char** argv) {
  cout.precision(5);
  cout.setf(std::ios::fixed,std::ios::floatfield);

  
  double p_left[Dim];
  double p_right[Dim];
    
  BasisFunctionType u[Dim];
  BasisFunctionType v[Dim];
  
  // Example: Constant Poisson
  //////////////////////////////////
  
  cout << " Poisson!!!!!!!!!!! " << endl;
  cout << " ---------------------" << endl;
  
  IntegratorPoisson<Dim> localStiffnessPoisson;

  for(unsigned i=0;i<Dim;++i) {
      p_left[i]  = 0.0;
      p_right[i] = 0.125;      
  }
  

  for(IteratorBasisFunction<Dim> iterU;iterU.hasNext();iterU.next()) {
      cout << endl;
      cout << " u = (";    
      for(int d=0;d<Dim;++d) {
	  cout << iterU.getBasisTypeNum(d);
	  if(d<Dim-1) cout << ", ";
      }
      cout << ")" << endl;
      cout << " ----------- " << endl;
      for(IteratorBasisFunction<Dim> iterV;iterV.hasNext();iterV.next()) {
          cout << " v = (";    
          for(int d=0;d<Dim;++d) {
	      cout << iterV.getBasisTypeNum(d);
	      if(d<Dim-1) cout << ", ";
          }
          
	  for(int d=0;d<Dim;++d) {
	    u[d] =iterU.getBasisTypeCoord(d);
	    v[d] =iterV.getBasisTypeCoord(d);	    
	  }
          
          cout << "), Integral =  "
	       << localStiffnessPoisson.stencil_integration(p_left,p_right,u,v)
	       << endl;
      }
  }
    
    
    
    
    
    
 // Example: Constant Helm
  //////////////////////////////////

  cout << endl << endl;
  cout << " Helm!!!!!!!!!!! " << endl;
  cout << " ---------------------" << endl;
  
  
  IntegratorHelmConstant<Dim> localStiffnessHelm(100.0);

  for(unsigned i=0;i<Dim;++i) {
      p_left[i]  = 0.0;
      p_right[i] = 0.125;      
  }
  
  for(IteratorBasisFunction<Dim> iterU;iterU.hasNext();iterU.next()) {
      cout << endl;
      cout << " u = (";    
      for(int d=0;d<Dim;++d) {
	  cout << iterU.getBasisTypeNum(d);
	  if(d<Dim-1) cout << ", ";
      }
      cout << ")" << endl;
      cout << " ----------- " << endl;
      for(IteratorBasisFunction<Dim> iterV;iterV.hasNext();iterV.next()) {
          cout << " v = (";    
          for(int d=0;d<Dim;++d) {
	      cout << iterV.getBasisTypeNum(d);
	      if(d<Dim-1) cout << ", ";
          }
          
	  for(unsigned d=0;d<Dim;++d) {
	    u[d] =iterU.getBasisTypeCoord(d);
	    v[d] =iterV.getBasisTypeCoord(d);	    
	  }
          
          cout << "), Integral =  "
	       << localStiffnessHelm.stencil_integration(p_left,p_right,u,v)
	       << endl;
      }
  }
    
    
    
    
    
    

    
    
    
    
    
    
  cout << " Ende! " << endl;
}


