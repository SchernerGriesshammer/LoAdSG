/**********************************************************************************
 * Copyright 2016 Christoph Pflaum, Tim Rheinfels
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

#ifndef INTERATORRHS_H
#define INTERATORRHS_H

#include <functional>
#include <limits>

#include "integral.hpp"
#include "../interfaceMatrices.h"

///
///  @brief  Integrator for the Helmholtz operator with variable coefficients
///
///  @tparam  F  Functor type for the variable coefficient
///  @tparam  N  Number of dimensions
///
///  @note  Object of type F need to be copyable and implement an operator() ( const std::array<double, N>& )
///
template <typename F, size_t N>
    class IntegratorRHS : public InterfaceLocalStiffnessMatrices<N>
{

    private:

        ///
        ///  @brief  Class for describing a Helmholtz term
        ///
        class TermRHS
        {
                
            private:
                
                std::array<double,N> h_u;              ///<  Meshwidth of Base function u
                std::array<double,N> h_v;              ///<  Meshwidth of Base function v
                std::array<double,N> p_u;              ///<  Position of Base function u;
                std::array<double,N> p_v;              ///<  Position of Base function v;
                std::array<int,N> d_u;              ///<  Depth of Base function u;
                std::array<int,N> d_v;              ///<  Depth of Base function v;
                F m_variable_coefficient;              ///<  Functor providing the variable coefficient
                double m_gamma;                        ///<  Cutoff parameter for singularities

                
            public:

                ///
                ///  @brief  Constructor
                ///
                ///  @param[in]  p_left                First point spanning the domain
                ///  @param[in]  p_right               Second point spanning the domain
                ///  @param[in]  u                     Base function tensor for the first d-dimensional base function
                ///  @param[in]  v                     Base function tensor for the second d-dimensional base function
                ///  @param[in]  variable_coefficient  Functor providing the variable coefficient
                ///  @param[in]  gamma                 Cutoff parameter for singularities
                ///
                inline TermRHS( const std::array<double, N>& h_u_, const std::array<double, N>&  h_v_, const std::array<double, N>& p_u_, const std::array<double, N>&  p_v_,const std::array<int, N>& d_u_, const std::array<int, N>&  d_v_, const F& variable_coefficient, double gamma )
                    : h_u(h_u_),
                      h_v(h_v_),
                      p_u(p_u_),
                      p_v(p_v_),
                      d_u(d_u_),
                      d_v(d_v_),
                      m_variable_coefficient(variable_coefficient),
                      m_gamma(gamma)
                {
                    
                }

                ///
                ///  @brief  Evaluation operator
                ///
                ///  @param[in]  x  Position to evaluate the RHS term at
                ///
                ///  @returns  The evaluation of the RHS term at \p x
                ///
                inline double operator() (std::array<double, N>& x ) const
                {


                    
                    double value = 1.0;
                    
                    for( size_t i = 0; i < N; i++ )
                    {
                        // Multiply by the first base function

                        if(d_u[i]!=0){
                            double suppleft=p_u[i]-h_u[i];
                            double suppright=p_u[i]+h_u[i];

                                if(suppleft <= x[i] && x[i] <= suppright){
                                    double slope = 1.0/h_u[i];
                                    if(x[i]<p_u[i]){
                                        double c = -(p_u[i]-h_u[i])/h_u[i];
                                        value *= (slope*x[i]+c);
                                    } else{
                                        double c = + (p_u[i]+h_u[i])/h_u[i];
                                        value *= (-1.0*slope*x[i]+c);
                                    }
                                }else{
                                    value = 0.0;
                                }
                        }else{
                            if(p_u[i]==0.0)
                                value *= -x[i]+1.0;
                            if(p_u[i]==1.0)
                                value *= x[i];
                        }
                        if(d_v[i]!=0) {
                            double suppleft = p_v[i] - h_v[i];
                            double suppright = p_v[i] + h_v[i];
                            if (suppleft <= x[i] && x[i] <= suppright) {
                                double slope = 1.0 / h_v[i];
                                if (x[i] < p_v[i]) {

                                    double c = -(p_v[i] - h_v[i]) / h_v[i];
                                    value *= (slope * x[i] + c);
                                } else {

                                    double c = +(p_v[i] + h_v[i]) / h_v[i];
                                    value *= (-1.0 * slope * x[i] + c);
                                }
                            } else {
                                value = 0.0;
                            }
                        }else{

                            if(p_v[i]==0.0){
                                value *= -x[i]+1.0;
                            }

                            if(p_v[i]==1.0){
                                value *= x[i];
                            }

                        }

                    }


                    
                    // Evaluate the variable coefficient replacing floating point inf by gamma
                    double c_eval = m_variable_coefficient(x);
                    
                    if( (std::abs(c_eval) > m_gamma) || (std::fpclassify(c_eval) == FP_INFINITE) )
                    {
                        c_eval = (std::signbit(c_eval) ? (-m_gamma) : m_gamma);
                    }
                    
                    // Mulitply by the variable coefficient
                    value *= c_eval;
                    
                    // Return Helmholtz function with variable coefficient evaluated at x
                    return value;

                }
                    
        };
        
        
        F m_variable_coefficient;  ///<  The variable coefficient
        double m_epsilon;          ///<  Refinement criterion on the hierarchical surplus
        size_t m_depth_min;        ///<  Minimum recursion depth
        size_t m_depth_max;        ///<  Maximum recursion depth
        double m_gamma;            ///<  Value to replace infinte values with

        
    public:
        
        ///
        ///  @brief  Constructor
        ///
        ///  @param[in]  c          The variable coefficient scalar field
        ///  @param[in]  epsilon    Refinement criterion of the hierarchical surplus
        ///  @param[in]  depth_max  Upper limit for the recursion depth (Neccessary for discontinious functions)
        ///  @param[in]  depth_min  Lower limit for the recursion depth (Might be useful for functions with small curvature in the domain's center)
        ///  @param[in]  gamma      Value by which an infinte value is replaced (For singularities inside of the domain)
        ///
        IntegratorRHS( const F& c, double epsilon, size_t depth_min = 0, size_t depth_max = std::numeric_limits<size_t>::max(), double gamma = std::numeric_limits<double>::max() )
            : m_variable_coefficient(c),
              m_epsilon(epsilon),
              m_depth_min(depth_min),
              m_depth_max(depth_max),
              m_gamma(gamma)
        {

            /* Nothing to do here... */

        }

        ///
        ///  @brief  Performs the integration of the FEM function defined by \p p_left, \p p_right, \p u and \p u with the additional variable coefficient passed to the constructor
        ///
        ///  @param[in]  p_left   First point of the rectilinear boundary
        ///  @param[in]  p_right  Second point of the rectilinear boundary
        ///  @param[in]  u        First set of base functions used in the bilinear form
        ///  @param[in]  v        Second set of base functions used in the bilinear form
        ///
        ///  @returns  The integral over the tensor product of (\p u * \p v) multiplied by the variable coefficient passed to the constructor on the domain spanned by \p p_left and \p p_right
        ///
        double integration(std::array<double, N> l_b, std::array<double, N> r_b, std::array<double, N> p_u, std::array<double, N> p_v, std::array<double, N> h_u, std::array<double, N> h_v,
                           array<int, N> d_u, std::array<int, N> d_v) const
        {

            TermRHS ht(h_u, h_v, p_u,p_v,d_u,d_v,m_variable_coefficient, m_gamma );





            // Integrate
            return integral<double, TermRHS, N>( ht, l_b, r_b, m_epsilon, m_depth_min, m_depth_max, m_gamma );

        }

        double stencil_integration(double p_left[N], double p_right[N],
                                   BasisFunctionType u[N], BasisFunctionType v[N]) const {return 0;};



};

#endif
