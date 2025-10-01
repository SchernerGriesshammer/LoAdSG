/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer, Tim Rheinfels
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#ifndef INTERATORHELM_H
#define INTERATORHELM_H

#include <cmath>
#include <functional>
#include <limits>

#include "integral.hpp"

///
///  @brief  Integrator for the Helmholtz operator with variable coefficients
///
///  @tparam  F  Functor type for the variable coefficient
///  @tparam  N  Number of dimensions
///
///  @note  Object of type F need to be copyable and implement an operator() ( const std::array<double, N>& )
///
template <typename F, size_t N>
    class IntegratorHelm : public InterfaceLocalStiffnessMatrices<N>
{

    private:

        ///
        ///  @brief  Class for describing a Helmholtz term
        ///
        class HelmholtzTerm
        {
                
            private:
                
                std::array<double, N> m_p_left;        ///<  First point spanning the domain
                std::array<double, N> m_p_right;       ///<  Second point spanning the domain
                std::array<BasisFunctionType, N> m_u;  ///<  Base function tensor for the first d-dimensional base function
                std::array<BasisFunctionType, N> m_v;  ///<  Base function tensor for the second d-dimensional base function
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
                inline HelmholtzTerm( const std::array<double, N>& p_left, const std::array<double, N>&  p_right,
                                      const std::array<BasisFunctionType, N>& u, const std::array<BasisFunctionType, N>& v,
                                      const F& variable_coefficient, double gamma )
                    : m_p_left(p_left),
                      m_p_right(p_right),
                      m_u(u),
                      m_v(v),
                      m_variable_coefficient(variable_coefficient),
                      m_gamma(gamma)
                {
                    
                }

                ///
                ///  @brief  Evaluation operator
                ///
                ///  @param[in]  x  Position to evaluate the Helmholtz term at
                ///
                ///  @returns  The evaluation of the Helmholtz term at \p x
                ///
                inline double operator() ( const std::array<double, N>& x ) const
                {
                    
                    double value = 1.0;
                    
                    for( size_t i = 0; i < N; ++i )
                    {
                        
                        // === Muliply by the first base function ===
                        switch( m_u[i] )
                        {
                            case BasisFunctionType::leftBasis:
                                // p(x) = (x-b)/(a-b)
                                value *= (x[i] - m_p_right[i]) / (m_p_left[i] - m_p_right[i]);
                                break;
                                
                            case BasisFunctionType::rightBasis:
                                // p(x) = (x-a)/(b-a)
                                value *= (x[i] - m_p_left[i]) / (m_p_right[i] - m_p_left[i]);
                                break;
                                
                            case BasisFunctionType::gradLeftBasis:
                                // p'(x) =  1/(b-a)
                                value /= (m_p_left[i] - m_p_right[i]);
                                break;
                                
                            case BasisFunctionType::gradRightBasis:
                                // p'(x) = 1/(a-b)
                                value /= (m_p_right[i] - m_p_left[i]);
                                break;
                        }
                        
                        // === Muliply by the second base function ===
                        switch( m_v[i] )
                        {
                            case BasisFunctionType::leftBasis:
                                // p(x) = (x-b)/(a-b)
                                value *= (x[i] - m_p_right[i]) / (m_p_left[i] - m_p_right[i]);
                                break;
                                
                            case BasisFunctionType::rightBasis:
                                // p(x) = (x-a)/(b-a)
                                value *= (x[i] - m_p_left[i]) / (m_p_right[i] - m_p_left[i]);
                                break;
                                
                            case BasisFunctionType::gradLeftBasis:
                                // p'(x) =  1/(b-a)
                                value /= (m_p_left[i] - m_p_right[i]);
                                break;
                                
                            case BasisFunctionType::gradRightBasis:
                                // p'(x) = 1/(a-b)
                                value /= (m_p_right[i] - m_p_left[i]);
                                break;
                        }

                    }
                    
                    // Evaluate the variable coefficient replacing floating point inf by gamma
                    double c_eval = m_variable_coefficient(x);
                    
                    if( (std::abs(c_eval) > m_gamma) || (std::fpclassify(c_eval) == FP_INFINITE ))
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
        IntegratorHelm( const F& c, double epsilon, size_t depth_min = 0, size_t depth_max = std::numeric_limits<size_t>::max(), double gamma = std::numeric_limits<double>::max() )
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
        double stencil_integration(double p_left[N], double p_right[N], BasisFunctionType u[N], BasisFunctionType v[N]) const
        {

            // Duplicate input data as std::array
            std::array<double, N> p_left_array;
            std::array<double, N> p_right_array;

            std::array<BasisFunctionType, N> u_array;
            std::array<BasisFunctionType, N> v_array;
            
            for( size_t i = 0; i < N; ++i )
            {

                p_left_array[i] = p_left[i];
                p_right_array[i] = p_right[i];

                u_array[i] = u[i];
                v_array[i] = v[i];
                
            }
            
            // Create local helmholtz term
            HelmholtzTerm ht(p_left_array, p_right_array, u_array, v_array, m_variable_coefficient, m_gamma);

            
            // Integrate
            return integral<double, HelmholtzTerm, N>( ht, p_left_array, p_right_array, m_epsilon, m_depth_min, m_depth_max, m_gamma );

        }

};

#endif
