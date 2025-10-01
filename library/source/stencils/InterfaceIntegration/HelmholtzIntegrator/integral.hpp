///
///  @file  integral.hpp
///         Provides functions for sparse grid quadrature on rectilinear domains
///
///  @date  2016-06-12
///
///  @author  Tim Rheinfels
///

#ifndef INTEGRAL_HPP
#define INTEGRAL_HPP

#include <functional>
#include <iostream>
#include <limits>

#include "util.hpp"

///
///  @brief  Auxiliary function for evaluating boundary base functions
///
///  @tparam  T  Datatype used for calculation
///  @tparam  F  Type of the functor providing the integrand
///  @tparam  N  Number of dimensions
///
///  @param[in]   f            Function to evaluate
///  @param[in]   f_evaluated  Previously evaluted function values
///  @param[in]   x1           First point of the rectilinear domain
///  @param[in]   x2           Second point of the rectilinear domain
///  @param[in]   volume       Volume of the rectilinear domain
///  @param[in]   epsilon      Refinement boundary for the hierarchical surplus
///  @param[in]   level        Multi-index containing the base function's level
///  @param[in]   index        Multi-index containing the base function's index
///  @param[in]   dimension    Dimension that is being traversed in the current function call
///  @param[out]  refine       Flag for propagating refinement the higher dimensions
///  @param[in]   depth        Current depth of recursion
///  @param[in]   depth_max    Upper limit for the recursion depth (Neccessary for discontinious functions)
///  @param[in]   depth_min    Lower limit for the recursion depth (Might be useful for functions with small curvature in the domain's center)
///  @param[in]   gamma        Value by which an infinte value is replaced (For singularities inside of the domain)
///
///  @returns  Portion of the integral over f on the rectilinear domain, that is generate by the base function described by \p level and \p index
///
template <typename T, typename F, size_t N>
    T integral_aux_boundary( const F& f, const std::array<T, Util::pow(3, N)>& f_evaluated, std::array<T, N>& x1, std::array<T, N> x2, const T& volume, const T& epsilon, std::array<size_t, N>& level, std::array<size_t, N>& index, size_t dimension, bool& refine, size_t depth, size_t depth_min, size_t depth_max, const T& gamma )
{

    ///
    ///  @brief  Lambda for perfoming refinements (Additional levels)
    ///
    auto perform_refinement = [&] () -> T
        {

            // Do not refine when maximum depth is reached or base function is a boundary for the dimension
            if( (depth >= depth_max) || (level[dimension] == 0) )
            {
                return (T) 0.0;
            }

            // Calculate stride for traversing f_evaluated
            size_t stride = 1;
            for( size_t i = 0; i < dimension; ++i )
            {
                stride *= 3;
            }

            // Create data structures to evaluate and store function values for the next hierarchical level
            std::array<T, Util::pow(3, N)> f_evaluated_l;
            std::array<T, Util::pow(3, N)> f_evaluated_r;

            size_t refinement_index = 0;
            std::array<T, N> position;
    
            // Traverse the N-dimensional rectilinear domain distributing f_evaluated along the new arrays, evaluating f when necessary
            for( size_t i = 0; i < Util::pow<size_t>(3, N); ++i )
            {

                // Variable for determining combinations
                size_t tmp = i;

                // Iterate over the dimensions setting position to the correct value
                for( size_t k = 0; k < N; ++k )
                {

                    // Skip refinement dimension
                    if( k == dimension )
                    {
                        refinement_index = (tmp%3);
                        tmp /= 3;
                        continue;
                    }

                    switch( (tmp%3) )
                    {
                        // Left
                        case 0:
                            position[k] = x1[k];
                            break;
                        
                            // Center
                        case 1:
                            position[k] = (x1[k] + x2[k]) * (T) 0.5;
                            break;
                        
                            // Right
                        case 2:
                            position[k] = x2[k];
                            break;
                        
                            // This label is only for compiler compatibility, it cannot be reached
                        default:
                            break;
                    }
                
                    // Shift the second digit to the first place
                    tmp /= 3;
                
                }

                // Adjust position in refinement dimension
                switch( refinement_index )
                {

                    // Left
                    case 0:
                        f_evaluated_l[i] = f_evaluated[i];
                        f_evaluated_r[i] = f_evaluated[i + stride];
                        break;

                        // Center
                    case 1:
                        // Interpolate position linearily
                        position[dimension] = x1[dimension] + (x2[dimension] - x1[dimension]) * (T) 0.25;
                        f_evaluated_l[i] = f(position);
                        // Interpolate position linearily
                        position[dimension] = x1[dimension] + (x2[dimension] - x1[dimension]) * (T) 0.75;
                        f_evaluated_r[i] = f(position);
                        break;

                        // Right
                    case 2:
                        f_evaluated_l[i] = f_evaluated[i - stride];
                        f_evaluated_r[i] = f_evaluated[i];
                        break;

                    default:
                        break;
                
                }

            }
        
            // Correction achieved by the higher level basis functions
            T correction = (T) 0.0;

            // Temporary refinement flag
            bool refine2 = false;

            // Calculate the center of the corresponding dimension
            T center = (x1[dimension] + x2[dimension]) * (T) 0.5;

            // Increase level vector in the current component
            ++level[dimension];

            // Calculate the index of the left base function
            index[dimension] *= 2;

            // Trim the domain to the left half
            T tmp = x2[dimension];
            x2[dimension] = center;

            // Perform correction using the left base function
            correction += integral_aux_boundary( f, f_evaluated_l, x1, x2, volume*(T)0.5, epsilon, level, index, dimension, refine2, depth+1, depth_min, depth_max, gamma );

            // Calculate the index of the right base function
            ++index[dimension];

            // Trim the domain to right half
            x2[dimension] = tmp;
            tmp = x1[dimension];
            x1[dimension] = center;

            // Perform correction using the right base function
            correction += integral_aux_boundary( f, f_evaluated_r, x1, x2, volume*(T)0.5, epsilon, level, index, dimension, refine2, depth+1, depth_min, depth_max, gamma );

            // Reset state
            x1[dimension] = tmp;

            // Restore old state for level and index
            --index[dimension];
            index[dimension] /= 2;

            --level[dimension];

            // Return correction
            return correction;

        };

    // Perform the actual integration
    if( dimension > 0 )
    {

        // Dimensions above zero will only call the next lower dimension and refine when needed
        T local_integral = integral_aux_boundary( f, f_evaluated, x1, x2, volume, epsilon, level, index, dimension-1, refine, depth, depth_min, depth_max, gamma );

        // If refinement was requests on dimension 0, we refine them all (octree refinement)...
        if( refine )
        {
            // ...by incrementing the local integral by the correction from the higher level base functions
            local_integral += perform_refinement();
        }

        // return the local integral
        return local_integral;

    }
    else
    {
        
        // Dimension 0 is where the actual work takes place
        
        // Accumulator for the hierarchical surplus
        T w = (T) 0.0;
        
        // Traverse the N-dimensional rectilinear domain at all combinations of corner and center points
        for( size_t i = 0; i < Util::pow<size_t>(3, N); ++i )
        {
	
            // Variable for determining combinations
            size_t tmp = i;
            
            // Variable storing the weight of the function evaluation
            T wcoeff = (T) 1.0;
            
            // Flag to indicate if the combination stored in i should be used (Boundary dimensions are evaluated at one point, instead of three)
            bool use = true;
            
            // Iterate over the dimensions setting position to the right value and calculating wcoeff
            for( size_t k = 0; k < N; ++k )
            {
                
                // If k is a boundary dimension and the corresponding index in i is not the index of the boundary function,
                // we skip the combination
                if( (level[k] == 0) && ( (tmp%3) != index[k] ) )
                {
                    use = false;
                    break;
                }
                
                if( level[k] > 0 )
                {
                    if( (tmp%3) != 1 )
                    {
                        wcoeff *= (T) -0.5;
                    }
                }
                
                // Shift the second digit to the first place
                tmp /= 3;
                
            }
            
            // Skip combination if dimension k is a boundary and position would not be centered in k
            if( !use )
            {
                continue;
            }
            
            // Increment the hierarchical surplus by the weighted function evaluation
            w += wcoeff * f_evaluated[i];
        }
        
        //  The base function's contribution is the product of two factors:
        //    2^(-N) * volume: The integral over the unscaled base function
        //    w:               The hierarchical surplus
        //
        T local_integral = (Util::pow<T>(0.5, N) * volume) * w;
        
        // Refine when the hierarchical surplus exceeds epsilon
        if( (depth < depth_min) || (Util::abs(w) > epsilon)  )
        {
            local_integral += perform_refinement();
            refine = true;
        }
        else
        {
            refine = false;
        }
        
        // Return local integral
        return local_integral;
        
    }

    // Just for compiler compatibility
    return (T) 0.0;

}

///
///  @brief  Auxiliary function for evaluating inner base functions
///
///  @tparam  T  Datatype used for calculation
///  @tparam  F  Type of the functor providing the integrand
///  @tparam  N  Number of dimensions
///
///  @param[in]   f            Function to evaluate
///  @param[in]   f_evaluated  Previously evaluted function values
///  @param[in]   x1           First point of the rectilinear domain
///  @param[in]   x2           Second point of the rectilinear domain
///  @param[in]   volume       Volume of the rectilinear domain
///  @param[in]   epsilon      Refinement boundary for the hierarchical surplus
///  @param[in]   dimension    Dimension that is being traversed in the current function call
///  @param[out]  refine       Flag for propagating refinement the higher dimensions
///  @param[in]   depth        Current depth of recursion
///  @param[in]   depth_max    Upper limit for the recursion depth (Neccessary for discontinious functions)
///  @param[in]   depth_min    Lower limit for the recursion depth (Might be useful for functions with small curvature in the domain's center)
///  @param[in]   gamma        Value by which an infinte value is replaced (For singularities inside of the domain)
///
///  @returns  Portion of the integral over f on the rectilinear domain, that is generate by the base function
///
///  @note  As opposed to integral_aux_boundary, we don't need to know which level or index the currently processed base function has, since we make no difference
///
template <typename T, typename F, size_t N>
    T integral_aux( const F& f, const std::array<T, Util::pow(3, N)>& f_evaluated, std::array<T, N>& x1, std::array<T, N>& x2, const T& volume, const T& epsilon, size_t dimension, bool& refine, size_t depth, size_t depth_min, size_t depth_max, const T& gamma )
{

    ///
    ///  @brief  Lambda for perfoming refinements (Additional levels)
    ///
    auto perform_refinement = [&] () -> T
        {

            // Do not refine when maximum depth is reached or base function is a boundary for the dimension
            if( depth >= depth_max )
            {
                return (T) 0.0;
            }
    
            // Calculate stride for traversing f_evaluated
            size_t stride = 1;
            for( size_t i = 0; i < dimension; ++i )
            {
                stride *= 3;
            }

            // Create data structures to evaluate and store function values for the next hierarchical level
            std::array<T, Util::pow(3, N)> f_evaluated_l;
            std::array<T, Util::pow(3, N)> f_evaluated_r;

            size_t refinement_index = 0;
            std::array<T, N> position;
  
            // Traverse the N-dimensional rectilinear domain distributing f_evaluated along the new arrays, evaluating f when necessary
            for( size_t i = 0; i < Util::pow<size_t>(3, N); ++i )
            {

                // Variable for determining combinations
                size_t tmp = i;

                // Iterate over the dimensions setting position to the right value and calculating wcoeff
                for( size_t k = 0; k < N; ++k )
                {

                    if( k == dimension )
                    {
                        refinement_index = (tmp%3);
                        tmp /= 3;
                        continue;
                    }
                    
                    switch( (tmp%3) )
                    {
                        // Left
                        case 0:
                            position[k] = x1[k];
                            break;
                        
                            // Center
                        case 1:
                            position[k] = (x1[k] + x2[k]) * (T) 0.5;
                            break;
                        
                            // Right
                        case 2:
                            position[k] = x2[k];
                            break;
                        
                            // This label is only for compiler compatibility, it cannot be reached
                        default:
                            break;
                    }
                
                    // Shift the second digit to the first place
                    tmp /= 3;
                
                }

                // Adjust position in refinement dimension
                switch( refinement_index )
                {

                    // Left
                    case 0:
                        f_evaluated_l[i] = f_evaluated[i];
                        f_evaluated_r[i] = f_evaluated[i + stride];
                        break;

                        // Center
                    case 1:
                        // Interpolate position linearily
                        position[dimension] = x1[dimension] + (x2[dimension] - x1[dimension]) * (T) 0.25;
                        f_evaluated_l[i] = f(position);
                        // Interpolate position linearily
                        position[dimension] = x1[dimension] + (x2[dimension] - x1[dimension]) * (T) 0.75;
                        f_evaluated_r[i] = f(position);
                        break;

                        // Right
                    case 2:
                        f_evaluated_l[i] = f_evaluated[i - stride];
                        f_evaluated_r[i] = f_evaluated[i];
                        break;

                    default:
                        break;
                
                }

            }
        
            // Correction achieved by the higher level basis functions
            T correction = (T) 0.0;

            // Temporary refinement flag
            bool refine2 = false;

            // Calculate the center of the corresponding dimension
            T center = (x1[dimension] + x2[dimension]) * (T) 0.5;

            // Trim the domain to the left half
            T tmp = x2[dimension];
            x2[dimension] = center;

            // Perform correction using the left base function
            correction += integral_aux( f, f_evaluated_l, x1, x2, volume*(T)0.5, epsilon, dimension, refine2, depth+1, depth_min, depth_max, gamma );

            // Trim the domain to right half
            x2[dimension] = tmp;
            tmp = x1[dimension];
            x1[dimension] = center;

            // Perform correction using the right base function
            correction += integral_aux( f, f_evaluated_r, x1, x2, volume*(T)0.5, epsilon, dimension, refine2, depth+1, depth_min, depth_max, gamma );
            x1[dimension] = tmp;

            // Return correction
            return correction;

        };

    // Perform the actual integration
    if( dimension > 0 )
    {

        // Dimensions above zero will only call the next lower dimension and refine when needed
        T local_integral = integral_aux( f, f_evaluated, x1, x2, volume, epsilon, dimension-1, refine, depth, depth_min, depth_max, gamma );

        // If refinement was requests on dimension 0, we refine them all (octree refinement)...
        if( refine )
        {
            // ...by incrementing the local integral by the correction from the higher level base functions
            local_integral += perform_refinement();
        }
        
        // return the local integral
        return local_integral;

    }
    else
    {
        
        // Dimension 0 is where the actual work takes place

        // Accumulator for the hierarchical surplus
        T w = (T) 0.0;
        
        // Traverse the N-dimensional rectilinear domain at all combinations of corner and center points
        for( size_t i = 0; i < Util::pow<size_t>(3, N); ++i )
        {
            
            // Variable for determining combinations
            size_t tmp = i;
	
            // Variable storing the weight of the function evaluation
            T wcoeff = (T) 1.0;
        
            // Iterate over the dimensions setting position to the right value and calculating wcoeff
            for( size_t k = 0; k < N; ++k )
            {

                if( (tmp%3) != 1 )
                {
                    wcoeff *= (T) -0.5;
                }
                
                // Shift the second digit to the first place
                tmp /= 3;
                
        }
            
            // Increment the hierarchical surplus by the weighted function evaluation
            w += wcoeff * f_evaluated[i];
            
        }
        
        //  The base function's contribution is the product of two factors:
        //    2^(-N) * volume: The integral over the unscaled base function
        //    w:               The hierarchical surplus
        //
        T local_integral = (Util::pow<T>(0.5, N)) * volume * w;
        
        // Refine when the hierarchical surplus exceeds epsilon
        if( (depth < depth_min) || (Util::abs(w) > epsilon)  )
        {
            local_integral += perform_refinement();
            refine = true;
        }
        else
        {
            refine = false;
        }
        
        // Return local integral
        return local_integral;

    }
    
    // Just for compiler compatibility
    return (T) 0.0;

}

///
///  @brief  Uses sparse grid quadrature to approximate the integral over \p f on the rectilinear domain spanned by the points \p x1 and \p x2 using \p epsilon as the refinement criterion on the hierarchical surplus
///
///  @tparam  T  Datatype used for calculation
///  @tparam  F  Type of the functor providing the integrand
///  @tparam  N  Number of dimensions
///
///  @param[in]   f          Function to evaluate
///  @param[in]   x1         First point of the rectilinear domain
///  @param[in]   x2         Second point of the rectilinear domain
///  @param[in]   epsilon    Refinement boundary for the hierarchical surplus
///  @param[in]   depth_max  Upper limit for the recursion depth (Neccessary for discontinious functions)
///  @param[in]   depth_min  Lower limit for the recursion depth (Might be useful for functions with small curvature in the domain's center)
///  @param[in]   gamma      Value by which an infinte value is replaced (For singularities inside of the domain)
///
///  @returns  Approximation of Integral over f on the rectilinear domain spanned by the points \p x1 and \p x2
///
template <typename T, typename F, size_t N>
    T integral( const F& f, std::array<T, N>& x1, std::array<T, N> x2, const T& epsilon, size_t depth_min = 0, size_t depth_max = std::numeric_limits<size_t>::max(), const T& gamma = (T) 1.0e12 )
{

    // Accumulator for the integral
    T integral = (T) 0.0;

    // Calculate the domain's volume
    T volume = (T) 1.0;

    for( size_t i = 0; i < N; ++i )
    {
        volume *= (x2[i] - x1[i]);
    }

    // === Determine the portion of the outermost 2^N boundary points ===

    // Array for storing point to evaluate f in
    std::array<T, N> position;
    std::array<T, Util::pow(3, N)> f_evaluated;
    
    // Iterate over all 2^N combinations of edges
    for( size_t i = 0; i < Util::pow<size_t>(3, N); ++i )
    {

        bool boundary = true;

        size_t tmp = i;
        
        // Iterate over the dimensions
        for( size_t k = 0; k < N; ++k )
        {

            // Adjust position to represent the corner point encoded by i
            switch( tmp % 3 )
            {
                
                case 0:
                    position[k] = x1[k];
                    break;                

                case 1:
                    position[k] = (x1[k] + x2[k]) * (T) 0.5;
                    boundary = false;
                    break;
                    
                case 2:
                    position[k] = x2[k];
                    break;

                default:
                    break;

            }

            tmp /= 3;
                    
        }


        f_evaluated[i] = f(position);

        if( boundary )
        {
            integral += volume * (Util::pow<T>(0.5, N)) * f_evaluated[i];
        }
        
    }

    // === Adaptively calculate the portion of all other boundary points ===

    // Refinement flag
    bool refine = false;
    
    // Allocate multi indices storing the level and index of the corresponding base functions
    std::array<size_t, N> level;
    std::array<size_t, N> index;

    // Iterate over all 2^N combinations of corner points creating the levels of the base functions
    for( size_t i = 0; i < Util::pow<size_t>(2, N); ++i )
    {

        // Counter for the number of non-boundary dimensions
        size_t ones = 0;

        // Adjust level to represent the base function level encoded by i
        for( size_t k = 0; k < N; ++k )
        {
            if( (i & (1<<k)) > 0 )
            {
                level[k] = 1;
                ++ones;
            }
            else
            {
                level[k] = 0;
            }
        }

        // Skip the outermost boundary base functions (0, 0, ..., 0) and the inner base functions (1, 1, ..., 1)
        if( (ones == 0) || (ones == N) )
        {
            continue;
        }

        // Again iterate over all 2^N combinations of corner points, but this time creating the indices of the base functions
        for( size_t j = 0; j < Util::pow<size_t>(2, N); ++j )
        {

            // Flag indicating, whether to use the generated index or not
            bool use = true;

            // Again iterate over all dimensions
            for( size_t k = 0; k < N; ++k )
            {

                // index[k] encodes the position in the domain
                //   0: Left
                //   1: Right
                //
                index[k] = (( j & (1<<k) ) >> k) * 2;

                // Base functions (1) (2) do not exist, therefore skip them
                if( (level[k] == 1) && (index[k] == 2) )
                {
                    use = false;
                    break;
                }

            }

            // Skip non-existant base functions
            if( !use )
            {
                continue;
            }

            // Add the portion of the base function described by level and index to the integral accumulator
            integral += integral_aux_boundary( f, f_evaluated, x1, x2, volume, epsilon, level, index, N-1, refine, 0, depth_min, depth_max, gamma );

        }

    }

    // ===  Calculate the portion of the inner point base functions ===
    integral += integral_aux( f, f_evaluated, x1, x2, volume, epsilon, N-1, refine, 0, depth_min, depth_max, gamma );
    
    // === Return the approximated integral ===
    return integral;

}

#endif /* INTEGRAL_HPP */ 
