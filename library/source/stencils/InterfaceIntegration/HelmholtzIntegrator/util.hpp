///
///  @file  util.hpp
///         Provides utility functions
///
///  @date  2016-06-12
///
///  @author  Tim Rheinfels
///

#ifndef UTIL_HPP
#define UTIL_HPP

namespace Util
{

    ///
    ///  @brief  Calculates the \p e-th power of \p b
    ///
    ///  @tparam  Datatype of the base
    ///
    ///  @param[in]  b  Base
    ///  @param[in]  e  Exponent
    ///
    ///  @returns  b^e
    ///
    template <typename T>
	inline constexpr T pow( const T& b, size_t e )
    {
	return ( (e==0) ? 1.0 : b * pow(b, e-1) );
    }
    
///
///  @brief  Calculates the absolute value of \p x
///
///  @tparam  T  Datatype of \p x
///
///  @param[in]  x  Operand
///
///  @returns  |x|
///
template <typename T>
inline T abs( const T& x )
{
    return ( (x>(T)0.0) ? x : (-x) );
}

}

#endif /* UTIL_HPP */
