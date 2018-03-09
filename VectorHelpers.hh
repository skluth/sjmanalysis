#ifndef VECTORHELPERS_HH
#define VECTORHELPERS_HH

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <math.h>

template <typename T>
std::vector<T> multiplyVectors( std::vector<T> lhs, const std::vector<T>& rhs ) {
  return lhs*rhs;
}

// Default behaviour appropriate for calculation of bin-by-bin correction factors
template <typename T>
class DivisorChecked {
  bool lthrow;
public:
  DivisorChecked( bool lthr=false ) : lthrow(lthr) {}
  T operator()( const T lhs, const T rhs ) const {
    if( rhs == 0.0 ) {
      if( lthrow ) throw std::runtime_error( "divide by zero" );
      else return 0.0;
    }
    return lhs/rhs;
  }
};
template <typename T>
std::vector<T> divideChecked( std::vector<T> lhs, const std::vector<T>& rhs,
			      bool lthrow=false ) {
  DivisorChecked<T> divisor( lthrow );
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
		  divisor );
  return lhs;
}

template <typename T>
std::vector<T> subtractVectors( const std::vector<T>& lhs,
				const std::vector<T>& rhs ) {
  return lhs-rhs;
}


template <typename T>
std::vector<T> operator+( std::vector<T> lhs,
			  const std::vector<T>& rhs ) {
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
		  std::plus<T>() );
  return lhs;
}
template <typename T>
std::vector<T> operator-( std::vector<T> lhs,
			  const std::vector<T>& rhs ) {
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
		  std::minus<T>() );
  return lhs;
}
template <typename T>
std::vector<T> operator*( std::vector<T> lhs,
			  const std::vector<T>& rhs ) {
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
		  std::multiplies<T>() );
  return lhs;
}
template <typename T>
std::vector<T> operator*( std::vector<T> lhs,
			  T rhs ) {
  std::transform( lhs.begin(), lhs.end(), lhs.begin(),
		  [rhs]( T lhs ) { return lhs*rhs; } );
  return lhs;
}
template <typename T>
std::vector<T> operator*(  T lhs, std::vector<T> rhs ) {
  std::transform( rhs.begin(), rhs.end(), rhs.begin(),
		  [lhs]( T rhs ) { return lhs*rhs; } );
  return rhs;
}
template <typename T>
std::vector<T> operator/( std::vector<T> lhs,
			  const std::vector<T>& rhs ) {
  std::transform( lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
		  std::divides<T>() );
  return lhs;
}
template <typename T>
std::vector<T> sqrt( std::vector<T> in ) {
  std::transform( in.begin(), in.end(), in.begin(),
		  []( T in ) { return sqrt( in ); } );
  return in;
}
template <typename T>
std::vector<T> square( std::vector<T> in ) {
  std::transform( in.begin(), in.end(), in.begin(),
		  []( T in ) { return in*in; } );
  return in;
}

#endif
