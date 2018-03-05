#ifndef NORMALISATION_HH
#define NORMALISATION_HH

#include "TObject.h"

// Helpers to set and get TObject bit to store normalisation status
inline void SetIsNormalised( TObject & obj, bool value ) {
  obj.SetBit( 1<<14, value );
}
inline bool IsNormalised( const TObject & obj ) {
  return obj.TestBit( 1<<14 );
}

#endif
