// Function-like platform macros must not alter the stored arithmetic header.
// Preload its dependencies so this control isolates this header's call sites.
#include <gp_XYZ.hxx>
#include <Standard_ConstructionError.hxx>
#include <algorithm>
#include <array>
#include <cfenv>
#include <cfloat>
#if defined(__linux__) && defined(__x86_64__) && defined(__SSE2__)
  #include <xmmintrin.h>
#endif
#include <cmath>
#include <limits>

#ifdef _gp_StoredRepresentation_HeaderFile
  #error The target header must not already have been included by a dependency
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

// Select the checkout header explicitly, not an installed predecessor.
#include "../../src/gp/gp_StoredRepresentation.hxx"

#if !defined(min) || !defined(max)
  #error The header must preserve the caller's platform macros
#endif
static_assert(min(3, 7) == 3 && max(3, 7) == 7, "caller macros remain usable");

int main()
{
  const auto negative = gp_StoredRepresentation::Abs({-3.0, -2.0});
  const auto crossing = gp_StoredRepresentation::Abs({-3.0, 2.0});
  if (negative.Lo != 2.0 || negative.Hi != 3.0 || crossing.Lo != 0.0 || crossing.Hi != 3.0)
    return 1;
  gp_StoredRepresentation::Orthogonal(gp_XYZ(1, 0, 0), gp_XYZ(0, 1, 0));
  return min(0, 1) + max(-1, 0);
}
