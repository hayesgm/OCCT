// Bounded binary64 stored-representation applicability checks.
#ifndef _gp_StoredRepresentation_HeaderFile
#define _gp_StoredRepresentation_HeaderFile

#include <gp_XYZ.hxx>
#include <Standard_ConstructionError.hxx>
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>

//! Arithmetic-domain checks for exact restoration, never an equivalence tolerance.
//! Inputs and nonzero arithmetic intermediates must be normal finite binary64 values.
//! Round-to-nearest is required. Overflow/underflow and unsupported domains refuse.
namespace gp_StoredRepresentation
{
struct Interval { double Lo, Hi; };
inline void Check(bool theOK, const char* theMessage)
{
  if (!theOK) Standard_ConstructionError::Raise(theMessage);
}
inline void Arithmetic()
{
  Check(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53
    && std::fegetround() == FE_TONEAREST, "stored geometry requires nearest binary64 arithmetic");
}
inline void Finite(double theValue)
{
  Check(std::isfinite(theValue), "stored geometry has a nonfinite coefficient");
}
inline void NormalOrZero(double theValue)
{
  Finite(theValue);
  Check(theValue == 0.0 || std::isnormal(theValue), "stored geometry arithmetic subnormal domain");
}
inline double AddEnd(double a, double b, bool theUp)
{
  const double r = a + b;
  NormalOrZero(r);
  if (r == 0.0) {
    Check(a == -b, "stored geometry addition underflow"); return r;
  }
  if (a == 0.0 || b == 0.0) return r;
  return std::nextafter(r, theUp ? INFINITY : -INFINITY);
}
inline double MulEnd(double a, double b, bool theUp)
{
  const double r = a * b;
  NormalOrZero(r);
  if (a == 0.0 || b == 0.0) return r;
  Check(r != 0.0, "stored geometry multiplication underflow");
  return std::nextafter(r, theUp ? INFINITY : -INFINITY);
}
inline Interval Add(Interval a, Interval b)
{
  return {AddEnd(a.Lo,b.Lo,false), AddEnd(a.Hi,b.Hi,true)};
}
inline Interval Neg(Interval a) { return {-a.Hi,-a.Lo}; }
inline Interval Mul(Interval a, Interval b)
{
  const std::array<double,4> lo = {MulEnd(a.Lo,b.Lo,false),MulEnd(a.Lo,b.Hi,false),
    MulEnd(a.Hi,b.Lo,false),MulEnd(a.Hi,b.Hi,false)};
  const std::array<double,4> hi = {MulEnd(a.Lo,b.Lo,true),MulEnd(a.Lo,b.Hi,true),
    MulEnd(a.Hi,b.Lo,true),MulEnd(a.Hi,b.Hi,true)};
  return {*std::min_element(lo.begin(),lo.end()),*std::max_element(hi.begin(),hi.end())};
}
inline Interval Point(double a) { NormalOrZero(a); return {a,a}; }
inline Interval Abs(Interval a)
{
  return {a.Lo <= 0.0 && a.Hi >= 0.0 ? 0.0 : std::min(std::abs(a.Lo),std::abs(a.Hi)),
    std::max(std::abs(a.Lo),std::abs(a.Hi))};
}
inline Interval Dot(const gp_XYZ& a, const gp_XYZ& b)
{
  Interval r = {0,0};
  for (int i=1;i<=3;++i) r=Add(r,Mul(Point(a.Coord(i)),Point(b.Coord(i))));
  return r;
}
inline Interval Gamma(unsigned n)
{
  const double nu = n * 0x1p-53; // Exact for these small integer operation counts.
  const double denominator = 1.0 - nu; // Exact binary64 here.
  const double r = nu / denominator;
  return {std::nextafter(r,-INFINITY),std::nextafter(r,INFINITY)};
}
inline Interval FrameAngularEnvelope()
{
  // This is an explicit algebraic stored-representation invariant, not
  // provenance about how any saved frame was produced. The largest admitted
  // single-axis primitive has twelve dependent binary64 roundings:
  // matrix-vector (3 products + 2 sums), squared norm (3 products + 2 sums),
  // sqrt, and component division. Cross-plus-normalize has ten and plain
  // normalize has seven. Two independently represented axes therefore use at
  // most 24 and their three-term dot uses 5. Gamma(32) is the fixed next
  // enclosing budget, with three additional dependent-rounding slots as a
  // deliberate compositional guard. The interval evaluation independently
  // encloses the exact stored dot and does not consume that error budget.
  // Inputs outside it refuse even if some unrecorded producer history might
  // otherwise explain them.
  return Gamma(32);
}
inline void Unit(const gp_XYZ& a)
{
  Arithmetic();
  // Existing normalization: positive square/sum (three rounding factors), sqrt
  // (squared factor), and division (squared factor) imply norm^2 within
  // [(1-u)^2/(1+u)^5, (1+u)^2/(1-u)^5], enclosed by 1 +/- gamma_7.
  const Interval norm = Dot(a,a);
  const Interval gamma = Gamma(7);
  const double lo = AddEnd(1.0,-gamma.Hi,false);
  const double hi = AddEnd(1.0,gamma.Hi,true);
  Check(norm.Lo >= lo && norm.Hi <= hi, "stored direction outside binary64 unit domain");
}
inline void Orthogonal(const gp_XYZ& a, const gp_XYZ& b)
{
  const Interval dot = Dot(a,b);
  // The finite absolute angular envelope is dimensionless. It never compares,
  // changes, or snaps geometric lengths and carries no producer-path claim.
  const Interval angular = FrameAngularEnvelope();
  Check(std::max(std::abs(dot.Lo),std::abs(dot.Hi)) <= angular.Hi,
    "stored frame outside binary64 orthogonality domain");
}
inline Interval FrameDeterminant(const gp_XYZ& theZ, const gp_XYZ& theX, const gp_XYZ& theY)
{
  std::array<Interval,3> cross;
  for (int i=0;i<3;++i) {
    const int j=(i+1)%3+1,k=(i+2)%3+1;
    cross[i]=Add(Mul(Point(theX.Coord(j)),Point(theY.Coord(k))),
      Neg(Mul(Point(theX.Coord(k)),Point(theY.Coord(j)))));
  }
  Interval determinant={0,0};
  for (int i=0;i<3;++i) determinant=Add(determinant,Mul(cross[i],Point(theZ.Coord(i+1))));
  return determinant;
}
inline void Frame(const gp_XYZ& theZ, const gp_XYZ& theX, const gp_XYZ& theY)
{
  Unit(theZ); Unit(theX); Unit(theY);
  Orthogonal(theZ,theX); Orthogonal(theZ,theY); Orthogonal(theX,theY);
  const Interval determinant=FrameDeterminant(theZ,theX,theY);
  Check(determinant.Lo > 0.0 || determinant.Hi < 0.0,
    "stored frame has unresolved or zero handedness");
}
inline void RightHandedFrame(const gp_XYZ& theZ, const gp_XYZ& theX, const gp_XYZ& theY)
{
  Frame(theZ,theX,theY);
  Check(FrameDeterminant(theZ,theX,theY).Lo > 0.0,
    "stored right-handed frame has negative or unresolved determinant");
}
}
#endif
