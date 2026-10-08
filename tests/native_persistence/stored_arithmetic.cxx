// Exact dyadic controls for stored-frame arithmetic; no coefficient reconstruction.
#include <gp_Ax3.hxx>
#include <gp_StoredRepresentation.hxx>
#include <Standard_Failure.hxx>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>

static void require(bool ok, const char* message)
{
  if (!ok) throw std::runtime_error(message);
}
static double number(std::uint64_t bits)
{
  double result;
  std::memcpy(&result, &bits, sizeof result);
  return result;
}
static std::uint64_t bits(double value)
{
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}
template<class Function> static void refused(Function function, const char* reason)
{
  try { function(); }
  catch (const Standard_Failure& failure)
  {
    const char* message = failure.GetMessageString();
    require(message != nullptr && std::strstr(message, reason) != nullptr, "wrong refusal reason");
    return;
  }
  throw std::runtime_error("missing refusal");
}
struct DyadicCase
{
  const char* Name;
  bool Multiply;
  std::uint64_t A, B, Rounded, Floor, Ceil;
  bool Refuse;
};
#include "stored_arithmetic_cases.hxx"
#include "stored_underflow_frame.hxx"

static void endpoints()
{
  using namespace gp_StoredRepresentation;
  Arithmetic();
  for (const auto& test : THE_CASES)
  {
    const double a=number(test.A), b=number(test.B);
    const auto end = [&](bool up) { return test.Multiply ? MulEnd(a,b,up) : AddEnd(a,b,up); };
    if (test.Refuse)
    {
      // Either endpoint may be finite when the opposite outward endpoint overflows.
      bool rejected=false;
      try { (void)end(false); (void)end(true); }
      catch (const Standard_Failure&) { rejected=true; }
      require(rejected,test.Name);
      continue;
    }
    const double lo=end(false), hi=end(true), r=number(test.Rounded);
    require(std::isfinite(lo) && std::isfinite(hi),test.Name);
    require(lo <= number(test.Floor) && hi >= number(test.Ceil),test.Name);
    require(lo >= std::nextafter(r,-INFINITY) && hi <= std::nextafter(r,INFINITY),test.Name);
  }
  for (double a : {0.0, -0.0})
    for (double b : {0.0, -0.0, 1.0, -1.0, number(1), -number(1)})
      for (bool up : {false, true})
      {
        const std::uint64_t expected = (bits(a) ^ bits(b)) & 0x8000000000000000ULL;
        require(bits(MulEnd(a,b,up))==expected,"genuine zero product sign");
        require(bits(MulEnd(b,a,up))==expected,"genuine zero product sign reversed");
      }
  for (double a : {0x1p-1022, -0x1p-1022})
  {
    const double lo=MulEnd(a,0x1p-1022,false), hi=MulEnd(a,0x1p-1022,true);
    require(bits(lo)==0x8000000000000001ULL && bits(hi)==1,
            "nonzero rounded-zero product outward endpoints");
  }
  require(bits(AddEnd(-0.0,-0.0,true))==0x8000000000000000ULL,"negative zero sum");
  require(bits(AddEnd(0.0,-0.0,true))==0,"opposite zero sum");
  require(bits(AddEnd(number(1),0.0,false))==1,"subnormal endpoint zero identity");
  refused([] { Point(number(1)); },"subnormal domain");
  refused([] { Point(INFINITY); },"nonfinite");
  refused([] { Point(number(0x7ff8000000000001ULL)); },"nonfinite");
  refused([] { MulEnd(INFINITY,0.0,false); },"nonfinite");
  refused([] { AddEnd(INFINITY,-INFINITY,false); },"nonfinite");
  refused([] { Unit(gp_XYZ(0.0,0.0,0.0)); },"unit domain");
  refused([] { Unit(gp_XYZ(1.0,0.01,0.0)); },"unit domain");
  refused([] { Orthogonal(gp_XYZ(1,0,0),gp_XYZ(1,0,0)); },"orthogonality domain");
  refused([] { RightHandedFrame(gp_XYZ(0,0,1),gp_XYZ(1,0,0),gp_XYZ(0,-1,0)); },"negative or unresolved determinant");
}

static void underflowFrame()
{
  const auto xyz = [](int i) {
    return gp_XYZ(THE_UNDERFLOW_FRAME[i],THE_UNDERFLOW_FRAME[i+1],THE_UNDERFLOW_FRAME[i+2]);
  };
  const gp_Ax3 frame=gp_Ax3::FromStoredFrame(gp_Pnt(xyz(0)),xyz(3),xyz(6),xyz(9));
  const gp_XYZ actual[] = {frame.Location().XYZ(),frame.Direction().XYZ(),
                          frame.XDirection().XYZ(),frame.YDirection().XYZ()};
  for(int i=0;i<12;++i)
    require(bits(actual[i/3].Coord(i%3+1))==bits(THE_UNDERFLOW_FRAME[i]),
            "underflow frame coefficient bits");
}

static void exactFrame()
{
  const std::uint64_t raw[] = {
    0x3fc3333333333333ULL,0x403ba503e3cd0362ULL,0xb93f3ffffffffffcULL,
    0x9f675c4e89d655afULL,0x32ec37186d10a9a4ULL,0x3ff0000000000000ULL,
    0x3ff0000000000000ULL,0x2c6a7e810355e07dULL,0x8000000000000000ULL,
    0xac6a7e810355e07dULL,0x3ff0000000000000ULL,0xb2ec37186d10a9a4ULL};
  const auto xyz = [&](int i) { return gp_XYZ(number(raw[i]),number(raw[i+1]),number(raw[i+2])); };
  const gp_Ax3 frame=gp_Ax3::FromStoredFrame(gp_Pnt(xyz(0)),xyz(3),xyz(6),xyz(9));
  const gp_XYZ actual[] = {frame.Location().XYZ(),frame.Direction().XYZ(),frame.XDirection().XYZ(),frame.YDirection().XYZ()};
  for(int i=0;i<12;++i) require(bits(actual[i/3].Coord(i%3+1))==raw[i],"stored frame coefficient bits");
}

// Only the test changes control modes. Always restore the entire caller environment.
struct Environment
{
  std::fenv_t Saved;
  Environment() { require(std::fegetenv(&Saved)==0,"save environment"); }
  ~Environment() { if (std::fesetenv(&Saved)!=0) std::terminate(); }
};
static void modes()
{
#if defined(__clang__)
#pragma STDC FENV_ACCESS ON
#endif
  using namespace gp_StoredRepresentation;
  std::fenv_t before,after;
  require(std::fegetenv(&before)==0,"initial environment");
  {
    Environment saved;
    require(std::fesetround(FE_DOWNWARD)==0,"set downward");
    refused([] { Unit(gp_XYZ(0,0,1)); },"nearest binary64");
    require(std::fesetenv(&saved.Saved)==0,"restore before flush control");
#if defined(__APPLE__) && defined(__arm64__)
    const auto control=saved.Saved.__fpcr;
    require(Binary64Control(control),"initial control");
    // New architecture controls are classified, never written on an older CPU.
    for (auto flag : {1ULL,2ULL,0x100ULL,0x200ULL,0x400ULL,0x800ULL,0x1000ULL,0x8000ULL,
                      0x00400000ULL,0x00800000ULL})
      require(!Binary64Control(control | flag),"FPCR forbidden control");
    auto flushed=saved.Saved;
    flushed.__fpcr |= __fpcr_flush_to_zero;
    require(std::fesetenv(&flushed)==0,"set supported FZ");
    refused([] { Unit(gp_XYZ(0,0,1)); },"gradual masked");
#elif defined(__linux__) && defined(__x86_64__) && defined(__SSE2__)
    const auto control=_mm_getcsr();
    require(Binary64Control(control),"initial control");
    require(!Binary64Control(control | 0x40),"DAZ control");
    require(!Binary64Control(control | 0x2000),"MXCSR rounding control");
    require(!Binary64Control(control | 0x4000),"MXCSR rounding control");
    for(unsigned flag=0x80;flag<=0x1000;flag<<=1)
      require(!Binary64Control(control & ~flag),"MXCSR trap control");
    _mm_setcsr(control | 0x8000);
    refused([] { Unit(gp_XYZ(0,0,1)); },"gradual masked");
#endif
    require(std::fesetenv(&saved.Saved)==0,"restore before sticky flag");
    const int sticky = FE_INVALID | FE_DIVBYZERO | FE_INEXACT;
    require(std::feraiseexcept(sticky)==0,"set sticky flags");
    Arithmetic();
    require((std::fetestexcept(sticky)&sticky)==sticky,"sticky flags preserved");
  }
  require(std::fegetenv(&after)==0,"final environment");
  require(std::memcmp(&before,&after,sizeof before)==0,"entire environment restored");
}
int main()
{
  try
  {
    underflowFrame(); endpoints(); exactFrame(); modes();
    std::cout << "stored arithmetic: " << sizeof(THE_CASES)/sizeof(THE_CASES[0])
              << " exact dyadic cases; raw frame bits and mode refusals passed\n";
    return 0;
  }
  catch(const Standard_Failure& e) { std::cerr<<e.GetMessageString()<<'\n'; }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; }
  return 1;
}
