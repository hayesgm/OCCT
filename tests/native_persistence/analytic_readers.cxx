#include <BinTools.hxx>
#include <BinTools_CurveSet.hxx>
#include <BinTools_SurfaceSet.hxx>
#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Standard_Failure.hxx>
#include <gp_Lin.hxx>
#include <gp_Circ.hxx>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

using Values = std::vector<double>;
using Policy = BinTools::GeometryReadPolicy;

static void require(bool condition, const char* message)
{
  if (!condition)
    throw std::runtime_error(message);
}

static bool exact(const Values& left, const Values& right)
{
  return left.size() == right.size()
         && std::memcmp(left.data(), right.data(), left.size() * sizeof(double)) == 0;
}

static Values payload(int family)
{
  if (family == 0)
    return {0., -0., 3., 1., std::ldexp(1., -40), 0.};
  Values result = {0., -0., 3., 0., 0., 1., 1., 0., 0., std::ldexp(1., -50), 1., 0.};
  if (family == 1 || family == 3)
    result.push_back(2.5);
  return result;
}

static Values decode(int family, const Values& input, Policy policy)
{
  std::stringstream   stream(std::ios::in | std::ios::out | std::ios::binary);
  const unsigned char tags[] = {1, 2, 1, 4};
  stream.put(tags[family]);
  for (double value : input)
    BinTools::PutReal(stream, value);
  stream.seekg(0);
  BinTools::GeometryReadScope scope(stream, policy);
  gp_Ax3                      frame;
  double                      radius = 0.;
  if (family < 2)
  {
    Handle(Geom_Curve) curve;
    BinTools_CurveSet::ReadCurve(stream, curve);
    if (family == 0)
    {
      const gp_Lin line = Handle(Geom_Line)::DownCast(curve)->Lin();
      return {line.Location().X(),
              line.Location().Y(),
              line.Location().Z(),
              line.Direction().X(),
              line.Direction().Y(),
              line.Direction().Z()};
    }
    const gp_Circ circle = Handle(Geom_Circle)::DownCast(curve)->Circ();
    frame                = gp_Ax3(circle.Position());
    radius               = circle.Radius();
  }
  else
  {
    Handle(Geom_Surface) surface;
    BinTools_SurfaceSet::ReadSurface(stream, surface);
    if (family == 2)
      frame = Handle(Geom_Plane)::DownCast(surface)->Position();
    else
    {
      const auto sphere = Handle(Geom_SphericalSurface)::DownCast(surface);
      frame             = sphere->Position();
      radius            = sphere->Radius();
    }
  }
  Values result = {frame.Location().X(),
                   frame.Location().Y(),
                   frame.Location().Z(),
                   frame.Direction().X(),
                   frame.Direction().Y(),
                   frame.Direction().Z(),
                   frame.XDirection().X(),
                   frame.XDirection().Y(),
                   frame.XDirection().Z(),
                   frame.YDirection().X(),
                   frame.YDirection().Y(),
                   frame.YDirection().Z()};
  if (family == 1 || family == 3)
    result.push_back(radius);
  return result;
}

static void refuses(int family, const Values& values, const std::string& reason)
{
  bool refused = false;
  try
  {
    (void)decode(family, values, Policy::PreserveStoredLinePlaneCircleSphere);
  }
  catch (const Standard_Failure& error)
  {
    refused = std::string(error.GetMessageString()).find(reason) != std::string::npos;
    if (!refused)
      std::cerr << "family=" << family << " expected=" << reason
                << " actual=" << error.GetMessageString() << '\n';
  }
  require(refused, "invalid payload escaped or produced an unrelated refusal");
}

int main()
{
  try
  {
    for (int family = 0; family != 4; ++family)
    {
      const Values input = payload(family);
      for (Policy policy : {Policy::PreserveStoredLinePlaneCircleSphere,
                            Policy::PreserveStoredLinePlaneCircleSphereRigidTransforms})
      {
        const Values first = decode(family, input, policy);
        require(exact(first, input), "stored reader changed payload bits");
        require(exact(decode(family, first, policy), input), "repeated stored read changed bits");
      }
      if (family != 0)
        require(!exact(decode(family, input, Policy::Reconstruct), input),
                "generated frame does not distinguish reconstruction");
      Values truncated = input;
      truncated.pop_back();
      refuses(family,
              truncated,
              family == 3 ? "incomplete stored Sphere payload" : "Storage_StreamTypeMismatchError");
      for (int index : {0, 3})
        for (double invalid :
             {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        {
          Values changed = input;
          changed[index] = invalid;
          refuses(family, changed, "nonfinite coefficient");
        }
      Values nonunit = input;
      nonunit[3]     = 2.;
      nonunit[4]     = 0.;
      nonunit[5]     = 0.;
      refuses(family, nonunit, "binary64 unit domain");
      if (family == 1 || family == 3)
      {
        Values negative = input;
        negative.back() = -1.;
        refuses(family,
                negative,
                family == 1 ? "negative stored Circle radius" : "negative stored Sphere radius");
      }
      if (family != 0)
      {
        Values indirect = input;
        for (int index = 9; index != 12; ++index)
          indirect[index] = -indirect[index];
        if (family == 1)
          refuses(family, indirect, "right-handed frame has negative or unresolved determinant");
        else
          require(
            exact(decode(family, indirect, Policy::PreserveStoredLinePlaneCircleSphere), indirect),
            "surface reader changed indirect frame");
      }
      std::cout << "analytic_family=" << family << " exact_repeat_refusal_controls=pass\n";
    }
    const Values nonunit       = {0., 0., 0., 2., 0., 0.};
    const Values reconstructed = {0., 0., 0., 1., 0., 0.};
    require(exact(decode(0, nonunit, Policy::Reconstruct), reconstructed),
            "default direction reconstruction changed");
    refuses(0, nonunit, "binary64 unit domain");
    require(exact(decode(0, nonunit, Policy::Reconstruct), reconstructed),
            "preserving refusal changed subsequent default");
    std::cout << "tests=5 passed=5 failed=0\n";
    return 0;
  }
  catch (const Standard_Failure& error)
  {
    std::cerr << error.GetMessageString() << '\n';
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
  }
  return 1;
}
