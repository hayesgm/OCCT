// Standalone conformance for stored rigid transforms and binary location readers.
#include <BinTools.hxx>
#include <BinTools_IStream.hxx>
#include <BinTools_LocationSet.hxx>
#include <BinTools_OStream.hxx>
#include <BinTools_ShapeReader.hxx>
#include <BinTools_ShapeWriter.hxx>
#include <Standard_Failure.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#ifndef BINTOOLS_STORED_RIGID_TRANSFORMS
  #error Required stored rigid-transform capability is unavailable
#endif

using Matrix            = std::array<double, 12>;
using Policy            = BinTools::GeometryReadPolicy;
constexpr Policy stored = Policy::PreserveStoredLinePlaneCircleSphereRigidTransforms;

static void require(bool value, const char* message)
{
  if (!value)
    throw std::runtime_error(message);
}

static std::uint64_t bits(double value)
{
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

static Matrix coefficients(const gp_Trsf& value)
{
  Matrix result;
  for (int i = 0; i < 12; ++i)
    result[i] = value.Value(i / 4 + 1, i % 4 + 1);
  return result;
}

static bool exact(const Matrix& a, const Matrix& b)
{
  for (int i = 0; i < 12; ++i)
    if (bits(a[i]) != bits(b[i]))
      return false;
  return true;
}

static gp_Trsf factory(const Matrix& m)
{
  return gp_Trsf::FromStoredRigidRepresentation(
    gp_Mat(m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9], m[10]),
    gp_XYZ(m[3], m[7], m[11]));
}

static std::stringstream encoded(const Matrix& m, bool quick, bool compound = false)
{
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  if (quick)
    stream.put(BinTools_ObjectType_SimpleLocation);
  else
  {
    stream << (compound ? "Locations 2\n" : "Locations 1\n");
    stream.put(1);
  }
  for (double value : m)
    BinTools::PutReal(stream, value);
  if (compound)
  {
    stream.put(2);
    BinTools::PutInteger(stream, 1);
    BinTools::PutInteger(stream, 2);
    BinTools::PutInteger(stream, 0);
  }
  stream.seekg(0);
  return stream;
}

static gp_Trsf read(std::stringstream& stream, bool quick, Policy policy, bool compound = false)
{
  BinTools::GeometryReadScope scope(stream, policy);
  if (quick)
  {
    BinTools_IStream     reader(stream);
    BinTools_ShapeReader shapes;
    return shapes.ReadLocation(reader)->Transformation();
  }
  BinTools_LocationSet locations;
  locations.Read(stream);
  return locations.Location(compound ? 2 : 1).Transformation();
}

static gp_Trsf roundtrip(const Matrix& m, bool quick, Policy policy, bool compound = false)
{
  auto stream = encoded(m, quick, compound);
  return read(stream, quick, policy, compound);
}

static void refused(const Matrix& m, bool quick)
{
  auto stream   = encoded(m, quick);
  bool rejected = false;
  try
  {
    (void)read(stream, quick, stored);
  }
  catch (const Standard_Failure&)
  {
    rejected = true;
  }
  require(rejected, "invalid stored transform was accepted");
  require(BinTools::GeometryPolicy(stream) == Policy::Reconstruct, "exception leaked policy");
}

// Encode finite elementary values before any location arithmetic. Separate
// power and multiplication cases make each admission boundary falsifiable.
static std::stringstream arithmetic_record(double translation, bool quick, bool multiply)
{
  const Matrix a = {1, 0, 0, translation, 0, 1, 0, 0, 0, 0, 1, 0};
  Matrix       b = a;
  b[7]           = 1;
  std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
  if (quick)
  {
    stream.put(BinTools_ObjectType_Location);
    stream.put(BinTools_ObjectType_SimpleLocation);
    for (double v : a)
      BinTools::PutReal(stream, v);
    BinTools::PutInteger(stream, multiply ? 1 : 2);
    if (multiply)
    {
      stream.put(BinTools_ObjectType_SimpleLocation);
      for (double v : b)
        BinTools::PutReal(stream, v);
      BinTools::PutInteger(stream, 1);
    }
    stream.put(BinTools_ObjectType_LocationEnd);
  }
  else
  {
    stream << (multiply ? "Locations 3\n" : "Locations 2\n");
    stream.put(1);
    for (double v : a)
      BinTools::PutReal(stream, v);
    if (multiply)
    {
      stream.put(1);
      for (double v : b)
        BinTools::PutReal(stream, v);
    }
    stream.put(2);
    BinTools::PutInteger(stream, 1);
    BinTools::PutInteger(stream, multiply ? 1 : 2);
    if (multiply)
    {
      BinTools::PutInteger(stream, 2);
      BinTools::PutInteger(stream, 1);
    }
    BinTools::PutInteger(stream, 0);
  }
  stream.seekg(0);
  return stream;
}

static gp_Trsf arithmetic_read(std::stringstream& stream, bool quick, bool multiply, Policy policy)
{
  BinTools::GeometryReadScope scope(stream, policy);
  if (quick)
  {
    BinTools_IStream     input(stream);
    BinTools_ShapeReader reader;
    return reader.ReadLocation(input)->Transformation();
  }
  BinTools_LocationSet locations;
  locations.Read(stream);
  return locations.Location(multiply ? 3 : 2).Transformation();
}

int main()
{
  try
  {
    // Generate a generic composed rotation through geometry owners. Independent
    // rounded axis rotations expose reconstruction without captured coefficients.
    gp_Trsf firstRotation, secondRotation;
    firstRotation.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 2, 3)), 0.7);
    secondRotation.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(3, 1, 2)), 1.1);
    const Matrix skew       = coefficients(firstRotation.Multiplied(secondRotation));
    const Matrix identity   = {1, 0, -0.0, 0, -0.0, 1, 0, -0.0, 0, -0.0, 1, 0};
    const Matrix reflection = {-1, 0, -0.0, 2, -0.0, 1, 0, 3, 0, -0.0, 1, 4};
    for (const Matrix& m : {skew, identity, reflection})
    {
      require(exact(m, coefficients(factory(m))), "factory changed admitted bits");
      for (bool quick : {false, true})
        require(exact(m, coefficients(roundtrip(m, quick, stored))),
                "reader changed admitted bits");
    }
    require(factory(reflection).IsNegative(), "reflection lost negative scale");
    require(!factory(skew).IsNegative(), "rotation gained negative scale");
    require(factory(identity).Form() == gp_CompoundTrsf, "form shortcut loses stored coefficients");
    std::cout << "exact_reader_and_representation: pass\n";

    for (bool quick : {false, true})
    {
      const Matrix before = coefficients(roundtrip(skew, quick, Policy::Reconstruct));
      require(!exact(skew, before), "fixture does not distinguish reconstruction");
      require(
        exact(before,
              coefficients(roundtrip(skew, quick, Policy::PreserveStoredLinePlaneCircleSphere))),
        "old policy transform behavior changed");
      (void)roundtrip(skew, quick, stored);
      require(exact(before, coefficients(roundtrip(skew, quick, Policy::Reconstruct))),
              "default changed after scoped read");
      Matrix scaled = {2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0};
      require(roundtrip(scaled, quick, Policy::Reconstruct).ScaleFactor() == 2,
              "default scale lost");
      require(roundtrip(scaled, quick, Policy::PreserveStoredLinePlaneCircleSphere).ScaleFactor()
                == 2,
              "old policy scale lost");
      refused(scaled, quick);
    }
    std::cout << "default_old_policy_and_nonunit_refusal: pass\n";

    for (bool quick : {false, true})
    {
      for (double bad : {std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::denorm_min(),
                         std::numeric_limits<double>::max()})
      {
        Matrix m = identity;
        m[0]     = bad;
        refused(m, quick);
      }
      Matrix shear = identity;
      shear[1]     = 0.25;
      refused(shear, quick);
      Matrix singular = identity;
      singular[0]     = 0;
      refused(singular, quick);
      Matrix translation = identity;
      translation[3]     = std::numeric_limits<double>::quiet_NaN();
      refused(translation, quick);
      for (int count = 0; count < 96; ++count)
      {
        auto complete = encoded(identity, quick);
        auto bytes    = complete.str();
        bytes.resize(bytes.size() - 96 + count);
        std::stringstream truncated(bytes, std::ios::in | std::ios::binary);
        bool              rejected = false;
        try
        {
          (void)read(truncated, quick, stored);
        }
        catch (const Standard_Failure&)
        {
          rejected = true;
        }
        require(rejected, "truncated transform was accepted");
        require(BinTools::GeometryPolicy(truncated) == Policy::Reconstruct,
                "truncation leaked policy");
      }
    }
    std::cout << "malformed_and_all_scalar_truncations: pass\n";

    std::stringstream first, second;
    {
      BinTools::GeometryReadScope outer(first, stored);
      {
        BinTools::GeometryReadScope nested(first, Policy::Reconstruct);
      }
      require(BinTools::GeometryPolicy(first) == stored, "nested scope did not restore");
      require(BinTools::GeometryPolicy(second) == Policy::Reconstruct,
              "policy leaked across streams");
    }
    require(BinTools::GeometryPolicy(first) == Policy::Reconstruct, "outer scope did not restore");
    const int rounding = std::fegetround();
    require(std::fesetround(FE_DOWNWARD) == 0, "rounding unavailable");
    bool rejected = false;
    try
    {
      (void)factory(identity);
    }
    catch (const Standard_Failure&)
    {
      rejected = true;
    }
    std::fesetround(rounding);
    require(rejected, "unsupported rounding accepted");
    std::cout << "scope_and_rounding_refusal: pass\n";

    const gp_Trsf mirror = factory(reflection);
    require(exact(reflection, coefficients(mirror.Powered(1))), "power1 altered representation");
    require(mirror.Powered(0).Form() == gp_Identity, "power0 semantics changed");
    const gp_Trsf inverse = mirror.Inverted();
    const gp_XYZ  point(5, 7, 9);
    gp_XYZ        transformed = point;
    mirror.Transforms(transformed);
    inverse.Transforms(transformed);
    require(transformed.X() == point.X() && transformed.Y() == point.Y()
              && transformed.Z() == point.Z(),
            "reflection inverse algebra failed");
    gp_XYZ powered = point;
    mirror.Powered(2).Transforms(powered);
    gp_XYZ repeated = point;
    mirror.Transforms(repeated);
    mirror.Transforms(repeated);
    require(powered.X() == repeated.X() && powered.Y() == repeated.Y()
              && powered.Z() == repeated.Z(),
            "power2 differs from composition");
    require(exact(coefficients(factory(skew).Powered(2)),
                  coefficients(roundtrip(skew, false, stored, true))),
            "ordinary referenced power2 chain changed");
    // Preserve datum structure, powers and serialized references through the
    // actual quick writer/reader; do not flatten a compound to its matrix.
    for (int power : {-1, 1, 2})
    {
      const TopLoc_Location a(factory(skew)), b(factory(reflection));
      const TopLoc_Location chain = a.Powered(power) * b;
      require(!chain.NextLocation().IsIdentity(), "compound fixture flattened");
      std::stringstream    stream(std::ios::in | std::ios::out | std::ios::binary);
      BinTools_OStream     output(stream);
      BinTools_ShapeWriter writer;
      writer.WriteLocation(output, chain);
      const std::uint64_t referencePosition = output.Position();
      writer.WriteLocation(output, chain);
      require(output.Position() - referencePosition < 96, "repeat was not a reference");
      stream.seekg(0);
      BinTools::GeometryReadScope scope(stream, stored);
      BinTools_IStream            input(stream);
      BinTools_ShapeReader        reader;
      const TopLoc_Location*      decoded = reader.ReadLocation(input);
      require(exact(coefficients(chain.Transformation()), coefficients(decoded->Transformation())),
              "quick compound arithmetic changed");
      TopLoc_Location expectedDatum = chain, actualDatum = *decoded;
      for (; !expectedDatum.IsIdentity();
           expectedDatum = expectedDatum.NextLocation(), actualDatum = actualDatum.NextLocation())
      {
        require(!actualDatum.IsIdentity(), "quick compound lost datum");
        require(expectedDatum.FirstPower() == actualDatum.FirstPower(),
                "quick compound changed power");
        require(exact(coefficients(expectedDatum.FirstDatum()->Transformation()),
                      coefficients(actualDatum.FirstDatum()->Transformation())),
                "quick compound changed elementary bits");
      }
      require(actualDatum.IsIdentity(), "quick compound added datum");
      require(reader.ReadLocation(input) == decoded, "cached location reference was not reused");
      stream.seekg(referencePosition);
      BinTools_IStream       uncachedInput(stream);
      BinTools_ShapeReader   uncachedReader;
      const TopLoc_Location* uncached = uncachedReader.ReadLocation(uncachedInput);
      require(exact(coefficients(chain.Transformation()), coefficients(uncached->Transformation())),
              "uncached location reference changed compound");
      require(uncachedInput.Position() == output.Position(),
              "reference read did not restore position");
    }
    const Matrix reread = coefficients(roundtrip(skew, true, stored));
    require(exact(skew, coefficients(roundtrip(reread, true, stored))),
            "second preserving read drifted");
    std::cout << "inverse_power_compound_and_repeated_read: pass\n";
    for (bool quick : {false, true})
      for (bool multiply : {false, true})
      {
        const double maximum = std::numeric_limits<double>::max();
        for (int policy = 0; policy <= 4; ++policy)
        {
          auto       finite = arithmetic_record(maximum / 8, quick, multiply);
          const auto result = arithmetic_read(finite, quick, multiply, static_cast<Policy>(policy));
          require(result.Value(1, 4) == maximum / 4, "finite location arithmetic changed");
          require(result.Value(2, 4) == (multiply ? 1 : 0), "finite compound datum lost");
          auto overflow = arithmetic_record(maximum, quick, multiply);
          if (policy < 4)
          {
            const auto old =
              arithmetic_read(overflow, quick, multiply, static_cast<Policy>(policy));
            require(!std::isfinite(old.Value(1, 4)), "old policy arithmetic behavior changed");
          }
          else
          {
            bool refusedOverflow = false;
            try
            {
              (void)arithmetic_read(overflow, quick, multiply, stored);
            }
            catch (const Standard_Failure& error)
            {
              refusedOverflow =
                std::string(error.GetMessageString()).find("nonfinite stored location arithmetic")
                != std::string::npos;
            }
            require(refusedOverflow, "nonfinite computed location escaped or wrong refusal");
            require(BinTools::GeometryPolicy(overflow) == Policy::Reconstruct,
                    "arithmetic refusal leaked scoped policy");
          }
        }
      }
    std::cout << "powered_and_composed_finite_arithmetic: pass\n";
    std::cout << "tests=6 passed=6 failed=0\n";
    return 0;
  }
  catch (const Standard_Failure& e)
  {
    std::cerr << "native refusal: " << e.GetMessageString() << '\n';
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }
  return 1;
}
