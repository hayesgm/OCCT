// Generated geometric controls for classifier projection and singular boundaries.
#include <IntTools_FClass2dProjection.pxx>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <BRepTools.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <ProjLib_HCompProjectedCurve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom2d_Line.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <CSLib_NormalPolyDef.hxx>

static void require(bool value, const char* message)
{
  if (!value)
    throw std::runtime_error(message);
}

static void projection(const gp_Pnt2d& a, const gp_Pnt2d& m, const gp_Pnt2d& b, double u, double v)
{
  const gp_Pnt2d result = IntTools_FClass2dProjection(a, m, b);
  require(std::abs(result.X() - m.X()) == u && std::abs(result.Y() - m.Y()) == v,
          "classifier projection deviation changed");
}

static void controls()
{
  projection({0, 0}, {1, 0}, {2, 0}, 0, 0);
  projection({0, 0}, {1, 2}, {2, 0}, 0, 2);
  projection({1, 2}, {1, 2}, {1, 2}, 0, 0);
  projection({0, 0}, {3, 4}, {0, 0}, 3, 4);
  projection({0, 0}, {1.e-12, 2.e-12}, {2.e-12, 0}, 0, 2.e-12);
  projection({0, 0}, {3.e-200, 4.e-200}, {1.e-200, 0}, 3.e-200, 4.e-200);
  std::cout << "classifier_projection_controls=6 passed=6\n";
  for (int degree : {0, 1, 2})
  {
    TColStd_Array1OfReal coefficients(0, degree);
    coefficients(0) = 2.;
    if (degree >= 1)
      coefficients(1) = -3.;
    if (degree == 2)
      coefficients(2) = 5.;
    CSLib_NormalPolyDef polynomial(degree, coefficients);
    for (double angle :
         {0., std::acos(-1.0) / 2, std::acos(-1.0), 3 * std::acos(-1.0) / 2, .37, -.91})
    {
      const double c = std::cos(angle), s = std::sin(angle);
      const double expected   = degree == 0   ? 2.
                                : degree == 1 ? 2 * s - 3 * c
                                              : 2 * s * s - 6 * c * s + 5 * c * c;
      const double derivative = degree == 0   ? 0.
                                : degree == 1 ? 2 * c + 3 * s
                                              : -6 * s * c - 6 * (c * c - s * s);
      double       v, d, v2, d2;
      polynomial.Value(angle, v);
      polynomial.Derivative(angle, d);
      polynomial.Values(angle, v2, d2);
      require(std::isfinite(v) && std::isfinite(d) && std::abs(v - expected) < 1e-12
                && std::abs(d - derivative) < 1e-12 && std::abs(v - v2) < 1e-12
                && std::abs(d - d2) < 1e-12,
              "normal polynomial identity changed");
    }
    std::cout << "CONTROL name=normal_polynomial degree=" << degree << " angles=6 pass=1\n";
  }
  Handle(Geom_SurfaceOfLinearExtrusion) plane =
    new Geom_SurfaceOfLinearExtrusion(new Geom_Line(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))),
                                      gp_Dir(0, 1, 0));
  Handle(GeomAdaptor_Surface) flat = new GeomAdaptor_Surface(plane, -1, 1, -1, 1);
  for (const bool point : {false, true})
  {
    Handle(Geom_Curve) c =
      point ? Handle(Geom_Curve)(new Geom_Line(gp_Ax1(gp_Pnt(.25, -.1, 0), gp_Dir(0, 0, 1))))
            : Handle(Geom_Curve)(new Geom_Line(gp_Ax1(gp_Pnt(0, .2, 0), gp_Dir(1, 0, 0))));
    const double                first = point ? .001 : -.5, last = point ? .002 : .5;
    Handle(GeomAdaptor_Curve)   curve = new GeomAdaptor_Curve(c, first, last);
    ProjLib_HCompProjectedCurve projected(flat, curve, 1e-7, 1e-7, .01);
    gp_Pnt2d                    single;
    require(projected.NbCurves() == 1 && projected.IsSinglePnt(1, single) == point,
            "point classification changed");
    for (int k = 0; k <= 16; ++k)
    {
      const double t  = first + (last - first) * k / 16.;
      const auto   uv = projected.Value(t);
      require(std::abs(uv.X() - (point ? .25 : t)) < 1e-7
                && std::abs(uv.Y() - (point ? -.1 : .2)) < 1e-7,
              "plane projection changed");
    }
    std::cout << "CONTROL name=" << (point ? "genuine_point" : "nonsingular_extrusion")
              << " pass=1 parts=" << projected.NbCurves() << " point=" << point << '\n';
  }
  Handle(Geom_Circle) circle = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1);
  Handle(Geom_SurfaceOfLinearExtrusion) cylinder =
    new Geom_SurfaceOfLinearExtrusion(circle, gp_Dir(0, 0, 1));
  Handle(GeomAdaptor_Surface) periodic =
    new GeomAdaptor_Surface(cylinder, 0, 2 * std::acos(-1.0), -1, 1);
  Handle(Geom_Circle)       shifted = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, .3), gp_Dir(0, 0, 1)), 1);
  Handle(GeomAdaptor_Curve) arc     = new GeomAdaptor_Curve(shifted, 6.1, 6.4);
  ProjLib_HCompProjectedCurve projected(periodic, arc, 1e-7, 1e-7, 1e-4);
  require(projected.NbCurves() > 0, "periodic projection missing");
  for (int k = 0; k <= 32; ++k)
  {
    const double t  = 6.1 + .3 * k / 32.;
    const auto   uv = projected.Value(t);
    require(std::abs(std::remainder(uv.X() - t, 2 * std::acos(-1.0))) < 1e-7
              && std::abs(uv.Y() - .3) < 1e-7,
            "periodic projection branch changed");
    require(periodic->Value(uv.X(), uv.Y()).Distance(arc->Value(t)) < 1e-7,
            "periodic 3D projection changed");
  }
  std::cout << "CONTROL name=periodic_extrusion_across_seam pass=1 parts=" << projected.NbCurves()
            << '\n';
}

static void normal_controls(bool interior)
{
  // Exact quadratic C(u)=((u-s)^2,u,0), extruded in +Y. At u=s,
  // N=0 exactly and N_u=(0,0,2): inward lower limit +Z, upper -Z.
  for (const double singular :
       interior ? std::initializer_list<double>{.5} : std::initializer_list<double>{0., 1.})
  {
    TColgp_Array1OfPnt poles(1, 3);
    poles(1) = gp_Pnt(singular * singular, 0, 0);
    poles(2) = gp_Pnt(singular * singular - singular, .5, 0);
    poles(3) = gp_Pnt((1 - singular) * (1 - singular), 1, 0);
    Handle(Geom_SurfaceOfLinearExtrusion) folded =
      new Geom_SurfaceOfLinearExtrusion(new Geom_BezierCurve(poles), gp_Dir(0, 1, 0));
    const auto          face = BRepBuilderAPI_MakeFace(folded, 0, 1, -1, 1, 1e-7).Face();
    Handle(Geom2d_Line) pc   = new Geom2d_Line(gp_Pnt2d(singular, -.5), gp_Dir2d(0, 1));
    const auto          edge = BRepBuilderAPI_MakeEdge(pc, folded, 0, 1).Edge();
    gp_Pnt              p;
    gp_Vec              du, dv, duu, dvv, duv;
    folded->D2(singular, 0, p, du, dv, duu, dvv, duv);
    const double cross = gp_Dir(du).XYZ().Crossed(gp_Dir(dv).XYZ()).Modulus();
    require(cross == 0 && cross <= gp::Resolution(),
            "analytic fixture does not hit exact singular trigger");
    double lo, hi, bottom, top;
    BRepTools::UVBounds(face, lo, hi, bottom, top);
    require(lo == 0 && hi == 1 && bottom == -1 && top == 1, "analytic fixture bounds changed");
    if (singular == .5)
    {
      bool refused = false;
      try
      {
        gp_Dir normal;
        BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(edge,
                                                    face,
                                                    .5,
                                                    normal,
                                                    Handle(IntTools_Context)());
      }
      catch (const Standard_ConstructionError&)
      {
        refused = true;
      }
      require(refused, "interior ambiguous normal was accepted");
      std::cout << "CONTROL name=analytic_interior_ambiguous_normal_refused cross=" << cross
                << " pass=1\n";
    }
    else
    {
      gp_Dir normal;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(edge,
                                                  face,
                                                  .5,
                                                  normal,
                                                  Handle(IntTools_Context)());
      const gp_Dir expected(0, 0, singular == 0 ? 1 : -1);
      require(normal.Dot(expected) > 1 - 1e-12, "analytic one-sided normal orientation changed");
      std::cout << "CONTROL name="
                << (singular == 0 ? "analytic_lower_normal" : "analytic_upper_normal")
                << " cross=" << cross << " normal=" << normal.X() << ',' << normal.Y() << ','
                << normal.Z() << " expected_z=" << expected.Z() << " pass=1\n";
    }
  }
}

int main(int argc, char** argv)
{
  try
  {
    require(argc == 1 || (argc == 2 && std::string(argv[1]) == "--interior-normal"),
            "usage: projection_controls [--interior-normal]");
    if (argc == 2)
    {
      normal_controls(true);
      return 0;
    }
    controls();
    normal_controls(false);
    std::cout << "projection_and_boundary_controls=pass\n";
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
