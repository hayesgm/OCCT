#include "sliver_fixture.hxx"
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <IntTools_Context.hxx>
#include <IntTools_Tools.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Pln.hxx>
#include <iostream>

int main()
{
  constexpr bool patched = true;
  BRep_Builder   b;
  auto           base = makeSliver();
  TopoDS_Face    f    = TopoDS::Face(BRepBuilderAPI_Copy(base).Shape());
  f.Orientation(TopAbs_FORWARD);
  auto                    surf = BRep_Tool::Surface(f);
  BRepBuilderAPI_MakeWire wire;
  gp_Pnt2d                points[] = {{.0075, -.1}, {.0095, -.1}, {.0095, .1}, {.0075, .1}};
  for (int i = 0; i < 4; ++i)
  {
    Handle(Geom2d_Curve) seg = GCE2d_MakeSegment(points[i], points[(i + 1) % 4]).Value();
    wire.Add(BRepBuilderAPI_MakeEdge(seg, surf).Edge());
  }
  auto hole = wire.Wire();
  hole.Reverse();
  b.Add(f, hole);
  BRepLib::BuildCurves3d(f);
  bool faceValid = BRepCheck_Analyzer(f).IsValid();
  std::cout << "holed_face_valid=" << faceValid << std::endl;
  if (!faceValid)
    return 6;
  Handle(IntTools_Context) c = new IntTools_Context;
  double                   u0, u1, v0, v1;
  c->UVBounds(f, u0, u1, v0, v1);
  for (int i = 0; i < 2; ++i)
  {
    double u = IntTools_Tools::IntermediatePoint(u0, u1);
    if (i)
      u = u1 - (u - u0);
    Handle(Geom2d_Line) line = new Geom2d_Line(gp_Pnt2d(u, 0), gp_Dir2d(0, 1));
    gp_Pnt              p;
    gp_Pnt2d            uv;
    int                 err = BOPTools_AlgoTools3D::PointInFace(f, line, p, uv, c);
    std::cout << "original_direction=" << i << " error=" << err << std::endl;
    if (err == 0)
      return 2;
  }
  gp_Pnt   p;
  gp_Pnt2d uv;
  int      err       = BOPTools_AlgoTools3D::PointInFace(f, p, uv, c);
  auto     holeState = c->StatePointFace(f, gp_Pnt2d(.0085, 0));
  std::cout << "fallback_error=" << err << " uv=" << uv.X() << "," << uv.Y()
            << " hole_state=" << holeState << std::endl;
  if (err == 0)
  {
    auto state = c->StatePointFace(f, uv);
    bool valid = c->IsValidPointForFace(p, f, 1e-7);
    std::cout << "returned_state=" << state << " returned_valid=" << valid << std::endl;
    if (state != TopAbs_IN || !valid)
      return 3;
  }
  if ((err == 0) != patched || holeState != TopAbs_OUT)
    return 4;
  auto ordinary =
    BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 0, 1, 0, 1).Face();
  int ordinaryErr = BOPTools_AlgoTools3D::PointInFace(ordinary, p, uv, c);
  std::cout << "ordinary_error=" << ordinaryErr
            << " ordinary_inside=" << c->IsPointInFace(ordinary, uv) << std::endl;
  if (ordinaryErr || !c->IsPointInFace(ordinary, uv))
    return 5;
}
