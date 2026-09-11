#pragma once
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepLib.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <IntTools_Context.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom2d_Line.hxx>
#include <TopoDS_Face.hxx>
#include <iostream>
#include <cmath>

inline TopoDS_Face makeSliver()
{
  const double         radius = .06, angle = .01, halfHeight = .35;
  Handle(Geom_Surface) s =
    new Geom_CylindricalSurface(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), radius);
  auto         a   = BRepBuilderAPI_MakeVertex(s->Value(0, -halfHeight)).Vertex();
  auto         z   = BRepBuilderAPI_MakeVertex(s->Value(angle, -halfHeight)).Vertex();
  auto         tip = BRepBuilderAPI_MakeVertex(s->Value(angle / 2, halfHeight)).Vertex();
  BRep_Builder b;
  b.UpdateVertex(tip, .0004);
  Handle(Geom2d_Curve) bottom = new Geom2d_Line(gp_Pnt2d(0, -halfHeight), gp_Dir2d(1, 0));
  Handle(Geom2d_Curve) left   = new Geom2d_Line(gp_Pnt2d(0, 0), gp_Dir2d(0, 1));
  Handle(Geom2d_Curve) right  = new Geom2d_Line(gp_Pnt2d(angle, 0), gp_Dir2d(0, 1));
  auto                 eb     = BRepBuilderAPI_MakeEdge(bottom, s, a, z, 0, angle).Edge();
  auto er = BRepBuilderAPI_MakeEdge(right, s, z, tip, -halfHeight, halfHeight).Edge();
  auto el = BRepBuilderAPI_MakeEdge(left, s, a, tip, -halfHeight, halfHeight).Edge();
  el.Reverse();
  BRepBuilderAPI_MakeWire w;
  w.Add(eb);
  w.Add(er);
  w.Add(el);
  auto f = BRepBuilderAPI_MakeFace(s, w.Wire(), true).Face();
  BRepLib::BuildCurves3d(f);
  return f;
}
