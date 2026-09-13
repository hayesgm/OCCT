#include "sliver_fixture.hxx"
#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS.hxx>

int main()
{
  constexpr bool           repaired = true;
  auto                     a = makeSliver(), b = makeSliver();
  Handle(IntTools_Context) c = new IntTools_Context;
  bool valid                 = BRepCheck_Analyzer(a).IsValid() && BRepCheck_Analyzer(b).IsValid();
  std::cout << "valid=" << valid << " distinct=" << !a.IsSame(b) << std::endl;
  if (!valid || a.IsSame(b))
    return 1;
  gp_Trsf tr;
  tr.SetTranslation(gp_Vec(1, 0, 0));
  auto moved       = TopoDS::Face(BRepBuilderAPI_Transform(b, tr, true).Shape());
  bool separatedAB = BOPTools_AlgoTools::AreFacesSameDomain(a, moved, c, 1e-7),
       separatedBA = BOPTools_AlgoTools::AreFacesSameDomain(moved, a, c, 1e-7);
  std::cout << "separated=" << separatedAB << "," << separatedBA << std::endl;
  if (separatedAB || separatedBA)
    return 4;
  for (auto f : {a, b})
  {
    gp_Pnt   p;
    gp_Pnt2d uv;
    int      err  = BOPTools_AlgoTools3D::PointInFace(f, p, uv, c);
    bool     same = BOPTools_AlgoTools::AreFacesSameDomain(f, f.IsSame(a) ? b : a, c, 1e-7);
    std::cout << "err=" << err << " same=" << same << " uv=" << uv.X() << "," << uv.Y()
              << std::endl;
    if ((err == 0) != repaired || same != repaired)
      return 2;
    if (err == 0 && !c->IsValidPointForFace(p, f, 1e-7))
      return 3;
  }
}
