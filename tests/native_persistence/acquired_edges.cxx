#include <BOPAlgo_PaveFiller.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_CommonBlock.hxx>
#include <BOPDS_PaveBlock.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TColStd_ListIteratorOfListOfInteger.hxx>
#include <iostream>
#include <iomanip>
#include <set>
#include <GC_MakeArcOfCircle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <BRepCheck_Analyzer.hxx>

int main()
{
  std::cout << std::setprecision(17);
  for (int variant = 0; variant < 4; ++variant)
  {
    double          d = variant == 2 ? .01 : .0001;
    BRep_Builder    b;
    auto            a  = BRepBuilderAPI_MakeVertex(gp_Pnt(0, 0, 0)).Vertex();
    auto            z  = BRepBuilderAPI_MakeVertex(gp_Pnt(1, 0, 0)).Vertex();
    auto            c  = (variant == 1 || variant == 3)
                           ? a
                           : BRepBuilderAPI_MakeVertex(gp_Pnt(0, d, 0)).Vertex();
    auto            w  = variant == 3 ? z : BRepBuilderAPI_MakeVertex(gp_Pnt(1, d, 0)).Vertex();
    auto            e1 = BRepBuilderAPI_MakeEdge(a, z).Edge();
    auto            e2 = BRepBuilderAPI_MakeEdge(c, w).Edge();
    if (variant == 3)
    {
      Handle(Geom_TrimmedCurve) curve1 =
        GC_MakeArcOfCircle(gp_Pnt(0, 0, 0), gp_Pnt(.5, .00001, 0), gp_Pnt(1, 0, 0));
      Handle(Geom_TrimmedCurve) curve2 =
        GC_MakeArcOfCircle(gp_Pnt(0, 0, 0), gp_Pnt(.5, -.00001, 0), gp_Pnt(1, 0, 0));
      e1 = BRepBuilderAPI_MakeEdge(curve1, a, z).Edge();
      e2 = BRepBuilderAPI_MakeEdge(curve2, c, w).Edge();
    }
    if (!BRepCheck_Analyzer(e1).IsValid() || !BRepCheck_Analyzer(e2).IsValid())
      return 4;
    TopoDS_Compound aa, bb;
    b.MakeCompound(aa);
    b.Add(aa, e1);
    b.Add(aa, e2);
    b.MakeCompound(bb);
    for (double x : {0., 1.})
    {
      auto v = BRepBuilderAPI_MakeVertex(
                 gp_Pnt(x, variant == 3 || (variant == 1 && x == 0) ? 0 : d / 2, 0))
                 .Vertex();
      b.UpdateVertex(v, .0001);
      b.Add(bb, v);
    }
    TopTools_ListOfShape args;
    args.Append(aa);
    args.Append(bb);
    BOPAlgo_PaveFiller pf;
    pf.SetArguments(args);
    pf.SetFuzzyValue(1e-7);
    pf.SetRunParallel(false);
    pf.SetUseOBB(false);
    pf.SetNonDestructive(true);
    std::cout << "CASE " << variant << std::endl;
    pf.Perform();
    auto ds = pf.PDS();
    int  n1 = ds->Index(e1), n2 = ds->Index(e2);
    std::cout << "errors=" << pf.HasErrors() << " warnings=" << pf.HasWarnings()
              << " indices=" << n1 << "," << n2 << " ranks=" << ds->Rank(n1) << "," << ds->Rank(n2)
              << " interf_ee=" << ds->InterfEE().Length() << std::endl;
    pf.DumpWarnings(std::cout);
    std::set<int> orig[2], canon[2];
    int           i = 0;
    for (int n : {n1, n2})
    {
      for (TColStd_ListIteratorOfListOfInteger it(ds->ShapeInfo(n).SubShapes()); it.More();
           it.Next())
      {
        orig[i].insert(it.Value());
        int cv = it.Value();
        ds->HasShapeSD(it.Value(), cv);
        canon[i].insert(cv);
        std::cout << "endpoint edge=" << n << " original=" << it.Value() << " canonical=" << cv
                  << " new=" << ds->IsNewShape(cv) << std::endl;
      }
      ++i;
    }
    int shared = 0;
    for (int n : orig[0])
      shared += orig[1].count(n);
    int sharedCanon = 0;
    for (int n : canon[0])
      sharedCanon += canon[1].count(n);
    bool acquired = false;
    for (const auto& pb : ds->PaveBlocks(n1))
    {
      auto cb = ds->CommonBlock(pb);
      if (cb.IsNull())
        continue;
      for (const auto& other : cb->PaveBlocks())
      {
        if (other->OriginalEdge() == n2)
          acquired = true;
      }
    }
    std::cout << "RESULT original_shared=" << shared << " canonical_shared=" << sharedCanon
              << " same_rank_common_block=" << acquired << std::endl;
    bool allNew = true;
    for (int cv : canon[0])
      allNew &= ds->IsNewShape(cv);
    if (pf.HasErrors() || ds->Rank(n1) != ds->Rank(n2)
        || ds->InterfEE().Length() != (variant < 2 ? 1 : 0) || sharedCanon != (variant != 2 ? 2 : 0)
        || allNew != (variant != 2) || pf.HasWarnings() != (variant < 2)
        || shared != (variant == 3 ? 2 : variant == 1 ? 1 : 0) || acquired != (variant < 2))
      return 2;
  }
}
