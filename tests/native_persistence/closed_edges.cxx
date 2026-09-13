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
#include <Geom_Circle.hxx>
#include <BRepCheck_Analyzer.hxx>

// Observe ordinary stage results without changing the data structure or order.
class ClosureStageObserver : public BOPAlgo_PaveFiller
{
public:
  TopoDS_Edge Edges[2];
  bool SawOpenRanges = false;
  bool SawAcquiredClosure = false;

protected:
  void PerformVE(const Message_ProgressRange& theRange) override
  {
    BOPAlgo_PaveFiller::PerformVE(theRange);
    const BOPDS_DS& aDS = *PDS();
    int anOpenCount = 0;
    for (const auto& anEdge : Edges)
    {
      const int anIndex = aDS.Index(anEdge);
      if (!aDS.HasPaveBlocks(anIndex))
        continue;
      for (const auto& aPB : aDS.PaveBlocks(anIndex))
      {
        double aFirst, aLast;
        aPB->Range(aFirst, aLast);
        if (aLast - aFirst > 6. && aPB->Pave1().Index() != aPB->Pave2().Index())
          ++anOpenCount;
      }
    }
    SawOpenRanges |= anOpenCount == 2;
  }

  void PerformEE(const Message_ProgressRange& theRange) override
  {
    // The ordinary caller has just updated existing paves with SD vertices.
    const BOPDS_DS& aDS = *PDS();
    int aClosedCount = 0;
    for (const auto& anEdge : Edges)
    {
      const int anIndex = aDS.Index(anEdge);
      if (!aDS.HasPaveBlocks(anIndex))
        continue;
      for (const auto& aPB : aDS.PaveBlocks(anIndex))
      {
        double aFirst, aLast;
        aPB->Range(aFirst, aLast);
        if (aLast - aFirst > 6. && aPB->Pave1().Index() == aPB->Pave2().Index())
          ++aClosedCount;
      }
    }
    SawAcquiredClosure |= SawOpenRanges && aClosedCount == 2;
    BOPAlgo_PaveFiller::PerformEE(theRange);
  }
};

static bool CheckAcquiredClosure()
{
  for (const bool shouldMerge : {true, false})
  {
    BRep_Builder aBuilder;
    const double aLast = 2 * acos(-1.) - .0004;
    auto aShared = BRepBuilderAPI_MakeVertex(gp_Pnt(1, 0, 0)).Vertex();
    Handle(Geom_Circle) aCircle1 =
      new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.);
    Handle(Geom_Circle) aCircle2 =
      new Geom_Circle(gp_Ax2(gp_Pnt(-.0001, 0, 0), gp_Dir(0, 0, 1)), 1.0001);
    auto anEnd1 = BRepBuilderAPI_MakeVertex(aCircle1->Value(aLast)).Vertex();
    auto anEnd2 = BRepBuilderAPI_MakeVertex(aCircle2->Value(aLast)).Vertex();
    auto anEdge1 = BRepBuilderAPI_MakeEdge(aCircle1, aShared, anEnd1, 0, aLast).Edge();
    auto anEdge2 = BRepBuilderAPI_MakeEdge(aCircle2, aShared, anEnd2, 0, aLast).Edge();
    TopoDS_Compound anArgument, aTool;
    aBuilder.MakeCompound(anArgument);
    aBuilder.Add(anArgument, anEdge1);
    aBuilder.Add(anArgument, anEdge2);
    aBuilder.MakeCompound(aTool);
    Handle(Geom_Circle) aToolCircle =
      new Geom_Circle(gp_Ax2(gp_Pnt(-98.9995, 0, 0), gp_Dir(0, 0, 1)), 100.);
    auto aToolArc = BRepBuilderAPI_MakeEdge(aToolCircle, -.01, .01).Edge();
    aBuilder.UpdateEdge(aToolArc, shouldMerge ? .001 : .0001);
    aBuilder.Add(aTool, aToolArc);
    // These disjoint interior vertices cause ordinary VE broad-phase setup to
    // initialize both target ranges, but neither splits nor coincides with them.
    Handle(Geom_Circle) anInnerCircle =
      new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), .9);
    aBuilder.Add(aTool, BRepBuilderAPI_MakeEdge(anInnerCircle, acos(-1.), acos(-1.) + .1).Edge());
    if (!BRepCheck_Analyzer(anArgument).IsValid() || !BRepCheck_Analyzer(aTool).IsValid())
    {
      std::cerr << "Acquired closure input invalid, merge=" << shouldMerge << std::endl;
      return false;
    }
    ClosureStageObserver aFiller;
    aFiller.Edges[0] = anEdge1;
    aFiller.Edges[1] = anEdge2;
    TopTools_ListOfShape anArgs;
    anArgs.Append(anArgument);
    anArgs.Append(aTool);
    aFiller.SetArguments(anArgs);
    aFiller.SetFuzzyValue(1e-7);
    aFiller.SetRunParallel(false);
    aFiller.SetUseOBB(false);
    aFiller.SetNonDestructive(true);
    aFiller.Perform();
    const auto aDS = aFiller.PDS();
    const int n1 = aDS->Index(anEdge1), n2 = aDS->Index(anEdge2);
    int aSharedCount = 0;
    for (TColStd_ListIteratorOfListOfInteger i(aDS->ShapeInfo(n1).SubShapes()); i.More(); i.Next())
      for (TColStd_ListIteratorOfListOfInteger j(aDS->ShapeInfo(n2).SubShapes()); j.More(); j.Next())
        aSharedCount += i.Value() == j.Value();
    bool hasCommon = false;
    int aClosedCount = 0;
    for (int n : {n1, n2})
    {
      if (aDS->ShapeInfo(n).SubShapes().Extent() != 2
          || (shouldMerge && !aDS->HasPaveBlocks(n))
          || (aDS->HasPaveBlocks(n) && aDS->PaveBlocks(n).Extent() != 1))
      {
        std::cerr << "Acquired closure range inventory changed, merge=" << shouldMerge << std::endl;
        return false;
      }
      // A fully separated tool may leave the target's lazy pave list absent.
      if (!aDS->HasPaveBlocks(n))
        continue;
      for (const auto& aPB : aDS->PaveBlocks(n))
      {
        double aFirst, aEnd;
        aPB->Range(aFirst, aEnd);
        if (aFirst != 0 || aEnd != aLast)
        {
          std::cerr << "Acquired closure range changed, merge=" << shouldMerge
                    << " actual=" << aFirst << "," << aEnd << " expected=" << aLast << std::endl;
          return false;
        }
        const int aV1 = aPB->Pave1().Index(), aV2 = aPB->Pave2().Index();
        aClosedCount += aV1 == aV2 && aDS->IsNewShape(aV1);
        auto aCB = aDS->CommonBlock(aPB);
        if (!aCB.IsNull())
          for (const auto& anOther : aCB->PaveBlocks())
            hasCommon |= anOther->OriginalEdge() == (n == n1 ? n2 : n1);
      }
    }
    std::cout << "ACQUIRED_CLOSURE merge=" << shouldMerge
              << " original_shared=" << aSharedCount << " open_ranges=" << aFiller.SawOpenRanges
              << " acquired_closure=" << aFiller.SawAcquiredClosure
              << " final_closed=" << aClosedCount << " common=" << hasCommon
              << " ee=" << aDS->InterfEE().Length() << std::endl;
    if (aFiller.HasErrors() || aSharedCount != 1 || aDS->Rank(n1) != aDS->Rank(n2)
        || (shouldMerge && !aFiller.SawOpenRanges) || aFiller.SawAcquiredClosure != shouldMerge
        || aClosedCount != (shouldMerge ? 2 : 0) || hasCommon != shouldMerge
        || aDS->InterfEE().Length() != (shouldMerge ? 1 : 0))
      return false;
  }
  return true;
}

int main()
{
  constexpr bool patched = true;
  std::cout << std::setprecision(17);
  for (int variant = 0; variant < 2; ++variant)
  {
    double              d = .0001;
    BRep_Builder        b;
    auto                v  = BRepBuilderAPI_MakeVertex(gp_Pnt(1, 0, 0)).Vertex();
    auto                w  = variant == 0 ? v : BRepBuilderAPI_MakeVertex(gp_Pnt(1, 0, 0)).Vertex();
    Handle(Geom_Circle) c1 = new Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 1.);
    Handle(Geom_Circle) c2 = new Geom_Circle(gp_Ax2(gp_Pnt(-d, 0, 0), gp_Dir(0, 0, 1)), 1. + d);
    auto                e1 = BRepBuilderAPI_MakeEdge(c1, v, v, 0, 2 * acos(-1.)).Edge();
    auto                e2 = BRepBuilderAPI_MakeEdge(c2, w, w, 0, 2 * acos(-1.)).Edge();
    TopoDS_Compound     aa, bb;
    b.MakeCompound(aa);
    b.Add(aa, e1);
    b.Add(aa, e2);
    b.MakeCompound(bb);
    auto anchor = BRepBuilderAPI_MakeVertex(gp_Pnt(1, 0, 0)).Vertex();
    b.UpdateVertex(anchor, .001);
    b.Add(bb, anchor);
    std::cout << "valid_input=" << BRepCheck_Analyzer(aa).IsValid() << std::endl;
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
    bool inputValid = BRepCheck_Analyzer(aa).IsValid();
    bool allNew     = true;
    for (int cv : canon[0])
      allNew &= ds->IsNewShape(cv);
    bool closed = true;
    for (int n : {n1, n2})
      for (const auto& pb : ds->PaveBlocks(n))
        closed &= pb->Pave1().Index() == pb->Pave2().Index();
    bool expected = !patched || variant == 1;
    if (!inputValid || !allNew || !closed || ds->Rank(n1) != ds->Rank(n2)
        || ds->InterfEE().Length() != (expected ? 1 : 0) || pf.HasWarnings() != expected
        || pf.HasErrors() || shared != (variant == 0 ? 1 : 0) || sharedCanon != 1
        || acquired != (!patched || variant == 1))
      return 2;
  }
  return CheckAcquiredClosure() ? 0 : 3;
}
