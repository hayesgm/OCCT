#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BOPAlgo_CheckerSI.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_Pair.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <GeomConvert.hxx>
#include <IntTools_FaceFace.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <iostream>
#include <stdexcept>
#include <string>

static void require(bool condition, const char* message)
{
  if (!condition) throw std::runtime_error(message);
}

static TopoDS_Shape readFixture(const char* file)
{
  TopoDS_Shape shape;
  BRep_Builder builder;
  const std::string path = std::string(SPLINE_FIXTURE_DIR) + "/" + file;
  require(BRepTools::Read(shape, path.c_str(), builder), "read fixture");
  require(!shape.IsNull() && BRepCheck_Analyzer(shape, true).IsValid(), "valid fixture");
  return shape;
}

static void requireEmptyCut(const TopoDS_Shape& a, const TopoDS_Shape& b)
{
  BRepAlgoAPI_Cut cut;
  TopTools_ListOfShape arguments, tools;
  arguments.Append(a); tools.Append(b);
  cut.SetArguments(arguments); cut.SetTools(tools);
  cut.SetFuzzyValue(1e-7); cut.SetRunParallel(false); cut.SetUseOBB(false);
  cut.SetNonDestructive(true); cut.SetCheckInverted(true);
  cut.Build();
  require(cut.IsDone() && !cut.HasErrors(), "cut completion");
  const auto result = cut.Shape();
  require(!result.IsNull() && BRepCheck_Analyzer(result, true).IsValid(), "valid cut result");
  for (auto kind : {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID})
    require(!TopExp_Explorer(result, kind).More(), "nonempty roundtrip difference");
}

static Handle(Geom_BSplineSurface) bilinearPatch(const gp_Pnt& a, const gp_Pnt& b,
                                               const gp_Pnt& c, const gp_Pnt& d)
{
  TColgp_Array2OfPnt poles(1,2,1,2);
  poles(1,1)=a; poles(2,1)=b; poles(1,2)=c; poles(2,2)=d;
  TColStd_Array1OfReal knots(1,2); knots(1)=0; knots(2)=1;
  TColStd_Array1OfInteger mult(1,2); mult.Init(2);
  return new Geom_BSplineSurface(poles,knots,knots,mult,mult,1,1);
}

static void requireIntersections(const TopoDS_Face& a, const TopoDS_Face& b, int expected)
{
  for (bool reverse : {false,true})
  {
    IntTools_FaceFace ff;
    ff.SetParameters(true,true,true,1e-7); ff.SetFuzzyValue(1e-7);
    ff.Perform(reverse?b:a,reverse?a:b,false);
    std::cout << "intersection reverse=" << reverse << " done=" << ff.IsDone()
              << " curves=" << ff.Lines().Length() << std::endl;
    require(ff.IsDone() && ff.Lines().Length()==expected, "genuine intersection changed");
    if (expected) require(!ff.TangentFaces(), "transverse faces classified coincident");
  }
}

static void crossingControls()
{
  const auto flat=bilinearPatch(gp_Pnt(0,0,0),gp_Pnt(1,0,0),gp_Pnt(0,1,0),gp_Pnt(1,1,0));
  for (bool separated : {false,true})
  {
    const auto other=bilinearPatch(gp_Pnt(0,0,separated?.001:-.001),gp_Pnt(1,0,.001),
                                  gp_Pnt(0,1,separated?.001:-.001),gp_Pnt(1,1,.001));
    requireIntersections(BRepBuilderAPI_MakeFace(flat,1e-7).Face(),
                         BRepBuilderAPI_MakeFace(other,1e-7).Face(), separated?0:1);
  }

  // A genuine intersection ending at a collapsed spherical pole must survive.
  const Handle(Geom_Surface) sphere=new Geom_SphericalSurface(gp_Ax3(),1);
  const Handle(Geom_Surface) octant=new Geom_RectangularTrimmedSurface(sphere,0,M_PI/2,0,M_PI/2);
  const auto curved=BRepBuilderAPI_MakeFace(GeomConvert::SurfaceToBSplineSurface(octant),1e-7).Face();
  const auto plane=bilinearPatch(gp_Pnt(-1,-1,-1),gp_Pnt(2,2,-1),gp_Pnt(-1,-1,2),gp_Pnt(2,2,2));
  requireIntersections(curved,BRepBuilderAPI_MakeFace(plane,1e-7).Face(),1);

  // A cone has a direction-dependent normal at its apex. Its meridian crossing
  // must not be discarded as a coincident neighborhood.
  const Handle(Geom_Surface) cone=new Geom_ConicalSurface(gp_Ax3(),M_PI/4,0);
  const Handle(Geom_Surface) conePart=new Geom_RectangularTrimmedSurface(cone,0.0,M_PI/2,0.0,1.0);
  const auto conical=BRepBuilderAPI_MakeFace(GeomConvert::SurfaceToBSplineSurface(conePart),1e-7).Face();
  requireIntersections(conical,BRepBuilderAPI_MakeFace(plane,1e-7).Face(),1);
}

static void distinctBranchControls()
{
  // A cubic planar loop extruded along Z crosses at u=.2 and u=.8. Each
  // restricted face is regular and injective; equal-parameter closeness of the
  // underlying surfaces does not authorize suppressing this intersection.
  const double x[4]={.16,.16-1./3.,.16-1./3.,.16};
  const double y[4]={-.08,.14,-.14,.08};
  TColgp_Array2OfPnt poles(1,4,1,2);
  for(int i=1;i<=4;++i) for(int j=1;j<=2;++j)
    poles(i,j)=gp_Pnt(x[i-1],y[i-1],j-1);
  TColStd_Array1OfReal knots(1,2); knots(1)=0; knots(2)=1;
  TColStd_Array1OfInteger um(1,2),vm(1,2); um.Init(4); vm.Init(2);
  const Handle(Geom_BSplineSurface) loop=new Geom_BSplineSurface(poles,knots,knots,um,vm,3,1);
  requireIntersections(BRepBuilderAPI_MakeFace(loop,0,.4,0,1,1e-7).Face(),
                       BRepBuilderAPI_MakeFace(loop,.6,1,0,1,1e-7).Face(),1);
  const auto face=BRepBuilderAPI_MakeFace(loop,1e-7).Face();
  BOPAlgo_CheckerSI checker;
  TopTools_ListOfShape arguments; arguments.Append(face);
  checker.SetArguments(arguments); checker.SetLevelOfCheck(5); checker.SetRunParallel(false);
  checker.Perform();
  require(!checker.HasErrors(), "self-intersection checker completion");
  const int index=checker.DS().Index(face);
  require(index>=0 && checker.DS().Interferences().Contains(BOPDS_Pair(index,index)),
          "lost genuine self-intersection");
}

int main(int argc, char** argv)
{
  try
  {
    crossingControls();
    distinctBranchControls();
    if (argc==2 && std::string(argv[1])=="controls") return 0;
    const auto warm=readFixture("warm.brep"), reopened=readFixture("reopened.brep");
    TopTools_IndexedMapOfShape a,b;
    TopExp::MapShapes(warm,TopAbs_FACE,a); TopExp::MapShapes(reopened,TopAbs_FACE,b);
    require(a.Extent()==44 && b.Extent()==44,"fixture topology drift");
    requireEmptyCut(a(8),b(8)); requireEmptyCut(b(8),a(8));
    std::cout << "corner pair passed" << std::endl;
    requireEmptyCut(warm,warm); requireEmptyCut(reopened,reopened);
    requireEmptyCut(warm,reopened); requireEmptyCut(reopened,warm);
    std::cout << "complete solids passed" << std::endl;
    return 0;
  }
  catch(const Standard_Failure& e) { std::cerr<<e.GetMessageString()<<'\n'; }
  catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; }
  return 1;
}
