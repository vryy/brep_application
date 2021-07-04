// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes

// Project includes
#include "includes/element.h"
#include "containers/array_1d.h"
#include "custom_python3/add_brep_and_level_set_to_python.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/section/section.h"
#include "custom_algebra/section/polygonal_section.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/and_brep.h"
#include "custom_algebra/or_brep.h"
#include "custom_algebra/not_brep.h"
#include "custom_algebra/natm_arc_brep.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_algebra/occ_brep.h"
#endif
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/level_set/nodal_level_set.h"
#include "custom_algebra/level_set/circular_level_set.h"
#include "custom_algebra/level_set/circular2_level_set.h"
#include "custom_algebra/level_set/doughnut_level_set.h"
#include "custom_algebra/level_set/spherical_level_set.h"
#include "custom_algebra/level_set/spherical2_level_set.h"
#include "custom_algebra/level_set/cylinder_level_set.h"
#include "custom_algebra/level_set/cylinder2_level_set.h"
#include "custom_algebra/level_set/cone_level_set.h"
#include "custom_algebra/level_set/linear_level_set.h"
#include "custom_algebra/level_set/planar_level_set.h"
#include "custom_algebra/level_set/product_level_set.h"
#include "custom_algebra/level_set/inverse_level_set.h"
#include "custom_algebra/level_set/union_level_set.h"
#include "custom_algebra/level_set/intersection_level_set.h"
#include "custom_algebra/level_set/difference_level_set.h"
#include "custom_algebra/level_set/closest_level_set.h"
#include "custom_algebra/level_set/distance_to_curve_level_set.h"
#include "custom_algebra/curve/curve.h"
#include "custom_algebra/curve/parametric_curve.h"
#include "custom_algebra/surface/parametric_surface.h"
#include "custom_algebra/volume/parametric_volume.h"


namespace Kratos
{

const int BRep::_CUT;
const int BRep::_IN;
const int BRep::_OUT;

namespace Python
{

using namespace pybind11;

bool BRep_IsInside2(BRep& rDummy, const double& x, const double& y)
{
    BRep::PointType P;
    P[0] = x;
    P[1] = y;
    P[2] = 0.0;
    return rDummy.IsInside(P);
}

bool BRep_IsInside3(BRep& rDummy, const double& x, const double& y, const double& z)
{
    BRep::PointType P;
    P[0] = x;
    P[1] = y;
    P[2] = z;
    return rDummy.IsInside(P);
}

bool BRep_IsInside4Element(BRep& rDummy, Element::Pointer pElement,
    const BRep::CoordinatesArrayType& local_coords, const int& configuration)
{
    return rDummy.IsInside(pElement->GetGeometry(), local_coords, configuration);
}

int BRep_CutStatusPoints(BRep& rDummy, pybind11::list& list_points)
{
    std::vector<BRep::PointType> points;

    for (std::size_t i = 0; i < pybind11::len(list_points); ++i)
        points.push_back(list_points[i].cast<BRep::PointType>());

    return rDummy.CutStatus(points);
}

BRep::PointType BRep_Bisect2(BRep& rDummy, const BRep::PointType& P1, const BRep::PointType& P2)
{
    BRep::PointType P;
    int error_code = rDummy.Bisect(P, P1, P2, rDummy.Tolerance());
    if (error_code == 0)
    {
        return P;
    }
    else
    {
        P[0] = 1.0e99;
        P[1] = 1.0e99;
        P[2] = 1.0e99;
        return P;
    }
}

BRep::PointType BRep_Bisect3(BRep& rDummy, const BRep::PointType& P1, const BRep::PointType& P2, const BRep::PointType& P3)
{
    BRep::PointType P;
    int error_code = rDummy.Bisect(P, P1, P2, P3, rDummy.Tolerance());
    if (error_code == 0)
    {
        return P;
    }
    else
    {
        P[0] = 1.0e99;
        P[1] = 1.0e99;
        P[2] = 1.0e99;
        return P;
    }
}

Section::Pointer BRep_IntersectElement(BRep& rDummy, Element::Pointer pElement)
{
    return rDummy.Intersect(pElement->GetGeometry());
}

BRep::PointType BRep_ProjectOnSurface(BRep& rDummy, const BRep::PointType& rPoint)
{
    BRep::PointType P;
    int error_code = rDummy.ProjectOnSurface(rPoint, P);
    if (error_code == 0)
    {
        return P;
    }
    else
    {
        P[0] = 1.0e99;
        P[1] = 1.0e99;
        P[2] = 1.0e99;
        return P;
    }
}

LevelSet::Pointer InverseLevelSet_GetLevelSet(InverseLevelSet& rDummy)
{
    return rDummy.pLeveSet();
}

void InverseLevelSet_SetLevelSet(InverseLevelSet& rDummy, LevelSet::Pointer p_level_set)
{
    rDummy.pSetLevelSet(p_level_set);
}

template<class TLevelSel>
pybind11::list LevelSet_CreateQ4Elements(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const int& nsampling_axial,
    const int& nsampling_radial,
    const double& start_radial_angle, // in degree
    const double& end_radial_angle) // in degree
{
    const double Pi = 3.1415926535897932384626433;
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4Elements(r_model_part, sample_element_name, pProperties,
            static_cast<std::size_t>(nsampling_axial), static_cast<std::size_t>(nsampling_radial),
            start_radial_angle/180.0*Pi, end_radial_angle/180*Pi);
    pybind11::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
pybind11::list LevelSet_CreateQ4ElementsClosedLoop(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const int& nsampling_axial,
    const int& nsampling_radial)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4ElementsClosedLoop(r_model_part, sample_element_name, pProperties,
            static_cast<std::size_t>(nsampling_axial), static_cast<std::size_t>(nsampling_radial));
    pybind11::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
pybind11::list LevelSet_CreateQ4ElementsClosedLoopWithRange(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const int& nsampling_axial,
    const int& nsampling_radial,
    const double& tmin,
    const double& tmax)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4ElementsClosedLoop(r_model_part, sample_element_name, pProperties,
            static_cast<std::size_t>(nsampling_axial), static_cast<std::size_t>(nsampling_radial), tmin, tmax);
    pybind11::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
pybind11::list LevelSet_CreateQ4ConditionsClosedLoopWithRange(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_condition_name,
    Properties::Pointer pProperties,
    const int& nsampling_axial,
    const int& nsampling_radial,
    const double& tmin,
    const double& tmax,
    const bool& reverse)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> Results
        = rDummy.CreateQ4ConditionsClosedLoop(r_model_part, sample_condition_name, pProperties,
            static_cast<std::size_t>(nsampling_axial), static_cast<std::size_t>(nsampling_radial), tmin, tmax, reverse);
    pybind11::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

pybind11::list LinearLevelSet_CreateLineConditions(
    LinearLevelSet& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_condition_name,
    Properties::Pointer pProperties,
    const BRep::PointType& StartPoint, const BRep::PointType& EndPoint,
    const int& type, const int& nsampling)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> Results
        = rDummy.CreateLineConditions(r_model_part, sample_condition_name, pProperties,
            StartPoint, EndPoint, type, nsampling);
    pybind11::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

void NodalLevelSet_InitializeFromNodes(NodalLevelSet& dummy,
    const ModelPart::NodesContainerType& rNodes, const int& configuration)
{
    dummy.Initialize(rNodes, configuration);
}

void NodalLevelSet_InitializeFromElements(NodalLevelSet& dummy,
    const ModelPart::ElementsContainerType& rElements, const int& configuration)
{
    dummy.Initialize(rElements, configuration);
}

pybind11::list Curve_ProjectOnCurve(Curve& rDummy, const Curve::PointType& rPoint)
{
    double t;
    Curve::PointType Proj;
    int stat = rDummy.ProjectOnCurve(rPoint, Proj, t);
    pybind11::list Output;
    Output.append(t);
    Output.append(Proj);
    Output.append(stat);
    return Output;
}

BRep::PointType Section_ComputeCenter(Section& rDummy)
{
    BRep::PointType Point;
    rDummy.ComputeCenter(Point);
    return Point;
}

pybind11::list Section_Triangulation(Section& rDummy)
{
    std::vector<BRep::PointType> Points;
    std::vector<std::vector<std::size_t> > Connectivities;

    rDummy.Triangulation(Points, Connectivities);

    pybind11::list Output;

    pybind11::list OutputP;
    for (std::size_t i = 0; i < Points.size(); ++i)
        OutputP.append(Points[i]);
    Output.append(OutputP);

    pybind11::list OutputC;
    for (std::size_t i = 0; i < Connectivities.size(); ++i)
    {
        pybind11::list tri;
        for (std::size_t j = 0; j < Connectivities[i].size(); ++j)
            tri.append(Connectivities[i][j]);
        OutputC.append(tri);
    }
    Output.append(OutputC);

    return Output;
}

void BRepApplication_AddBRepAndLevelSetToPython(pybind11::module& m)
{
    /**************************************************************/
    /************** EXPORT INTERFACE FOR SECTION ******************/
    /**************************************************************/

    class_<Section, Section::Pointer>
    ( m, "Section" )
    .def( init<>() )
    .def("ComputeCenter", &Section_ComputeCenter)
    .def("Triangulation", &Section_Triangulation)
    .def("__str__", &PrintObject<Section>)
    ;

    class_<PolygonalSection, PolygonalSection::Pointer, Section>
    ( m, "PolygonalSection" )
    .def( init<>() )
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR BREP **********************/
    /**************************************************************/

    bool(BRep::*pointer_to_IsInside)(const BRep::PointType&) const = &BRep::IsInside;
    int(BRep::*pointer_to_CutStatusElement)(Element::Pointer, const int&) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusGeometry)(Element::GeometryType::Pointer, const int&) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusBySamplingElement)(Element::Pointer, const std::size_t&, const int&) const = &BRep::CutStatusBySampling;
    int(BRep::*pointer_to_CutStatusBySamplingGeometry)(Element::GeometryType::Pointer, const std::size_t&, const int&) const = &BRep::CutStatusBySampling;

    auto pybrep = class_<BRep, BRep::Pointer>
    ( m, "BRep" )
    .def( init<>() )
    .def_property("Name", &BRep::Name, &BRep::SetName)
    .def_property("Tolerance", &BRep::GetTolerance, &BRep::SetTolerance)
    .def("IsInside", &BRep_IsInside2)
    .def("IsInside", &BRep_IsInside3)
    .def("IsInside", pointer_to_IsInside)
    .def("IsInside", &BRep_IsInside4Element)
    .def("CutStatus", pointer_to_CutStatusElement)
    .def("CutStatus", pointer_to_CutStatusGeometry)
    .def("CutStatus", &BRep_CutStatusPoints)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingElement)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingGeometry)
    .def("Bisect", &BRep_Bisect2)
    .def("Bisect", &BRep_Bisect3)
    .def("Intersect", &BRep_IntersectElement)
    .def("ProjectOnSurface", &BRep_ProjectOnSurface)
    .def("Clone", &BRep::CloneBRep)
    .def("__str__", &PrintObject<BRep>)
    ;
    pybrep.attr("_CUT") = pybind11::int_(BRep::_CUT);
    pybrep.attr("_IN") = pybind11::int_(BRep::_IN);
    pybrep.attr("_OUT") = pybind11::int_(BRep::_OUT);

    double(Curve::*pointer_to_ComputeDistance)(const Curve::PointType&) const = &Curve::ComputeDistance;
    Curve::PointType(Curve::*pointer_to_ComputeProjection)(const Curve::PointType&) const = &Curve::ComputeProjection;
    Curve::PointType(Curve::*pointer_to_ComputeIntersection)(const LevelSet&) const = &Curve::ComputeIntersection;

    class_<Curve, Curve::Pointer, FunctionR1R3, DataValueContainer>
    (m, "Curve")
    .def(init<>())
    .def("ComputeDistance", pointer_to_ComputeDistance)
    .def("ComputeProjection", pointer_to_ComputeProjection)
    .def("ComputeIntersection", pointer_to_ComputeIntersection)
    .def("ProjectOnCurve", &Curve_ProjectOnCurve)
    ;

    class_<ParametricCurve, ParametricCurve::Pointer, Curve>
    (m, "ParametricCurve")
    .def(init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    .def("Export", &ParametricCurve::Export)
    ;

    class_<ParametricSurface, ParametricSurface::Pointer, FunctionR2R3>
    (m, "ParametricSurface")
    .def(init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    class_<ParametricVolume, ParametricVolume::Pointer, FunctionR3R3>
    (m, "ParametricVolume")
    .def(init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    double(LevelSet::*LevelSet_pointer_to_GetValue)(const double&, const double&, const double&) const = &LevelSet::GetValue;
    double(LevelSet::*LevelSet_pointer_to_GetValueAtPoint)(const LevelSet::PointType&) const = &LevelSet::GetValue;

    class_<LevelSet, LevelSet::Pointer, FunctionR3R1, BRep>
    ( m, "LevelSet" )
    .def( init<>() )
    .def("GetValue", LevelSet_pointer_to_GetValue)
    .def("GetValue", LevelSet_pointer_to_GetValueAtPoint)
    .def("__str__", &PrintObject<LevelSet>)
    ;

    class_<CircularLevelSet, CircularLevelSet::Pointer, LevelSet>
    ( m, "CircularLevelSet" )
    .def( init<const double&, const double&, const double&>() )
    ;

    class_<Circular2LevelSet, Circular2LevelSet::Pointer, CircularLevelSet>
    ( m, "Circular2LevelSet" )
    .def( init<const double&, const double&, const double&>() )
    ;

    class_<SphericalLevelSet, SphericalLevelSet::Pointer, LevelSet>
    ( m, "SphericalLevelSet" )
    .def( init<const double&, const double&, const double&, const double&>() )
    ;

    class_<Spherical2LevelSet, Spherical2LevelSet::Pointer, SphericalLevelSet>
    ( m, "Spherical2LevelSet" )
    .def( init<const double&, const double&, const double&, const double&>() )
    ;

    class_<DoughnutLevelSet, DoughnutLevelSet::Pointer, LevelSet>
    ( m, "DoughnutLevelSet" )
    .def( init<const double&, const double&>() )
    ;

    class_<CylinderLevelSet, CylinderLevelSet::Pointer, LevelSet>
    ( m, "CylinderLevelSet" )
    .def( init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    .def("CreateQ4Elements", &LevelSet_CreateQ4Elements<CylinderLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoop<CylinderLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoopWithRange<CylinderLevelSet>)
    ;

    class_<Cylinder2LevelSet, Cylinder2LevelSet::Pointer, CylinderLevelSet>
    ( m, "Cylinder2LevelSet" )
    .def( init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    ;

    class_<ConeLevelSet, ConeLevelSet::Pointer, LevelSet>
    ( m, "ConeLevelSet" )
    .def( init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    ;

    class_<LinearLevelSet, LinearLevelSet::Pointer, LevelSet>
    ( m, "LinearLevelSet" )
    .def( init<const double&, const double&, const double&>() )
    .def("CreateLineConditions", &LinearLevelSet_CreateLineConditions)
    ;

    class_<PlanarLevelSet, PlanarLevelSet::Pointer, LevelSet>
    ( m, "PlanarLevelSet" )
    .def( init<const double&, const double&, const double&, const double&>() )
    .def(init<const BRep::PointType&, const BRep::PointType&>())
    ;

    class_<ProductLevelSet, ProductLevelSet::Pointer, LevelSet>
    ( m, "ProductLevelSet" )
    .def( init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    ;

    class_<InverseLevelSet, InverseLevelSet::Pointer, LevelSet>
    ( m, "InverseLevelSet" )
    .def( init<const LevelSet::Pointer>() )
    .def_property("LevelSet", &InverseLevelSet_GetLevelSet, &InverseLevelSet_SetLevelSet)
    ;

    class_<UnionLevelSet, UnionLevelSet::Pointer, LevelSet>
    ( m, "UnionLevelSet" )
    .def( init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    ;

    class_<IntersectionLevelSet, IntersectionLevelSet::Pointer, LevelSet>
    ( m, "IntersectionLevelSet" )
    .def(init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    ;

    class_<DifferenceLevelSet, DifferenceLevelSet::Pointer, LevelSet>
    ( m, "DifferenceLevelSet" )
    .def(init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    ;

    class_<ClosestLevelSet, ClosestLevelSet::Pointer, LevelSet>
    ( m, "ClosestLevelSet" )
    .def(init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    ;

    class_<DistanceToCurveLevelSet, DistanceToCurveLevelSet::Pointer, LevelSet>
    ( m, "DistanceToCurveLevelSet" )
    .def(init<const Curve::Pointer, const double&>() )
    .def("CreateQ4Elements", &LevelSet_CreateQ4Elements<DistanceToCurveLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoop<DistanceToCurveLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoopWithRange<DistanceToCurveLevelSet>)
    .def("CreateQ4ConditionsClosedLoop", &LevelSet_CreateQ4ConditionsClosedLoopWithRange<DistanceToCurveLevelSet>)
    ;

    double(NodalLevelSet::*NodalLevelSet_pointer_to_GetValueAtNode)(const NodalLevelSet::NodeType&) const = &NodalLevelSet::GetValue;

    class_<NodalLevelSet, NodalLevelSet::Pointer, LevelSet>
    ( m, "NodalLevelSet" )
    .def(init<LevelSet::Pointer>() )
    .def_property("Postfix", &NodalLevelSet::Postfix, &NodalLevelSet::SetPostfix)
    .def_property_readonly("IsNodalLevelSet", &NodalLevelSet::IsNodalLevelSet)
    .def_property("OperationMode", &NodalLevelSet::OperationMode, &NodalLevelSet::SetOperationMode)
    .def("Initialize", &NodalLevelSet_InitializeFromNodes)
    .def("Initialize", &NodalLevelSet_InitializeFromElements)
    .def("GetValue", NodalLevelSet_pointer_to_GetValueAtNode)
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR PARTICULAR BREP ***********/
    /**************************************************************/

    class_<AndBRep, AndBRep::Pointer, BRep>
    ( m, "AndBRep" )
    .def(init<BRep::Pointer, BRep::Pointer>() )
    ;

    class_<OrBRep, OrBRep::Pointer, BRep>
    ( m, "OrBRep" )
    .def(init<BRep::Pointer, BRep::Pointer>() )
    ;

    class_<NotBRep, NotBRep::Pointer, BRep>
    ( m, "NotBRep" )
    .def(init<BRep::Pointer>() )
    ;

    #ifdef BREP_APPLICATION_USE_OPENCASCADE
    class_<OCCBRep, OCCBRep::Pointer, BRep>
    ( m, "OCCBRep" )
    .def(init<>() )
    .def("SetShape", &OCCBRep::SetShape)
    ;
    #endif

    class_<NATMArcBRep, NATMArcBRep::Pointer, BRep>
    ( m, "NATMArcBRep" )
    .def(init<>() )
    .def("SetReferenceCenter", &NATMArcBRep::SetReferenceCenter)
    .def("SetReferencePoint", &NATMArcBRep::SetReferencePoint)
    .def("AddArc", &NATMArcBRep::AddArc)
    ;

}
}  // namespace Python.
}  // namespace Kratos.
