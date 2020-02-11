// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Feb 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "containers/array_1d.h"
#include "custom_python/add_brep_and_level_set_to_python.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/and_brep.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_algebra/occ_brep.h"
#endif
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/level_set/circular_level_set.h"
#include "custom_algebra/level_set/doughnut_level_set.h"
#include "custom_algebra/level_set/spherical_level_set.h"
#include "custom_algebra/level_set/cylinder_level_set.h"
#include "custom_algebra/level_set/cone_level_set.h"
#include "custom_algebra/level_set/linear_level_set.h"
#include "custom_algebra/level_set/planar_level_set.h"
#include "custom_algebra/level_set/product_level_set.h"
#include "custom_algebra/level_set/inverse_level_set.h"
#include "custom_algebra/level_set/union_level_set.h"
#include "custom_algebra/level_set/intersection_level_set.h"
#include "custom_algebra/level_set/difference_level_set.h"
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

using namespace boost::python;

LevelSet::Pointer InverseLevelSet_GetLevelSet(InverseLevelSet& rDummy)
{
    return rDummy.pLeveSet();
}

void InverseLevelSet_SetLevelSet(InverseLevelSet& rDummy, LevelSet::Pointer p_level_set)
{
    rDummy.pSetLevelSet(p_level_set);
}

template<class TLevelSel>
boost::python::list LevelSet_CreateQ4Elements(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const std::size_t& nsampling_axial,
    const std::size_t& nsampling_radial,
    const double& start_radial_angle, // in degree
    const double& end_radial_angle) // in degree
{
    const double Pi = 3.1415926535897932384626433;
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4Elements(r_model_part, sample_element_name, pProperties, nsampling_axial, nsampling_radial,
                start_radial_angle/180.0*Pi, end_radial_angle/180*Pi);
    boost::python::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
boost::python::list LevelSet_CreateQ4ElementsClosedLoop(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const std::size_t& nsampling_axial,
    const std::size_t& nsampling_radial)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4ElementsClosedLoop(r_model_part, sample_element_name, pProperties, nsampling_axial, nsampling_radial);
    boost::python::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
boost::python::list LevelSet_CreateQ4ElementsClosedLoopWithRange(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_element_name,
    Properties::Pointer pProperties,
    const std::size_t& nsampling_axial,
    const std::size_t& nsampling_radial,
    const double& tmin,
    const double& tmax)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> Results
        = rDummy.CreateQ4ElementsClosedLoop(r_model_part, sample_element_name, pProperties, nsampling_axial, nsampling_radial, tmin, tmax);
    boost::python::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

template<class TLevelSel>
boost::python::list LevelSet_CreateQ4ConditionsClosedLoopWithRange(
    TLevelSel& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_condition_name,
    Properties::Pointer pProperties,
    const std::size_t& nsampling_axial,
    const std::size_t& nsampling_radial,
    const double& tmin,
    const double& tmax,
    const bool& reverse)
{
    std::pair<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> Results
        = rDummy.CreateQ4ConditionsClosedLoop(r_model_part, sample_condition_name, pProperties, nsampling_axial, nsampling_radial, tmin, tmax, reverse);
    boost::python::list Output;
    Output.append(Results.first);
    Output.append(Results.second);
    return Output;
}

boost::python::list Curve_ProjectOnCurve(Curve& rDummy, const Curve::PointType& rPoint)
{
    double t;
    Curve::PointType Proj;
    rDummy.ProjectOnCurve(rPoint, Proj, t);
    boost::python::list Output;
    Output.append(t);
    Output.append(Proj);
    return Output;
}

void BRepApplication_AddBRepAndLevelSetToPython()
{
    /**************************************************************/
    /************* EXPORT INTERFACE FOR BREP **********************/
    /**************************************************************/

    int(BRep::*pointer_to_CutStatusElement)(Element::Pointer, const int&) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusGeometry)(Element::GeometryType::Pointer, const int&) const = &BRep::CutStatus;
    int(BRep::*pointer_to_CutStatusBySamplingElement)(Element::Pointer, const std::size_t&, const int&) const = &BRep::CutStatusBySampling;
    int(BRep::*pointer_to_CutStatusBySamplingGeometry)(Element::GeometryType::Pointer, const std::size_t&, const int&) const = &BRep::CutStatusBySampling;

    class_<BRep, BRep::Pointer, boost::noncopyable>
    ( "BRep", init<>() )
    .def("SetTolerance", &BRep::SetTolerance)
    .def("GetTolerance", &BRep::GetTolerance)
    .def("IsInside", &BRep::IsInside)
    .def("CutStatus", pointer_to_CutStatusElement)
    .def("CutStatus", pointer_to_CutStatusGeometry)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingElement)
    .def("CutStatusBySampling", pointer_to_CutStatusBySamplingGeometry)
    .def("Clone", &BRep::CloneBRep)
    .def_readonly("_CUT", &BRep::_CUT)
    .def_readonly("_IN", &BRep::_IN)
    .def_readonly("_OUT", &BRep::_OUT)
    ;

    double(Curve::*pointer_to_ComputeDistance)(const Curve::PointType&) const = &Curve::ComputeDistance;
    Curve::PointType(Curve::*pointer_to_ComputeProjection)(const Curve::PointType&) const = &Curve::ComputeProjection;

    class_<Curve, Curve::Pointer, boost::noncopyable, bases<FunctionR1R3, DataValueContainer> >
    ("Curve", init<>())
    .def("ComputeDistance", pointer_to_ComputeDistance)
    .def("ComputeProjection", pointer_to_ComputeProjection)
    .def("ProjectOnCurve", Curve_ProjectOnCurve)
    ;

    class_<ParametricCurve, ParametricCurve::Pointer, boost::noncopyable, bases<Curve> >
    ("ParametricCurve", init<const FunctionR1R1::Pointer, const FunctionR1R1::Pointer, const FunctionR1R1::Pointer>())
    .def("Export", &ParametricCurve::Export)
    ;

    class_<ParametricSurface, ParametricSurface::Pointer, boost::noncopyable, bases<FunctionR2R3> >
    ("ParametricSurface", init<const FunctionR2R1::Pointer, const FunctionR2R1::Pointer, const FunctionR2R1::Pointer>())
    ;

    class_<ParametricVolume, ParametricVolume::Pointer, boost::noncopyable, bases<FunctionR3R3> >
    ("ParametricVolume", init<const FunctionR3R1::Pointer, const FunctionR3R1::Pointer, const FunctionR3R1::Pointer>())
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR LEVEL SET *****************/
    /**************************************************************/

    double(LevelSet::*LevelSet_pointer_to_GetValue)(const LevelSet::PointType&) const = &LevelSet::GetValue;

    class_<LevelSet, LevelSet::Pointer, boost::noncopyable, bases<FunctionR3R1, BRep> >
    ( "LevelSet", init<>() )
    .def("GetValue", LevelSet_pointer_to_GetValue)
    .def(self_ns::str(self))
    ;

    class_<CircularLevelSet, CircularLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "CircularLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<SphericalLevelSet, SphericalLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "SphericalLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<DoughnutLevelSet, DoughnutLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "DoughnutLevelSet", init<const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<CylinderLevelSet, CylinderLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "CylinderLevelSet", init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    .def("CreateQ4Elements", &LevelSet_CreateQ4Elements<CylinderLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoop<CylinderLevelSet>)
    .def(self_ns::str(self))
    ;

    class_<ConeLevelSet, ConeLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "ConeLevelSet", init<const double&, const double&, const double&, const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<LinearLevelSet, LinearLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "LinearLevelSet", init<const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<PlanarLevelSet, PlanarLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "PlanarLevelSet", init<const double&, const double&, const double&, const double&>() )
    .def(self_ns::str(self))
    ;

    class_<ProductLevelSet, ProductLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "ProductLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<InverseLevelSet, InverseLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "InverseLevelSet", init<const LevelSet::Pointer>() )
    .add_property("LevelSet", &InverseLevelSet_GetLevelSet, &InverseLevelSet_SetLevelSet)
    .def(self_ns::str(self))
    ;

    class_<UnionLevelSet, UnionLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "UnionLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<IntersectionLevelSet, IntersectionLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "IntersectionLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<DifferenceLevelSet, DifferenceLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "DifferenceLevelSet", init<const LevelSet::Pointer, const LevelSet::Pointer>() )
    .def(self_ns::str(self))
    ;

    class_<DistanceToCurveLevelSet, DistanceToCurveLevelSet::Pointer, boost::noncopyable, bases<LevelSet> >
    ( "DistanceToCurveLevelSet", init<const Curve::Pointer, const double&>() )
    .def("CreateQ4Elements", &LevelSet_CreateQ4Elements<DistanceToCurveLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoop<DistanceToCurveLevelSet>)
    .def("CreateQ4ElementsClosedLoop", &LevelSet_CreateQ4ElementsClosedLoopWithRange<DistanceToCurveLevelSet>)
    .def("CreateQ4ConditionsClosedLoop", &LevelSet_CreateQ4ConditionsClosedLoopWithRange<DistanceToCurveLevelSet>)
    .def(self_ns::str(self))
    ;

    /**************************************************************/
    /************* EXPORT INTERFACE FOR PARTICULAR BREP ***********/
    /**************************************************************/

    class_<AndBRep, AndBRep::Pointer, bases<BRep>, boost::noncopyable>
    ( "AndBRep", init<BRep::Pointer, BRep::Pointer>() )
    .def(self_ns::str(self))
    ;

    #ifdef BREP_APPLICATION_USE_OPENCASCADE
    class_<OCCBRep, OCCBRep::Pointer, bases<BRep>, boost::noncopyable>
    ( "OCCBRep", init<>() )
    .def("SetShape", &OCCBRep::SetShape)
    .def(self_ns::str(self))
    ;
    #endif

}
}  // namespace Python.
}  // namespace Kratos.
