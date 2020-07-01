// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Aug 2019 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes


// Project includes
#include "custom_python/add_utilities_to_python.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_utilities/occ_utility.h"
#endif
#include "custom_utilities/brep_utility.h"
#include "custom_utilities/brep_mesh_utility.h"
#include "custom_utilities/delaunay.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

std::size_t BRepUtility_GetLastNodeId(BRepUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t BRepUtility_GetLastElementId(BRepUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t BRepUtility_GetLastConditionId(BRepUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t BRepUtility_GetLastPropertiesId(BRepUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

boost::python::list BRepMeshUtility_CreateTriangleConditions(BRepMeshUtility& rDummy,
    ModelPart& r_model_part,
    const std::string& sample_condition_name,
    const int& type, // if 1: generate T3 elements; 2: T6 elements;
    const array_1d<double, 3>& rCenter,
    const array_1d<double, 3>& rNormal,
    const double& radius, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
    const int& activation_level,
    Properties::Pointer pProperties)
{
    Element::GeometryType::PointType::PointType C, N;
    noalias(C) = rCenter;
    noalias(N) = rNormal;
    BRepMeshUtility::ConditionMeshInfoSimpleType Results = rDummy.CreateTriangleConditions(r_model_part,
        sample_condition_name, type, C, N, radius, nsampling_axial, nsampling_radial,
        activation_level, pProperties);

    boost::python::list Output;
    Output.append(std::get<0>(Results));
    Output.append(std::get<1>(Results));
    return Output;
}

void BRepApplication_AddUtilitiesToPython()
{

    #ifdef BREP_APPLICATION_USE_OPENCASCADE
    class_<OCCUtility, OCCUtility::Pointer, boost::noncopyable>
    ( "OCCUtility", init<>() )
    .def("MakeBottle", &OCCUtility::MakeBottle)
    .def("MakeSphere", &OCCUtility::MakeSphere)
    .def("ReadSTEP", &OCCUtility::ReadSTEP)
    .def("WriteSTEP", &OCCUtility::WriteSTEP)
    .def(self_ns::str(self))
    ;
    #endif

    class_<BRepUtility, BRepUtility::Pointer, boost::noncopyable>
    ("BRepUtility", init<>())
    .def("GetLastNodeId", &BRepUtility_GetLastNodeId)
    .def("GetLastElementId", &BRepUtility_GetLastElementId)
    .def("GetLastConditionId", &BRepUtility_GetLastConditionId)
    .def("GetLastPropertiesId", &BRepUtility_GetLastPropertiesId)
    ;

    class_<BRepMeshUtility, BRepMeshUtility::Pointer, boost::noncopyable>
    ("BRepMeshUtility", init<>())
    .def("CreateTriangleConditions", &BRepMeshUtility_CreateTriangleConditions)
    ;

    void(Delaunay::*pointer_to_addPoint)(const double&, const double&) = &Delaunay::addPoint;
    class_<Delaunay, boost::shared_ptr<Delaunay>, boost::noncopyable>
    ("Delaunay", init<const double&, const double&, const double&, const double&>())
    .def("AddPoint", pointer_to_addPoint)
    .def("Print", &Delaunay::Print)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

