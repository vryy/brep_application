// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Aug 2019 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_utilities_to_python.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_utilities/occ_utility.h"
#endif
#include "custom_utilities/brep_utility.h"
#include "custom_utilities/brep_mesh_utility.h"
#include "custom_utilities/brep_intersection_utility.h"
#include "custom_utilities/delaunay.h"
#include "custom_utilities/tube_mesher.h"

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

template<int TDim, typename TConnectivityType>
boost::python::list BRepUtility_Swap(BRepUtility& rDummy, boost::python::list list_entities)
{
    TConnectivityType input, output;

    input.resize(boost::python::len(list_entities));
    for (std::size_t i = 0; i < boost::python::len(list_entities); ++i)
    {
        input[i] = boost::python::extract<typename TConnectivityType::value_type>(list_entities[i]);
    }

    rDummy.SwapConnectivity<TDim, TConnectivityType>(input, output);

    boost::python::list list_output;
    for (std::size_t i = 0; i < boost::python::len(list_entities); ++i)
    {
        list_output.append(output[i]);
    }

    return list_output;
}

Element::GeometryType::PointType::PointType BRepUtility_ComputerCenterElement(BRepUtility& rDummy, Element::Pointer pElement)
{
    return rDummy.ComputeCenter(pElement->GetGeometry());
}

Condition::GeometryType::PointType::PointType BRepUtility_ComputerCenterCondition(BRepUtility& rDummy, Condition::Pointer pCondition)
{
    return rDummy.ComputeCenter(pCondition->GetGeometry());
}

Condition::GeometryType::PointType::PointType BRepUtility_ComputerCenterConditions(BRepUtility& rDummy, ModelPart::ConditionsContainerType& rConditions)
{
    return rDummy.ComputeCenter(rConditions);
}

boost::python::list BRepMeshUtility_CreateTriangleConditions1(BRepMeshUtility& rDummy,
        ModelPart& r_model_part,
        const std::string& sample_condition_name,
        const int type, // if 1: generate T3 elements; 2: T6 elements;
        const array_1d<double, 3>& rCenter,
        const array_1d<double, 3>& rNormal,
        const double radius, const std::size_t nsampling_axial, const std::size_t nsampling_radial,
        const int activation_level,
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

boost::python::list BRepMeshUtility_CreateTriangleConditions2(BRepMeshUtility& rDummy,
        ModelPart& r_model_part,
        const Section& rSection,
        const std::string& sample_condition_name,
        Properties::Pointer pProperties)
{
    BRepMeshUtility::ConditionMeshInfoSimpleType Results = rDummy.CreateTriangleConditions(r_model_part,
            rSection, sample_condition_name, pProperties);

    boost::python::list Output;
    Output.append(std::get<0>(Results));
    Output.append(std::get<1>(Results));
    return Output;
}

ModelPart::ConditionsContainerType BRepMeshUtility_CreateConditionsOnSurface(BRepMeshUtility& rDummy,
        ModelPart& r_model_part,
        const ModelPart::ConditionsContainerType& rConditions, const BRep& r_brep,
        const std::string& sample_condition_name, Properties::Pointer pProperties,
        const bool& add_to_model_part)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_cond_id = BRepUtility::GetLastConditionId(r_model_part);

    Condition const& rCloneCondition = KratosComponents<Condition>::Get(sample_condition_name);

    return BRepMeshUtility::CreateConditionsOnSurface(r_model_part, rConditions, r_brep,
            last_node_id, last_cond_id, rCloneCondition, pProperties, add_to_model_part);
}

boost::python::list BRepMeshUtility_CreateElementsByProjectingOnSurface(BRepMeshUtility& rDummy,
        ModelPart& r_model_part,
        const ModelPart::ConditionsContainerType& rConditions, const BRep& r_brep,
        const std::string& sample_condition_name,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const bool& create_condition,
        const bool& add_to_model_part)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_condition_id = BRepUtility::GetLastConditionId(r_model_part);
    std::size_t last_element_id = BRepUtility::GetLastElementId(r_model_part);

    Condition const& rCloneCondition = KratosComponents<Condition>::Get(sample_condition_name);
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);

    auto Output = BRepMeshUtility::CreateElementsByProjectingOnSurface(r_model_part, rConditions, r_brep,
                  last_node_id, last_condition_id, last_element_id, rCloneCondition, rCloneElement, pProperties,
                  create_condition, add_to_model_part);

    boost::python::list list_output;
    list_output.append(Output.first);
    list_output.append(Output.second);
    return list_output;
}

template<typename TPointType>
boost::python::list TubeMesher_GetPoints(TubeMesher<TPointType>& dummy)
{
    boost::python::list point_list;
    for (std::size_t i = 0; i < dummy.GetPoints().size(); ++i)
    {
        point_list.append(dummy.GetPoints()[i]);
    }
    return point_list;
}

template<typename TPointType>
boost::python::list TubeMesher_GetElements(TubeMesher<TPointType>& dummy)
{
    boost::python::list element_list;
    for (std::size_t i = 0; i < dummy.GetElements().size(); ++i)
    {
        boost::python::list tmp1;
        for (std::size_t j = 0; j < dummy.GetElements()[i].size(); ++j)
        {
            boost::python::list tmp2;
            for (std::size_t k = 0; k < dummy.GetElements()[i][j].size(); ++k)
            {
                boost::python::list tmp3;
                for (std::size_t l = 0; l < dummy.GetElements()[i][j][k].size(); ++l)
                {
                    boost::python::list tmp4;
                    for (std::size_t m = 0; m < dummy.GetElements()[i][j][k][l].size(); ++m)
                    {
                        tmp4.append(dummy.GetElements()[i][j][k][l][m]);
                    }
                    tmp3.append(tmp4);
                }
                tmp2.append(tmp3);
            }
            tmp1.append(tmp2);
        }
        element_list.append(tmp1);
    }
    return element_list;
}

template<typename TPointType>
boost::python::list TubeMesher_GetConditions(TubeMesher<TPointType>& dummy)
{
    boost::python::list condition_list;
    for (std::size_t i = 0; i < dummy.GetConditions().size(); ++i)
    {
        boost::python::list tmp1;
        for (std::size_t j = 0; j < dummy.GetConditions()[i].size(); ++j)
        {
            boost::python::list tmp2;
            for (std::size_t k = 0; k < dummy.GetConditions()[i][j].size(); ++k)
            {
                boost::python::list tmp3;
                for (std::size_t l = 0; l < dummy.GetConditions()[i][j][k].size(); ++l)
                {
                    tmp3.append(dummy.GetConditions()[i][j][k][l]);
                }
                tmp2.append(tmp3);
            }
            tmp1.append(tmp2);
        }
        condition_list.append(tmp1);
    }
    return condition_list;
}

template<typename TPointType>
boost::python::list TubeMesher_GetSlice1(TubeMesher<TPointType>& dummy, const std::size_t slice,
        const std::size_t layer, const std::size_t sub_layer)
{
    std::vector<std::vector<std::size_t> > conditions;
    dummy.GetSlice(conditions, slice, layer, sub_layer);

    boost::python::list condition_list;
    for (std::size_t i = 0; i < conditions.size(); ++i)
    {
        boost::python::list tmp1;
        for (std::size_t j = 0; j < conditions[i].size(); ++j)
        {
            tmp1.append(conditions[i][j]);
        }
        condition_list.append(tmp1);
    }
    return condition_list;
}

template<typename TPointType>
boost::python::list TubeMesher_GetSlice2(TubeMesher<TPointType>& dummy, const std::size_t slice,
        const std::size_t layer)
{
    std::vector<std::vector<std::vector<std::size_t> > > conditions;
    dummy.GetSlice(conditions, slice, layer);

    boost::python::list condition_list;
    for (std::size_t i = 0; i < conditions.size(); ++i)
    {
        boost::python::list tmp1;
        for (std::size_t j = 0; j < conditions[i].size(); ++j)
        {
            boost::python::list tmp2;
            for (std::size_t k = 0; k < conditions[i][j].size(); ++k)
            {
                tmp2.append(conditions[i][j][k]);
            }
            tmp1.append(tmp2);
        }
        condition_list.append(tmp1);
    }
    return condition_list;
}

template<typename TPointType>
boost::python::list TubeMesher_GetRing(TubeMesher<TPointType>& dummy, const std::size_t ring,
                                       const std::size_t layer)
{
    std::vector<std::vector<std::vector<std::size_t> > > elements;
    dummy.GetRing(elements, ring, layer);

    boost::python::list element_list;
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        boost::python::list tmp1;
        for (std::size_t j = 0; j < elements[i].size(); ++j)
        {
            boost::python::list tmp2;
            for (std::size_t k = 0; k < elements[i][j].size(); ++k)
            {
                tmp2.append(elements[i][j][k]);
            }
            tmp1.append(tmp2);
        }
        element_list.append(tmp1);
    }
    return element_list;
}

template<typename TPointType>
struct TubeMesherWrapper
{
    static typename TubeMesher<TPointType>::Pointer initWrapper(const Curve::Pointer pCurve, boost::python::list r_list,
                                           boost::python::list nsampling_layers,
                                           const std::size_t nsampling_axial, const std::size_t nsampling_radial,
                                           const double rotate_angle, const double start_angle, const double end_angle,
                                           const double tmin, const double tmax,
                                           const int type, const std::size_t last_node_id)
    {
        std::vector<double> r_vec;
        for (int i = 0; i < boost::python::len(r_list); ++i)
        {
            r_vec.push_back(boost::python::extract<double>(r_list[i]));
        }

        std::vector<std::size_t> nsampling_layers_vec;
        for (int i = 0; i < boost::python::len(nsampling_layers); ++i)
        {
            nsampling_layers_vec.push_back(static_cast<std::size_t>(boost::python::extract<int>(nsampling_layers[i])));
        }

        return typename TubeMesher<TPointType>::Pointer(new TubeMesher<TPointType>(pCurve, r_vec, nsampling_layers_vec, nsampling_axial, nsampling_radial,
                                   rotate_angle, start_angle, end_angle, tmin, tmax, type, last_node_id));
    }
};

boost::python::list BRepIntersectionUtility_Intersect(BRepIntersectionUtility& rDummy,
        BRep::Pointer pBRep1, BRep::Pointer pBRep2, const std::size_t nsampling)
{
    std::vector<BRep::PointType> Points;
    rDummy.Intersect(Points, *pBRep1, *pBRep2, nsampling);

    boost::python::list Output;
    for (std::size_t i = 0; i < Points.size(); ++i)
    {
        Output.append(Points[i]);
    }

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
    .def("SwapConnectivityOfSurface", &BRepUtility_Swap<2, std::vector<std::size_t> >)
    .def("SwapConnectivityOfVolume", &BRepUtility_Swap<3, std::vector<std::size_t> >)
    .def("ComputeCenter", &BRepUtility_ComputerCenterElement)
    .def("ComputeCenter", &BRepUtility_ComputerCenterCondition)
    .def("ComputeCenter", &BRepUtility_ComputerCenterConditions)
    ;

    class_<BRepMeshUtility, BRepMeshUtility::Pointer, boost::noncopyable>
    ("BRepMeshUtility", init<>())
    .def("CreateTriangleConditions", &BRepMeshUtility_CreateTriangleConditions1)
    .def("CreateTriangleConditions", &BRepMeshUtility_CreateTriangleConditions2)
    .def("CreateConditionsOnSurface", &BRepMeshUtility_CreateConditionsOnSurface)
    .def("CreateElementsByProjectingOnSurface", &BRepMeshUtility_CreateElementsByProjectingOnSurface)
    ;

    typedef typename Element::GeometryType::PointType::PointType PointType;

    class_<TubeMesher<PointType>, typename TubeMesher<PointType>::Pointer, boost::noncopyable>
    ("TubeMesher", no_init)
    .def("__init__", make_constructor(&TubeMesherWrapper<PointType>::initWrapper))
    .def("GetPoints", &TubeMesher_GetPoints<PointType>)
    .def("GetElements", &TubeMesher_GetElements<PointType>)
    .def("GetConditions", &TubeMesher_GetConditions<PointType>)
    .def("GetSlice", &TubeMesher_GetSlice1<PointType>)
    .def("GetSlice", &TubeMesher_GetSlice2<PointType>)
    .def("GetRing", &TubeMesher_GetRing<PointType>)
    .def("NumberOfLayers", &TubeMesher<PointType>::NumberOfLayers)
    .def("NumberOfSubLayers", &TubeMesher<PointType>::NumberOfSubLayers)
    .def("NumberOfRings", &TubeMesher<PointType>::NumberOfRings)
    .def("NumberOfSegments", &TubeMesher<PointType>::NumberOfSegments)
    ;

    void(Delaunay::*pointer_to_addPoint)(const double, const double) = &Delaunay::addPoint;
    class_<Delaunay, boost::shared_ptr<Delaunay>, boost::noncopyable>
    ("Delaunay", init<const double, const double, const double, const double>())
    .def("AddPoint", pointer_to_addPoint)
    .def("Print", &Delaunay::Print)
    ;

    class_<BRepIntersectionUtility, BRepIntersectionUtility::Pointer, boost::noncopyable>
    ("BRepIntersectionUtility", init<>())
    .def("Intersect", &BRepIntersectionUtility_Intersect)
    ;

}
}  // namespace Python.
}  // namespace Kratos.
