//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         brep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            3 Feb 2018
//



// Project includes
#include "includes/deprecated_variables.h"
#include "custom_utilities/delaunay.h"
#include "custom_utilities/brep_utility.h"
#include "custom_utilities/brep_mesh_utility.h"


namespace Kratos
{

template<int TFrame>
void BRepMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling)
{
    Matrix DeltaPosition;

    if (TFrame == 0)
    {
        DeltaPosition.resize(r_geom.size(), 3, false);
        for ( unsigned int node = 0; node < r_geom.size(); ++node )
            noalias( row( DeltaPosition, node ) ) = r_geom[node].Coordinates() - r_geom[node].GetInitialPosition();
    }

    if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Triangle )
    {
        double xi_min = 0.0, xi_max = 1.0;
        double eta_min = 0.0, eta_max = 1.0;

        double dxi = (xi_max - xi_min) / nsampling;
        double deta = (eta_max - eta_min) / nsampling;

        SamplingPoints.clear();
        CoordinatesArrayType loc;
        PointType P;
        for(std::size_t i = 0; i < nsampling+1; ++i)
        {
            loc[0] = xi_min + i*dxi;
            for(std::size_t j = 0; j < nsampling+1; ++j)
            {
                loc[1] = eta_min + j*deta;
                if ( (loc[0] + loc[1]) < 1.0 + 1.0e-10 )
                {
                    r_geom.GlobalCoordinates(P, loc, DeltaPosition);
                    SamplingPoints.push_back(P);
                }
            }
        }
    }
    else if( r_geom.GetGeometryFamily() == GeometryData::Kratos_Quadrilateral
        || (r_geom.GetGeometryFamily() == GeometryData::Kratos_NURBS && r_geom.GetGeometryType() == GeometryData::Kratos_Bezier2D) )
    {
        double xi_min, xi_max, eta_min, eta_max;

        if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Quadrilateral)
        {
            xi_min = -1.0; xi_max = 1.0;
            eta_min = -1.0; eta_max = 1.0;
        }
        else
        {
            xi_min = 0.0; xi_max = 1.0;
            eta_min = 0.0; eta_max = 1.0;
        }

        double dxi = (xi_max - xi_min) / nsampling;
        double deta = (eta_max - eta_min) / nsampling;

        SamplingPoints.reserve((nsampling+1) * (nsampling+1));
        CoordinatesArrayType loc;
        PointType P;
        for(std::size_t i = 0; i < nsampling+1; ++i)
        {
            loc[0] = xi_min + i*dxi;
            for(std::size_t j = 0; j < nsampling+1; ++j)
            {
                loc[1] = eta_min + j*deta;
                if (TFrame == 0)
                    r_geom.GlobalCoordinates(P, loc, DeltaPosition);
                else if (TFrame == 1)
                    r_geom.GlobalCoordinates(P, loc);
                SamplingPoints.push_back(P);
            }
        }
    }
    else if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Tetrahedra )
    {
        double xi_min = 0.0, xi_max = 1.0;
        double eta_min = 0.0, eta_max = 1.0;
        double zeta_min = 0.0, zeta_max = 1.0;

        double dxi = (xi_max - xi_min) / nsampling;
        double deta = (eta_max - eta_min) / nsampling;
        double dzeta = (zeta_max - zeta_min) / nsampling;

        SamplingPoints.clear();
        CoordinatesArrayType loc;
        PointType P;
        for(std::size_t i = 0; i < nsampling+1; ++i)
        {
            loc[0] = xi_min + i*dxi;
            for(std::size_t j = 0; j < nsampling+1; ++j)
            {
                loc[1] = eta_min + j*deta;
                for(std::size_t k = 0; k < nsampling+1; ++k)
                {
                    loc[2] = zeta_min + k*dzeta;
                    if ( (loc[0] + loc[1] + loc[2]) < 1.0 + 1.0e-10 )
                    {
                        if (TFrame == 0)
                            r_geom.GlobalCoordinates(P, loc, DeltaPosition);
                        else if (TFrame == 1)
                            r_geom.GlobalCoordinates(P, loc);
                        SamplingPoints.push_back(P);
                    }
                }
            }
        }
    }
    else if( r_geom.GetGeometryFamily() == GeometryData::Kratos_Hexahedra
        ||  (r_geom.GetGeometryFamily() == GeometryData::Kratos_NURBS && r_geom.GetGeometryType() == GeometryData::Kratos_Bezier3D) )
    {
        double xi_min, xi_max, eta_min, eta_max, zeta_min, zeta_max;

        if(r_geom.GetGeometryFamily() == GeometryData::Kratos_Hexahedra)
        {
            xi_min = -1.0; xi_max = 1.0;
            eta_min = -1.0; eta_max = 1.0;
            zeta_min = -1.0; zeta_max = 1.0;
        }
        else
        {
            xi_min = 0.0; xi_max = 1.0;
            eta_min = 0.0; eta_max = 1.0;
            zeta_min = 0.0; zeta_max = 1.0;
        }

        double dxi = (xi_max - xi_min) / nsampling;
        double deta = (eta_max - eta_min) / nsampling;
        double dzeta = (zeta_max - zeta_min) / nsampling;

        SamplingPoints.reserve((nsampling+1) * (nsampling+1) * (nsampling+1));
        CoordinatesArrayType loc;
        PointType P;
        for(std::size_t i = 0; i < nsampling+1; ++i)
        {
            loc[0] = xi_min + i*dxi;
            for(std::size_t j = 0; j < nsampling+1; ++j)
            {
                loc[1] = eta_min + j*deta;
                for(std::size_t k = 0; k < nsampling+1; ++k)
                {
                    loc[2] = zeta_min + k*dzeta;
                    if (TFrame == 0)
                        r_geom.GlobalCoordinates(P, loc, DeltaPosition);
                    else if (TFrame == 1)
                        r_geom.GlobalCoordinates(P, loc);
                    SamplingPoints.push_back(P);
                }
            }
        }
    }
    else
    {
        std::stringstream ss;
        ss << "Geometry " << r_geom.GetGeometryType() << " is not supported";
        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
    }
}

void BRepMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
            const PointType& rCenter, const PointType& rNormal,
            const double& radius, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial)
{
    PointType zvec;
    zvec[0] = 0.0;
    zvec[1] = 0.0;
    zvec[2] = 1.0;

    PointType T1 = MathUtils<double>::CrossProduct(rNormal, zvec);
    PointType T2 = MathUtils<double>::CrossProduct(rNormal, T1);

    GenerateSamplingPoints(SamplingPoints, rCenter, T1, T2, radius, nsampling_axial, nsampling_radial);
}

void BRepMeshUtility::GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
            const PointType& rCenter, const PointType& rTangent1, const PointType& rTangent2,
            const double& radius, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial)
{
    const double pi = 4.0*std::atan(1.0);

    PointType T1 = rTangent1 / norm_2(rTangent1);
    PointType T2 = rTangent2 / norm_2(rTangent2);

    SamplingPoints.push_back(rCenter);

    for (std::size_t i = 0; i < nsampling_axial; ++i)
    {
        double r = (i+1)*radius / nsampling_axial;
        for (std::size_t j = 0; j < nsampling_radial; ++j)
        {
            double theta = (j+1)*2.0*pi / nsampling_radial;

            PointType P = rCenter + r*(std::cos(theta)*T1 + std::sin(theta)*T2);
            SamplingPoints.push_back(P);
        }
    }
}

BRepMeshUtility::ElementMeshInfoType BRepMeshUtility::CreateLineElements(ModelPart& r_model_part,
        const std::vector<PointType>& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2: L3 elements;
        const bool& close, // if false: open loop; true: close loop
        Properties::Pointer pProperties)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1 = sampling_points.size() - 1;
    // KRATOS_WATCH(last_node_id)

    std::size_t num_1, num_2;
    if (close)
    {
        num_1 = num_division_1 + 1;
    }
    else
    {
        num_1 = num_division_1;
    }

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_1; ++i)
    {
        NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                sampling_points[i][0], sampling_points[i][1], sampling_points[i][2]);
        // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
        NewNodes.push_back(pNewNode);
    }

    // secondly create elements
    std::size_t last_element_id = BRepUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(2);
    else if (type == 2)
        node.resize(3);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_1; ++i)
    {
        temp_element_nodes.clear();

        if (type == 1)
        {
            node[0] = last_node_id_old + i + 1;
            if (i < num_division_1)
            {
                node[1] = last_node_id_old + i + 2;
            }
            else
            {
                node[1] = last_node_id_old + 1;
            }

            temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
            temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
        }
        else if (type == 2)
        {
            // TODO
            KRATOS_THROW_ERROR(std::logic_error, "type == 2", "is not yet implemented")
        }

        Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
        // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
        pNewElement->Set(ACTIVE, true);
        pNewElement->SetValue(IS_INACTIVE, false);
        NewElements.push_back(pNewElement);
    }

    for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        r_model_part.Elements().push_back(*it);
    }

    r_model_part.Elements().Unique();

    std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}


BRepMeshUtility::ConditionMeshInfoSimpleType BRepMeshUtility::CreateTriangleConditions(ModelPart& r_model_part,
    const std::string& sample_condition_name,
    const int& type, // if 1: generate T3 elements; 2: T6 elements;
    const PointType& rCenter, const PointType& rNormal,
    const double& radius, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
    const int& activation_level,
    Properties::Pointer pProperties)
{
    std::vector<PointType> SamplingPoints;
    GenerateSamplingPoints(SamplingPoints, rCenter, rNormal, radius, nsampling_axial, nsampling_radial);

    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_cond_id = BRepUtility::GetLastConditionId(r_model_part);

    return CreateTriangleEntities<Condition, ModelPart::ConditionsContainerType>(r_model_part, r_model_part.Conditions(), SamplingPoints,
        sample_condition_name, last_node_id, last_cond_id, type, activation_level, pProperties);
}


BRepMeshUtility::ElementMeshInfoType BRepMeshUtility::CreateQuadElements(ModelPart& r_model_part,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    Properties::Pointer pProperties)
{
    int initial_activation_level = 0;
    return CreateQuadElements(r_model_part, sampling_points, sample_element_name, type, close_dir, activation_dir, initial_activation_level, pProperties);
}


BRepMeshUtility::ElementMeshInfoType BRepMeshUtility::CreateQuadElements(ModelPart& r_model_part,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    const int& initial_activation_level,
    Properties::Pointer pProperties)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_element_id = BRepUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    bool reverse = false;
    BRepMeshUtility::ElementMeshInfoType Info = CreateQuadEntities<Element, ModelPart::ElementsContainerType>(r_model_part, r_model_part.Elements(),
        sampling_points, sample_element_name, last_node_id, last_element_id, type, close_dir, activation_dir, initial_activation_level, reverse, pProperties);
    std::cout << std::get<1>(Info).size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;
    return Info;
}


BRepMeshUtility::ConditionMeshInfoType BRepMeshUtility::CreateQuadConditions(ModelPart& r_model_part,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_condition_name,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    const int& initial_activation_level,
    const bool& reverse,
    Properties::Pointer pProperties)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_condition_id = BRepUtility::GetLastConditionId(r_model_part);
    // KRATOS_WATCH(last_condition_id)
    BRepMeshUtility::ConditionMeshInfoType Info = CreateQuadEntities<Condition, ModelPart::ConditionsContainerType>(r_model_part, r_model_part.Conditions(),
        sampling_points, sample_condition_name, last_node_id, last_condition_id, type, close_dir, activation_dir, initial_activation_level, reverse, pProperties);
    std::cout << std::get<1>(Info).size() << " " << sample_condition_name << " conditions are created and added to the model_part" << std::endl;
    return Info;
}


template<class TEntityType, class TEntitiesContainerType>
std::tuple<ModelPart::NodesContainerType, TEntitiesContainerType,
    BRepMeshUtility::BoundaryNodesInfoType, BRepMeshUtility::BoundaryLayerInfoType>
BRepMeshUtility::CreateQuadEntities(ModelPart& r_model_part,
    TEntitiesContainerType& rEntities,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    const int& initial_activation_level,
    const bool& reverse,
    Properties::Pointer pProperties)
{
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1 = sampling_points.size() - 1;
    std::size_t num_division_2 = sampling_points[0].size() - 1;
    // KRATOS_WATCH(last_node_id)

    std::size_t num_1, num_2;
    if (close_dir == 1)
    {
        num_1 = num_division_1 + 1;
        num_2 = num_division_2;
    }
    else if (close_dir == 2)
    {
        num_1 = num_division_1;
        num_2 = num_division_2 + 1;
    }
    else
    {
        num_1 = num_division_1;
        num_2 = num_division_2;
    }

    Variable<int>& ACTIVATION_LEVEL_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("ACTIVATION_LEVEL"));

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_division_1+1; ++i)
    {
        for (std::size_t j = 0; j < num_division_2+1; ++j)
        {
            NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                    sampling_points[i][j][0], sampling_points[i][j][1], sampling_points[i][j][2]);
            // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
            NewNodes.push_back(pNewNode);
        }
    }

    // secondly create elements
    TEntityType const& rCloneElement = KratosComponents<TEntityType>::Get(sample_element_name);
    typename TEntityType::NodesArrayType temp_element_nodes;
    TEntitiesContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(4);
    else if (type == 2)
        node.resize(8);
    else if (type == 3)
        node.resize(9);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

    if (activation_dir == 1) activation_level = -num_division_1;

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_1; ++i)
    {
        if (activation_dir == 2) activation_level = -num_division_2;
        for (std::size_t j = 0; j < num_2; ++j)
        {
            temp_element_nodes.clear();

            if (type == 1)
            {
                node[0] = last_node_id_old + i * (num_division_2 + 1) + j + 1;
                node[2] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 1;
                if (j < num_division_2)
                {
                    node[1] = last_node_id_old + i * (num_division_2 + 1) + j + 2;
                    node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + j + 2;
                }
                else
                {
                    node[1] = last_node_id_old + i * (num_division_2 + 1) + 1;
                    node[3] = last_node_id_old + (i + 1) * (num_division_2 + 1) + 1;
                }
                // std::cout << node[0] << " " << node[1] << " " << node[2] << " " << node[3] << std::endl;

                if (reverse == false)
                {
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[3], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[2], NodeKey).base()));
                }
                else
                {
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[0], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[2], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[3], NodeKey).base()));
                    temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[1], NodeKey).base()));
                }
            }
            else if (type == 2)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, "type == 2", "is not yet implemented")
            }
            else if (type == 3)
            {
                // TODO
                KRATOS_THROW_ERROR(std::logic_error, "type == 3", "is not yet implemented")
            }

            typename TEntityType::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
            // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
            pNewElement->Set(ACTIVE, true);
            pNewElement->SetValue(IS_INACTIVE, false);
            if (activation_dir != 0)
                pNewElement->SetValue(ACTIVATION_LEVEL_var, activation_level + initial_activation_level);
            NewElements.push_back(pNewElement);

            if (activation_dir == 2) ++activation_level;
        }
        if (activation_dir == 1) ++activation_level;
    }

    for (typename TEntitiesContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        rEntities.push_back(*it);
    }

    rEntities.Unique();

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}


BRepMeshUtility::ElementMeshInfoType BRepMeshUtility::CreateHexElements(ModelPart& r_model_part,
    const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
    const std::string& sample_element_name,
    const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
    Properties::Pointer pProperties)
{
    std::size_t last_node_id = BRepUtility::GetLastNodeId(r_model_part);
    std::size_t last_node_id_old = last_node_id;
    std::size_t num_division_1, num_division_2, num_division_3;
    std::size_t num_1, num_2, num_3;

    if (type == 1)
    {
        num_division_1 = sampling_points.size() - 1;
        num_division_2 = sampling_points[0].size() - 1;
        num_division_3 = sampling_points[0][0].size() - 1;
        num_1 = num_division_1;
        num_2 = num_division_2;
        num_3 = num_division_3;
    }
    else if (type == 2)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
    }
    else if (type == 3)
    {
        num_division_1 = (sampling_points.size() - 1) / 2;
        num_division_2 = (sampling_points[0].size() - 1) / 2;
        num_division_3 = (sampling_points[0][0].size() - 1) / 2;
        num_1 = 2 * num_division_1;
        num_2 = 2 * num_division_2;
        num_3 = 2 * num_division_3;
    }
    // KRATOS_WATCH(last_node_id)

    // firstly create nodes and add to model_part
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < num_1+1; ++i)
    {
        for (std::size_t j = 0; j < num_2+1; ++j)
        {
            for (std::size_t k = 0; k < num_3+1; ++k)
            {
                NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                        sampling_points[i][j][k][0], sampling_points[i][j][k][1], sampling_points[i][j][k][2]);
                // std::cout << "node " << last_node_id << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
                NewNodes.push_back(pNewNode);
            }
        }
    }

    // secondly create elements
    std::size_t last_element_id = BRepUtility::GetLastElementId(r_model_part);
    // KRATOS_WATCH(last_element_id)
    Element const& rCloneElement = KratosComponents<Element>::Get(sample_element_name);
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType NewElements;
    const std::string NodeKey("Node");
    std::vector<std::size_t> node;
    int activation_level;

    if (type == 1)
        node.resize(8);
    else if (type == 2)
        node.resize(20);
    else if (type == 3)
        node.resize(27);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

//    KRATOS_WATCH(last_node_id_old)
    KRATOS_WATCH(num_division_1)
    KRATOS_WATCH(num_division_2)
    KRATOS_WATCH(num_division_3)

    BoundaryLayerInfoType boundary_layers;
    BoundaryNodesInfoType boundary_nodes;

    for (std::size_t i = 0; i < num_division_1; ++i)
    {
        for (std::size_t j = 0; j < num_division_2; ++j)
        {
            for (std::size_t k = 0; k < num_division_3; ++k)
            {
                temp_element_nodes.clear();

                if (type == 1)
                {
                    node[0] = last_node_id_old + (i * (num_2 + 1) + j) * (num_3 + 1) + k + 1;
                    node[1] = last_node_id_old + (i * (num_2 + 1) + j + 1) * (num_3 + 1) + k + 1;
                    node[2] = last_node_id_old + ((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k + 1;
                    node[3] = last_node_id_old + ((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k + 1;
                    node[4] = node[0] + 1;
                    node[5] = node[1] + 1;
                    node[6] = node[2] + 1;
                    node[7] = node[3] + 1;

                    for (int n = 0; n < 8; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }
                else if (type == 2)
                {
                    node[0] = last_node_id_old + 2*((i * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[1] = last_node_id_old + 2*((i * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[2] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[3] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[4] = node[0] + 2;
                    node[5] = node[1] + 2;
                    node[6] = node[2] + 2;
                    node[7] = node[3] + 2;

                    node[8] = (node[0] + node[1]) / 2;
                    node[9] = (node[1] + node[2]) / 2;
                    node[10] = (node[2] + node[3]) / 2;
                    node[11] = (node[0] + node[3]) / 2;

                    node[12] = (node[0] + node[4]) / 2;
                    node[13] = (node[1] + node[5]) / 2;
                    node[14] = (node[2] + node[6]) / 2;
                    node[15] = (node[3] + node[7]) / 2;

                    node[16] = (node[4] + node[5]) / 2;
                    node[17] = (node[5] + node[6]) / 2;
                    node[18] = (node[6] + node[7]) / 2;
                    node[19] = (node[4] + node[7]) / 2;

                    for (int n = 0; n < 20; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

//                    std::cout << "element " << i << " " << j << " " << k << ":" << std::endl;
//                    for (int n = 0; n < 27; ++n)
//                        std::cout << " " << node[n];
//                    std::cout << std::endl;

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1], node[11], node[10], node[9], node[8]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7], node[16], node[17], node[18], node[19]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3], node[12], node[19], node[15], node[11]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5], node[9], node[14], node[17], node[13]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4], node[8], node[13], node[16], node[12]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6], node[10], node[15], node[18], node[14]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }
                else if (type == 3)
                {
                    node[0] = last_node_id_old + 2*((i * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[1] = last_node_id_old + 2*((i * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[2] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j + 1) * (num_3 + 1) + k) + 1;
                    node[3] = last_node_id_old + 2*(((i + 1) * (num_2 + 1) + j) * (num_3 + 1) + k) + 1;
                    node[4] = node[0] + 2;
                    node[5] = node[1] + 2;
                    node[6] = node[2] + 2;
                    node[7] = node[3] + 2;

                    node[8] = (node[0] + node[1]) / 2;
                    node[9] = (node[1] + node[2]) / 2;
                    node[10] = (node[2] + node[3]) / 2;
                    node[11] = (node[0] + node[3]) / 2;

                    node[12] = (node[0] + node[4]) / 2;
                    node[13] = (node[1] + node[5]) / 2;
                    node[14] = (node[2] + node[6]) / 2;
                    node[15] = (node[3] + node[7]) / 2;

                    node[16] = (node[4] + node[5]) / 2;
                    node[17] = (node[5] + node[6]) / 2;
                    node[18] = (node[6] + node[7]) / 2;
                    node[19] = (node[4] + node[7]) / 2;

                    node[20] = (node[0] + node[1] + node[2] + node[3]) / 4;
                    node[21] = (node[0] + node[1] + node[4] + node[5]) / 4;
                    node[22] = (node[1] + node[2] + node[5] + node[6]) / 4;
                    node[23] = (node[2] + node[3] + node[6] + node[7]) / 4;
                    node[24] = (node[0] + node[3] + node[4] + node[7]) / 4;
                    node[25] = (node[4] + node[5] + node[6] + node[7]) / 4;

                    node[26] = (node[0] + node[1] + node[2] + node[3] + node[4] + node[5] + node[6] + node[7]) / 8;

                    for (int n = 0; n < 27; ++n)
                    {
//                        std::cout << "node " << n << ": " << node[n] << std::endl;
                        temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), node[n], NodeKey).base()));
                    }

//                    std::cout << "element " << i << " " << j << " " << k << ":" << std::endl;
//                    for (int n = 0; n < 27; ++n)
//                        std::cout << " " << node[n];
//                    std::cout << std::endl;

                    // extract the layer information
                    if (k == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[3], node[2], node[1], node[11], node[10], node[9], node[8], node[20]};
                        boundary_layers["xmin"].push_back(layer_cond);
                        boundary_nodes["xmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (k == num_division_3-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[4], node[5], node[6], node[7], node[16], node[17], node[18], node[19], node[25]};
                        boundary_layers["xmax"].push_back(layer_cond);
                        boundary_nodes["xmax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[4], node[7], node[3], node[12], node[19], node[15], node[11], node[24]};
                        boundary_layers["ymin"].push_back(layer_cond);
                        boundary_nodes["ymin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (j == num_division_2-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[1], node[2], node[6], node[5], node[9], node[14], node[17], node[13], node[22]};
                        boundary_layers["ymax"].push_back(layer_cond);
                        boundary_nodes["ymax"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == 0)
                    {
                        std::vector<std::size_t> layer_cond = {node[0], node[1], node[5], node[4], node[8], node[13], node[16], node[12], node[21]};
                        boundary_layers["zmin"].push_back(layer_cond);
                        boundary_nodes["zmin"].insert(layer_cond.begin(), layer_cond.end());
                    }

                    if (i == num_division_1-1)
                    {
                        std::vector<std::size_t> layer_cond = {node[2], node[3], node[7], node[6], node[10], node[15], node[18], node[14], node[23]};
                        boundary_layers["zmax"].push_back(layer_cond);
                        boundary_nodes["zmax"].insert(layer_cond.begin(), layer_cond.end());
                    }
                }

                Element::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
                // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
                pNewElement->Set(ACTIVE, true);
                pNewElement->SetValue(IS_INACTIVE, false);
                NewElements.push_back(pNewElement);
            }
        }
    }

    for (ModelPart::ElementsContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        r_model_part.Elements().push_back(*it);
    }

    r_model_part.Elements().Unique();

    std::cout << NewElements.size() << " " << sample_element_name << " elements are created and added to the model_part" << std::endl;

    return std::make_tuple(NewNodes, NewElements, boundary_nodes, boundary_layers);
}

template<class TEntityType, class TEntitiesContainerType>
std::tuple<ModelPart::NodesContainerType, TEntitiesContainerType> BRepMeshUtility::CreateTriangleEntities(
    ModelPart& r_model_part,
    TEntitiesContainerType& rEntities,
    const std::vector<PointType>& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate T3 elements; 2: T6 elements;
    const int& activation_level,
    Properties::Pointer pProperties)
{
    if (sampling_points.size() < 3)
        KRATOS_THROW_ERROR(std::logic_error, "The given point list is not sufficient for triangulation", "")

    PointType T1 = sampling_points[1] - sampling_points[0];
    T1 /= norm_2(T1);
    PointType Tmp = sampling_points[2] - sampling_points[0];

    PointType N = MathUtils<double>::CrossProduct(T1, Tmp);
    N /= norm_2(N);

    PointType T2 = MathUtils<double>::CrossProduct(N, T1);
    T2 /= norm_2(T2);

    // KRATOS_WATCH(T1)
    // KRATOS_WATCH(T2)
    // KRATOS_WATCH(N)

    // KRATOS_WATCH(sampling_points.size())

    std::vector<double> XYlist;
    // std::vector<std::vector<unsigned int> > Connectivities;

    double xmin = 1.0e99, xmax = -1.0e99, ymin = 1.0e99, ymax = -1.0e99;
    for (std::size_t i = 0; i < sampling_points.size(); ++i)
    {
        double x = inner_prod(sampling_points[i] - sampling_points[0], T1);
        double y = inner_prod(sampling_points[i] - sampling_points[0], T2);
        XYlist.push_back(x);
        XYlist.push_back(y);
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    // KRATOS_WATCH(XYlist.size())
    // for (int i = 0; i < XYlist.size()/2; ++i)
    // {
    //     std::cout << XYlist[2*i] << "," << XYlist[2*i+1] << std::endl;
    // }

    // KRATOS_WATCH(xmin)
    // KRATOS_WATCH(xmax)
    // KRATOS_WATCH(ymin)
    // KRATOS_WATCH(ymax)
    double dx = fabs(xmax - xmin);
    double dy = fabs(ymax - ymin);

    Delaunay D(xmin-dx, xmax+dx, ymin-dy, ymax+dy);
    for (std::size_t i = 0; i < XYlist.size()/2; ++i)
        D.addPoint(XYlist[2*i], XYlist[2*i+1]);
    // D.Print();

    auto triangles = D.getTriangles();
    // std::cout << "list of triangles" << std::endl;
    // KRATOS_WATCH(triangles.size())
    // for (auto t = triangles.begin(); t != triangles.end(); ++t)
    // {
    //     std::cout << t->pi.id << " " << t->pj.id << " " << t->pk.id << std::endl;
    // }

    Variable<int>& ACTIVATION_LEVEL_var = static_cast<Variable<int>&>(KratosComponents<VariableData>::Get("ACTIVATION_LEVEL"));

    // firstly create nodes and add to model_part
    std::size_t old_last_node_id = last_node_id;
    ModelPart::NodesContainerType NewNodes;
    for (std::size_t i = 0; i < sampling_points.size(); ++i)
    {
        NodeType::Pointer pNewNode = r_model_part.CreateNewNode(++last_node_id,
                sampling_points[i][0], sampling_points[i][1], sampling_points[i][2]);
        // std::cout << "node " << pNewNode->Id() << " is created at " << pNewNode->X0() << " " << pNewNode->Y0() << " " << pNewNode->Z0() << std::endl;
        NewNodes.push_back(pNewNode);
    }

    // secondly create elements
    TEntityType const& rCloneElement = KratosComponents<TEntityType>::Get(sample_element_name);
    typename TEntityType::NodesArrayType temp_element_nodes;
    TEntitiesContainerType NewElements;
    const std::string NodeKey("Node");

    for (auto t = triangles.begin(); t != triangles.end(); ++t)
    {
        temp_element_nodes.clear();

        if (type == 1)
        {
            temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), old_last_node_id + t->pi.id-4, NodeKey).base()));
            temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), old_last_node_id + t->pj.id-4, NodeKey).base()));
            temp_element_nodes.push_back(*(BRepUtility::FindKey(r_model_part.Nodes(), old_last_node_id + t->pk.id-4, NodeKey).base()));
        }
        else if (type == 2)
        {
            // TODO
            KRATOS_THROW_ERROR(std::logic_error, "type == 2", "is not yet implemented")
        }

        typename TEntityType::Pointer pNewElement = rCloneElement.Create(++last_element_id, temp_element_nodes, pProperties);
        // std::cout << "element " << pNewElement->Id() << " is created" << std::endl;
        pNewElement->Set(ACTIVE, true);
        pNewElement->SetValue(IS_INACTIVE, false);
        pNewElement->SetValue(ACTIVATION_LEVEL_var, activation_level);
        NewElements.push_back(pNewElement);
    }

    for (typename TEntitiesContainerType::ptr_iterator it = NewElements.ptr_begin(); it != NewElements.ptr_end(); ++it)
    {
        rEntities.push_back(*it);
    }

    rEntities.Unique();

    return std::make_tuple(NewNodes, NewElements);
}

// template instantiation

template std::tuple<ModelPart::NodesContainerType, ModelPart::ElementsContainerType, BRepMeshUtility::BoundaryNodesInfoType, BRepMeshUtility::BoundaryLayerInfoType>
BRepMeshUtility::CreateQuadEntities<Element, ModelPart::ElementsContainerType>(ModelPart& r_model_part,
    ModelPart::ElementsContainerType& rEntities,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    const int& initial_activation_level,
    const bool& reverse,
    Properties::Pointer pProperties);

template std::tuple<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType, BRepMeshUtility::BoundaryNodesInfoType, BRepMeshUtility::BoundaryLayerInfoType>
BRepMeshUtility::CreateQuadEntities<Condition, ModelPart::ConditionsContainerType>(ModelPart& r_model_part,
    ModelPart::ConditionsContainerType& rEntities,
    const std::vector<std::vector<PointType> >& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
    const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
    const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
    const int& initial_activation_level,
    const bool& reverse,
    Properties::Pointer pProperties);

template
std::tuple<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> BRepMeshUtility::CreateTriangleEntities<Element, ModelPart::ElementsContainerType>(
    ModelPart& r_model_part,
    ModelPart::ElementsContainerType& rEntities,
    const std::vector<PointType>& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate T3 elements; 2: T6 elements;
    const int& activation_level,
    Properties::Pointer pProperties);

template
std::tuple<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> BRepMeshUtility::CreateTriangleEntities<Condition, ModelPart::ConditionsContainerType>(
    ModelPart& r_model_part,
    ModelPart::ConditionsContainerType& rEntities,
    const std::vector<PointType>& sampling_points,
    const std::string& sample_element_name,
    std::size_t& last_node_id,
    std::size_t& last_element_id,
    const int& type, // if 1: generate T3 elements; 2: T6 elements;
    const int& activation_level,
    Properties::Pointer pProperties);

template void BRepMeshUtility::GenerateSamplingPoints<0>(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling);

template void BRepMeshUtility::GenerateSamplingPoints<1>(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling);

}  // namespace Kratos.
