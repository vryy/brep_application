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

#if !defined(KRATOS_BREP_MESH_UTILITY_H_INCLUDED )
#define  KRATOS_BREP_MESH_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/section/section.h"

namespace Kratos
{
///@addtogroup BRepApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** class for auxiliary mesh routines
*/
class KRATOS_API(BREP_APPLICATION) BRepMeshUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRepMeshUtility
    KRATOS_CLASS_POINTER_DEFINITION(BRepMeshUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::ElementsContainerType ElementsContainerType;

    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    typedef std::map<std::string, std::set<std::size_t> > BoundaryNodesInfoType;

    typedef std::map<std::string, std::vector<std::vector<std::size_t> > > BoundaryLayerInfoType;

    typedef std::tuple<NodesContainerType, ElementsContainerType, BoundaryNodesInfoType, BoundaryLayerInfoType> ElementMeshInfoType;

    typedef std::tuple<NodesContainerType, ConditionsContainerType> ConditionMeshInfoSimpleType;

    typedef std::tuple<NodesContainerType, ConditionsContainerType, BoundaryNodesInfoType, BoundaryLayerInfoType> ConditionMeshInfoType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRepMeshUtility() {}

    /// Destructor.
    virtual ~BRepMeshUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Generate the sampling local points on a geometry in the reference/current configuration
    static void GenerateSamplingLocalPoints(std::vector<CoordinatesArrayType>& SamplingLocalPoints,
                                            const GeometryType& r_geom, const std::size_t nsampling);

    /// Generate the sampling points on a geometry in the reference/current configuration
    template<int TFrame>
    static void GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
                                       const GeometryType& r_geom, const std::size_t nsampling)
    {
        std::vector<CoordinatesArrayType> SamplingLocalPoints;

        GenerateSamplingLocalPoints(SamplingLocalPoints, r_geom, nsampling);

        SamplingPoints.resize(SamplingLocalPoints.size());
        PointType P;
        for (std::size_t i = 0; i < SamplingLocalPoints.size(); ++i)
        {
            if constexpr (TFrame == 0)
            {
                GlobalCoordinates0(r_geom, P, SamplingLocalPoints[i]);
            }
            else if constexpr (TFrame == 1)
            {
                GlobalCoordinates(r_geom, P, SamplingLocalPoints[i]);
            }
            SamplingPoints[i] = P;
        }
    }

    /// Generate the sampling points along a line
    static void GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
                                       const PointType& StartPoint, const PointType& EndPoint, const std::size_t nsampling);

    /// Generate the sampling points on a circle surface
    static void GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
                                       const PointType& rCenter, const PointType& rNormal,
                                       const double radius, const std::size_t nsampling_axial, const std::size_t nsampling_radial);

    /// Generate the sampling points on a circle surface
    static void GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
                                       const PointType& rCenter, const PointType& rTangent1, const PointType& rTangent2,
                                       const double radius, const std::size_t nsampling_axial, const std::size_t nsampling_radial);

    /// Create the line elements based on given points list
    static ElementMeshInfoType CreateLineElements(ModelPart& r_model_part,
            const std::vector<PointType>& sampling_points,
            const std::string& sample_element_name,
            const int type, // if 1: generate L2 elements; 2: L3 elements;
            const bool close, // if false: open loop; true: close loop
            Properties::Pointer pProperties);

    /// Create the line conditions based on given points list
    static ConditionMeshInfoType CreateLineConditions(ModelPart& r_model_part,
            const std::vector<PointType>& sampling_points,
            const std::string& sample_condition_name,
            const int type, // if 1: generate L2 conditions; 2: L3 conditions;
            const bool close, // if false: open loop; true: close loop
            Properties::Pointer pProperties);

    /// Create the triangle elements on the cylinder
    static ConditionMeshInfoSimpleType CreateTriangleConditions(ModelPart& r_model_part,
            const std::string& sample_condition_name,
            const int type, // if 1: generate T3 elements; 2: T6 elements;
            const PointType& rCenter, const PointType& rNormal,
            const double radius, const std::size_t nsampling_axial, const std::size_t nsampling_radial,
            const int activation_level,
            Properties::Pointer pProperties);

    /// Create the conditions out from the section
    static ConditionMeshInfoSimpleType CreateTriangleConditions(ModelPart& r_model_part,
            const Section& rSection, const std::string& sample_condition_name,
            Properties::Pointer pProperties);

    /// Create the quad elements based on given points list
    static ElementMeshInfoType CreateQuadElements(ModelPart& r_model_part,
            const std::vector<std::vector<PointType> >& sampling_points,
            const std::string& sample_element_name,
            const int type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
            const int close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
            const int activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
            Properties::Pointer pProperties);

    /// Create the quad elements based on given points list
    static ElementMeshInfoType CreateQuadElements(ModelPart& r_model_part,
            const std::vector<std::vector<PointType> >& sampling_points,
            const std::string& sample_element_name,
            const int type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
            const int close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
            const int activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
            const int initial_activation_level,
            Properties::Pointer pProperties);

    /// Create the quad conditions based on given points list
    static ConditionMeshInfoType CreateQuadConditions(ModelPart& r_model_part,
            const std::vector<std::vector<PointType> >& sampling_points,
            const std::string& sample_condition_name,
            const int type, // if 1: generate Q4 conditions; 2: Q8 conditions; 3: Q9 conditions
            const int close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
            const int activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
            const int initial_activation_level,
            const bool reverse,
            Properties::Pointer pProperties);

    /// Create a block of hex elements based on given points list
    static ElementMeshInfoType CreateHexElements(ModelPart& r_model_part,
            const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
            const std::string& sample_element_name,
            const int type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
            Properties::Pointer pProperties);

    /// Helper function to compute the global coordinates in the reference frame
    static CoordinatesArrayType& GlobalCoordinates0( const GeometryType& rGeometry, CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates )
    {
        if (rResult.size() != 3)
        {
            rResult.resize(3, false);
        }
        noalias( rResult ) = ZeroVector( 3 );

        Vector N( rGeometry.size() );
        rGeometry.ShapeFunctionsValues( N, LocalCoordinates );

        for ( std::size_t i = 0 ; i < rGeometry.size() ; ++i )
        {
            noalias( rResult ) += N[i] * rGeometry[i].GetInitialPosition();
        }

        return rResult;
    }

    /// Helper function to compute the global coordinates in the current frame
    static CoordinatesArrayType& GlobalCoordinates( const GeometryType& rGeometry, CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates )
    {
        if (rResult.size() != 3)
        {
            rResult.resize(3, false);
        }
        noalias( rResult ) = ZeroVector( 3 );

        Vector N( rGeometry.size() );
        rGeometry.ShapeFunctionsValues( N, LocalCoordinates );

        for ( std::size_t i = 0 ; i < rGeometry.size() ; ++i )
        {
            noalias( rResult ) += N[i] * rGeometry[i];
        }

        return rResult;
    }

    // /// Providing a list of nodes, project those nodes on the surface and create new nodes
    // static void CreateNodesOnSurface(ModelPart& r_model_part, const std::vector<std::size_t>& nodes, const BRep& r_brep,
    //         std::size_t& last_node_id)
    // {
    //     PointType Proj;

    //     for (std::size_t i = 0; i < nodes.size(); ++i)
    //     {
    //         NodeType& rNode = r_model_part.GetNode(nodes[i]);

    //         r_brep.ProjectOnSurface(rNode, Proj);

    //         r_model_part.CreateNewNode(++last_node_id, Proj[0], Proj[1], Proj[2]);
    //     }
    // }

    /// Providing a list of conditions, project those conditions on the surface and create new conditions
    static ConditionsContainerType CreateConditionsOnSurface(ModelPart& r_model_part,
            const ConditionsContainerType& rConditions, const BRep& r_brep,
            std::size_t& last_node_id, std::size_t& last_cond_id,
            const Condition& rCloneCondition, Properties::Pointer pProperties,
            const bool add_to_model_part);

    /// Providing a list of conditions, project those conditions on the surface and create new conditions
    static std::pair<ConditionsContainerType, ElementsContainerType>
    CreateElementsByProjectingOnSurface(ModelPart& r_model_part,
                                        const ConditionsContainerType& rSourceConditions,
                                        const BRep& r_brep,
                                        std::size_t& last_node_id, std::size_t& last_condition_id, std::size_t& last_element_id,
                                        const Condition& rCloneCondition, const Element& rCloneElement,
                                        Properties::Pointer pProperties,
                                        const bool create_condition,
                                        const bool add_to_model_part);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BRep Mesh Utility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Create the line entities based on given points list
    template<class TEntityType, class TEntitiesContainerType>
    static std::tuple<NodesContainerType, TEntitiesContainerType> CreateLineEntities(
        ModelPart& r_model_part,
        TEntitiesContainerType& rEntities,
        const std::vector<PointType>& sampling_points,
        const TEntityType& rCloneElement,
        std::size_t& last_node_id,
        std::size_t& last_element_id,
        const int type, // if 1: generate L2 elements; 2, 3: L3 elements;
        const bool close, // if false: open loop; true: close loop
        Properties::Pointer pProperties);

    /// Create the mesh of triangle elements based on given points list
    /// Internally, a Delaunay triangulation is used to generate the triangle mesh
    template<class TEntityType, class TEntitiesContainerType>
    static std::tuple<NodesContainerType, TEntitiesContainerType> CreateTriangleEntities(
        ModelPart& r_model_part,
        TEntitiesContainerType& rEntities,
        const std::vector<PointType>& sampling_points,
        const TEntityType& rCloneElement,
        std::size_t& last_node_id,
        std::size_t& last_element_id,
        const int type, // if 1: generate T3 elements; 2: T6 elements;
        const int activation_level,
        Properties::Pointer pProperties);

    /// Create the mesh of triangle elements based on a triangulation
    template<class TEntityType, class TEntitiesContainerType>
    static std::tuple<NodesContainerType, TEntitiesContainerType> CreateTriangleEntities(
        ModelPart& r_model_part,
        TEntitiesContainerType& rEntities,
        const std::vector<PointType>& points,
        const std::vector<std::vector<std::size_t> >& connectivities,
        const TEntityType& rCloneElement,
        std::size_t& last_node_id,
        std::size_t& last_element_id,
        const int activation_level,
        Properties::Pointer pProperties);

    /// Create the patch of quad elements based on given points list
    template<class TEntityType, class TEntitiesContainerType>
    static std::tuple<NodesContainerType, TEntitiesContainerType, BoundaryNodesInfoType, BoundaryLayerInfoType> CreateQuadEntities(
        ModelPart& r_model_part,
        TEntitiesContainerType& rEntities,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        std::size_t& last_node_id,
        std::size_t& last_element_id,
        const int type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        const int initial_activation_level,
        const bool reverse,
        Properties::Pointer pProperties);

    /// Providing a list of surface entities, project those entities on the BRep surface and create new surface entities
    template<class TEntityType, class TEntitiesContainerType>
    static void CreateEntitiesOnSurface(TEntitiesContainerType& rOutputEntities,
                                        ModelPart& r_model_part,
                                        const TEntitiesContainerType& rSourceEntities,
                                        const BRep& r_brep,
                                        std::size_t& last_node_id, std::size_t& last_cond_id,
                                        const TEntityType& rCloneCondition, Properties::Pointer pProperties);

    /// Providing a list of surface entities, project those entities on the surface and create new hexahedra entities
    /// It is assumed that the surface entities must have the same number of nodes
    /// If the create_condition flag is set to true, the surface conditions will also be created
    template<class TSurfaceEntityType, class TSurfaceEntitiesContainerType, class TEntityType, class TEntitiesContainerType>
    static void CreateVolumetricEntitiesByProjectingOnSurface(TSurfaceEntitiesContainerType& rOutputSurfaceEntities,
            TEntitiesContainerType& rOutputEntities,
            ModelPart& r_model_part,
            const TSurfaceEntitiesContainerType& rSurfaceEntities,
            const BRep& r_brep,
            std::size_t& last_node_id, std::size_t& last_cond_id, std::size_t& last_elem_id,
            const TSurfaceEntityType& rCloneSurfaceEntity, const TEntityType& rCloneEntity,
            const bool create_condition, Properties::Pointer pProperties);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    BRepMeshUtility& operator=(BRepMeshUtility const& rOther);

    /// Copy constructor.
    BRepMeshUtility(BRepMeshUtility const& rOther);

    ///@}

}; // Class BRepMeshUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRepMeshUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRepMeshUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BREP_MESH_UTILITY_H_INCLUDED  defined
