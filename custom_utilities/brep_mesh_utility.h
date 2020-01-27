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
class BRepMeshUtility
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

    typedef std::map<std::string, std::set<std::size_t> > BoundaryNodesInfoType;

    typedef std::map<std::string, std::vector<std::vector<std::size_t> > > BoundaryLayerInfoType;

    typedef std::tuple<ModelPart::NodesContainerType, ModelPart::ElementsContainerType,
        BoundaryNodesInfoType, BoundaryLayerInfoType> MeshInfoType;

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

    /// Generate the sampling points on a geometry in the reference configuration
    static void GenerateSamplingPoints0(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling);

    /// Generate the sampling points on a geometry in the current configuration
    static void GenerateSamplingPoints(std::vector<PointType>& SamplingPoints,
            GeometryType& r_geom, const std::size_t& nsampling);


    /// Create the line elements based on given points list
    static MeshInfoType CreateLineElements(ModelPart& r_model_part,
        const std::vector<PointType>& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate L2 elements; 2: L3 elements;
        const bool& close, // if false: open loop; true: close loop
        Properties::Pointer pProperties);


    /// Create the quad elements based on given points list
    static MeshInfoType CreateQuadElements(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        Properties::Pointer pProperties);


    /// Create the quad elements based on given points list
    static MeshInfoType CreateQuadElements(ModelPart& r_model_part,
        const std::vector<std::vector<PointType> >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate Q4 elements; 2: Q8 elements; 3: Q9 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir
        const int& initial_activation_level,
        Properties::Pointer pProperties);


    /// Create the hex elements based on given points list
    static MeshInfoType CreateHexElements(ModelPart& r_model_part,
        const std::vector<std::vector<std::vector<PointType> > >& sampling_points,
        const std::string& sample_element_name,
        const int& type, // if 1: generate H8 elements; 2: H20 elements; 3: H27 elements
        const int& close_dir, // if 0: open loop; 1: close on 1st dir; 2: close on 2nd dir; 3: close on 3rd dir
        const int& activation_dir, // if 0: no activation; 1: activation on 1st dir; 2: activation on 2nd dir; r: activation on 3rd dir
        Properties::Pointer pProperties);

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
        return "BRep Utility";
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
