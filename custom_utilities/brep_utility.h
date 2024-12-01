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
//  Date:            29 Mar 2017
//

#if !defined(KRATOS_BREP_UTILITY_H_INCLUDED )
#define  KRATOS_BREP_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "containers/pointer_vector_set.h"
#include "brep_application/custom_algebra/function/function.h"

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
/** class for auxilliary routines
*/
class BRepUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRepUtility
    KRATOS_CLASS_POINTER_DEFINITION(BRepUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    struct Edge2
    {
        Edge2(const std::size_t& i, const std::size_t& j) : mi(i), mj(j) {}
        const std::size_t& V1() const {return mi;}
        const std::size_t& V2() const {return mj;}
        bool operator==(const Edge2& e) const
        {
            return ((mi == e.mi) && (mj == e.mj)) || ((mi == e.mj) && (mj == e.mi)) ;
        }
        bool operator<(const Edge2& e) const
        {
            if (mi == e.mi)
            {
                return mj < e.mj;
            }
            else
            {
                return mi < e.mi;
            }
        }
        Edge2& operator=(const Edge2& e)
        {
            mi = e.mi;
            mj = e.mj;
            return *this;
        }
        std::size_t mi, mj;
        struct HashFunction // is necessary for unordered_set and unordered_map
        {
            std::size_t operator()(const Edge2& edge) const
            {
                std::size_t Hash1 = std::hash<std::size_t>()(edge.V1());
                std::size_t Hash2 = std::hash<std::size_t>()(edge.V2());
                return Hash1 ^ Hash2;
            }
        };
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRepUtility() {}

    /// Destructor.
    virtual ~BRepUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Find the element in the KRATOS container with specific key
    template<class TContainerType, class TKeyType>
    static typename TContainerType::iterator FindKey(TContainerType& ThisContainer, TKeyType ThisKey, std::string ComponentName)
    {
        typename TContainerType::iterator i_result;
        if ((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            std::stringstream buffer;
            buffer << ComponentName << " #" << ThisKey << " is not found.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }

    /// Get the last node id of the model part
    static std::size_t GetLastNodeId(ModelPart& r_model_part)
    {
        std::size_t lastNodeId = 0;
        for (typename ModelPart::NodesContainerType::iterator it = r_model_part.Nodes().begin();
                it != r_model_part.Nodes().end(); ++it)
        {
            if (it->Id() > lastNodeId)
            {
                lastNodeId = it->Id();
            }
        }

        return lastNodeId;
    }

    /// Get the last element id of the model part
    static std::size_t GetLastElementId(ModelPart& r_model_part)
    {
        std::size_t lastElementId = 0;
        for (typename ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
                it != r_model_part.Elements().ptr_end(); ++it)
        {
            if ((*it)->Id() > lastElementId)
            {
                lastElementId = (*it)->Id();
            }
        }

        return lastElementId;
    }

    /// Get the last condition id of the model part
    static std::size_t GetLastConditionId(ModelPart& r_model_part)
    {
        std::size_t lastCondId = 0;
        for (typename ModelPart::ConditionsContainerType::ptr_iterator it = r_model_part.Conditions().ptr_begin();
                it != r_model_part.Conditions().ptr_end(); ++it)
        {
            if ((*it)->Id() > lastCondId)
            {
                lastCondId = (*it)->Id();
            }
        }

        return lastCondId;
    }

    /// Get the last properties id of the model_part
    static std::size_t GetLastPropertiesId(ModelPart& r_model_part)
    {
        std::size_t lastPropId = 0;
        for (typename ModelPart::PropertiesContainerType::ptr_iterator it = r_model_part.rProperties().ptr_begin();
                it != r_model_part.rProperties().ptr_end(); ++it)
        {
            if ((*it)->Id() > lastPropId)
            {
                lastPropId = (*it)->Id();
            }
        }

        return lastPropId;
    }

    template<int TDim, typename TConnectivityType>
    static void SwapConnectivity(const TConnectivityType& input, TConnectivityType& output)
    {
        if (output.size() != input.size())
        {
            output.resize(input.size());
        }

        if (TDim == 2) // swapping for surface element/condition
        {
            if (input.size() == 3) // tri3
            {
                output[0] = input[0];
                output[1] = input[2];
                output[2] = input[1];
            }
            else if (input.size() == 6) // tri6
            {
                output[0] = input[0];
                output[1] = input[2];
                output[2] = input[1];
                output[3] = input[5];
                output[4] = input[4];
                output[5] = input[3];
            }
            else if (input.size() == 4) // quad4
            {
                output[0] = input[0];
                output[1] = input[3];
                output[2] = input[2];
                output[3] = input[1];
            }
            else if (input.size() == 8) // quad8
            {
                output[0] = input[0];
                output[1] = input[3];
                output[2] = input[2];
                output[3] = input[1];
                output[4] = input[7];
                output[5] = input[6];
                output[6] = input[5];
                output[7] = input[4];
            }
            else if (input.size() == 9) // quad9
            {
                output[0] = input[0];
                output[1] = input[3];
                output[2] = input[2];
                output[3] = input[1];
                output[4] = input[7];
                output[5] = input[6];
                output[6] = input[5];
                output[7] = input[4];
                output[8] = input[8];
            }
            else
            {
                KRATOS_ERROR << "The input size " << input.size() << " is invalid";
            }
        }
        else if (TDim == 3) // swapping for volume element/condition
        {
            if (input.size() == 4) // tetrahedra
            {
                output[0] = input[0];
                output[1] = input[2];
                output[2] = input[1];
                output[3] = input[3];
            }
            else
            {
                KRATOS_ERROR << "The input size " << input.size() << " is invalid";
            }
        }
    }

    /// Compute the center of a geometry in the reference configuration
    static PointType ComputeCenter(const GeometryType& rGeometry)
    {
        PointType C;
        noalias(C) = ZeroVector(3);

        for (std::size_t i = 0; i < rGeometry.size(); ++i)
        {
            noalias(C) += rGeometry[i].GetInitialPosition();
        }

        C /= rGeometry.size();
        return C;
    }

    /// Compute the center of a list of conditions
    template<typename TConditionsContainerType>
    static PointType ComputeCenter(const TConditionsContainerType& rConditions)
    {
        // get the boundary edges and nodal coordinates
        // std::unordered_set<Edge2> boundary_edges;
        std::unordered_set<Edge2, Edge2::HashFunction> boundary_edges;
        std::vector<Edge2> local_edges;
        std::unordered_map<std::size_t, PointType> nodes;

        for (typename TConditionsContainerType::const_iterator it = rConditions.begin(); it != rConditions.end(); ++it)
        {
            local_edges.clear();
            GetEdges(local_edges, it->GetGeometry());

            for (std::size_t i = 0; i < local_edges.size(); ++i)
            {
                auto it2 = boundary_edges.find(local_edges[i]);
                if (it2 == boundary_edges.end())
                {
                    boundary_edges.insert(local_edges[i]);
                }
                else
                {
                    boundary_edges.erase(local_edges[i]);
                }
            }

            for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
            {
                nodes[it->GetGeometry()[i].Id()] = it->GetGeometry()[i].GetInitialPosition();
            }
        }

        // get the nodes on boundary edges
        std::unordered_set<std::size_t> boundary_nodes;
        for (auto it = boundary_edges.begin(); it != boundary_edges.end(); ++it)
        {
            boundary_nodes.insert(it->V1());
            boundary_nodes.insert(it->V2());
        }

        // compute the center
        PointType P;
        noalias(P) = ZeroVector(3);
        for (auto it = boundary_nodes.begin(); it != boundary_nodes.end(); ++it)
        {
            noalias(P) += nodes[*it];
        }
        P /= boundary_nodes.size();

        return P;
    }

    /// Get the edges of a finite element
    /// Refer to GiD manual for numbering on edges
    static void GetEdges(std::vector<Edge2>& edges, const GeometryType& rGeometry)
    {
        const GeometryData::KratosGeometryType& Type = rGeometry.GetGeometryType();

        if ( Type == GeometryData::KratosGeometryType::Kratos_Triangle2D3 || Type == GeometryData::KratosGeometryType::Kratos_Triangle3D3
                || Type == GeometryData::KratosGeometryType::Kratos_Triangle2D6 || Type == GeometryData::KratosGeometryType::Kratos_Triangle3D6 )
        {
            edges.push_back(Edge2(rGeometry[0].Id(), rGeometry[1].Id()));
            edges.push_back(Edge2(rGeometry[1].Id(), rGeometry[2].Id()));
            edges.push_back(Edge2(rGeometry[2].Id(), rGeometry[0].Id()));
        }
        else if ( Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4
                  || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8
                  || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9 )
        {
            edges.push_back(Edge2(rGeometry[0].Id(), rGeometry[1].Id()));
            edges.push_back(Edge2(rGeometry[1].Id(), rGeometry[2].Id()));
            edges.push_back(Edge2(rGeometry[2].Id(), rGeometry[3].Id()));
            edges.push_back(Edge2(rGeometry[3].Id(), rGeometry[0].Id()));
        }
        else
        {
            KRATOS_ERROR << "Not a 2D geometry type";
        }
    }

    /// Get the edges of a finite element
    /// Refer to GiD manual for numbering on edges
    static std::vector<std::vector<std::size_t> > GetEdges(const GeometryData::KratosGeometryType& Type)
    {
        std::vector<std::vector<std::size_t> > edges;

        if (Type == GeometryData::KratosGeometryType::Kratos_Triangle2D3 || Type == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
        {
            edges.push_back(std::vector<std::size_t> {0, 1});
            edges.push_back(std::vector<std::size_t> {1, 2});
            edges.push_back(std::vector<std::size_t> {2, 0});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
        {
            edges.push_back(std::vector<std::size_t> {0, 1});
            edges.push_back(std::vector<std::size_t> {1, 2});
            edges.push_back(std::vector<std::size_t> {2, 3});
            edges.push_back(std::vector<std::size_t> {3, 0});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Triangle2D6 || Type == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
        {
            edges.push_back(std::vector<std::size_t> {0, 1, 3});
            edges.push_back(std::vector<std::size_t> {1, 2, 4});
            edges.push_back(std::vector<std::size_t> {2, 0, 5});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8
                 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9 || Type == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
        {
            edges.push_back(std::vector<std::size_t> {0, 1, 4});
            edges.push_back(std::vector<std::size_t> {1, 2, 5});
            edges.push_back(std::vector<std::size_t> {2, 3, 6});
            edges.push_back(std::vector<std::size_t> {3, 0, 7});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
        {
            edges.push_back(std::vector<std::size_t> {0, 1});
            edges.push_back(std::vector<std::size_t> {1, 2});
            edges.push_back(std::vector<std::size_t> {2, 0});
            edges.push_back(std::vector<std::size_t> {0, 3});
            edges.push_back(std::vector<std::size_t> {1, 3});
            edges.push_back(std::vector<std::size_t> {2, 3});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10)
        {
            edges.push_back(std::vector<std::size_t> {0, 1, 4});
            edges.push_back(std::vector<std::size_t> {1, 2, 5});
            edges.push_back(std::vector<std::size_t> {2, 0, 6});
            edges.push_back(std::vector<std::size_t> {0, 3, 7});
            edges.push_back(std::vector<std::size_t> {1, 3, 8});
            edges.push_back(std::vector<std::size_t> {2, 3, 9});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8)
        {
            edges.push_back(std::vector<std::size_t> {0, 1});
            edges.push_back(std::vector<std::size_t> {1, 2});
            edges.push_back(std::vector<std::size_t> {2, 3});
            edges.push_back(std::vector<std::size_t> {3, 0});
            edges.push_back(std::vector<std::size_t> {4, 5});
            edges.push_back(std::vector<std::size_t> {5, 6});
            edges.push_back(std::vector<std::size_t> {6, 7});
            edges.push_back(std::vector<std::size_t> {7, 4});
            edges.push_back(std::vector<std::size_t> {0, 4});
            edges.push_back(std::vector<std::size_t> {1, 5});
            edges.push_back(std::vector<std::size_t> {2, 6});
            edges.push_back(std::vector<std::size_t> {3, 7});
        }
        else if (Type == GeometryData::KratosGeometryType::Kratos_Hexahedra3D20 || Type == GeometryData::KratosGeometryType::Kratos_Hexahedra3D27)
        {
            edges.push_back(std::vector<std::size_t> {0, 1, 8});
            edges.push_back(std::vector<std::size_t> {1, 2, 9});
            edges.push_back(std::vector<std::size_t> {2, 3, 10});
            edges.push_back(std::vector<std::size_t> {3, 0, 11});
            edges.push_back(std::vector<std::size_t> {4, 5, 16});
            edges.push_back(std::vector<std::size_t> {5, 6, 17});
            edges.push_back(std::vector<std::size_t> {6, 7, 18});
            edges.push_back(std::vector<std::size_t> {7, 4, 19});
            edges.push_back(std::vector<std::size_t> {0, 4, 12});
            edges.push_back(std::vector<std::size_t> {1, 5, 13});
            edges.push_back(std::vector<std::size_t> {2, 6, 14});
            edges.push_back(std::vector<std::size_t> {3, 7, 15});
        }

        return edges;
    }

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
    BRepUtility& operator=(BRepUtility const& rOther);

    /// Copy constructor.
    BRepUtility(BRepUtility const& rOther);

    ///@}

}; // Class BRepUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRepUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRepUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BREP_UTILITY_H_INCLUDED  defined
