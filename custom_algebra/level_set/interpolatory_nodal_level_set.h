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
//  Date:            25 Feb 2021
//

#if !defined(KRATOS_INTERPOLATORY_NODAL_LEVEL_SET_H_INCLUDED )
#define  KRATOS_INTERPOLATORY_NODAL_LEVEL_SET_H_INCLUDED

// System includes
#include <unordered_map>
#include <fstream>

// External includes

// Project includes
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/brep_mesh_utility.h"
#include "brep_application_variables.h"

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

/** As name suggested, this level set approximates the level set value by using interpolation.
 */
class InterpolatoryNodalLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterpolatoryNodalLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(InterpolatoryNodalLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef Element::IntegrationMethod IntegrationMethod;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InterpolatoryNodalLevelSet()
        : BaseType()
    {}

    /// Constructor with level set
    InterpolatoryNodalLevelSet(LevelSet::Pointer pLevelSet)
        : BaseType()
    {}

    /// Copy constructor.
    InterpolatoryNodalLevelSet(InterpolatoryNodalLevelSet const& rOther)
        : BaseType(rOther)
    {}

    /// Destructor.
    ~InterpolatoryNodalLevelSet() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// A utility function to use in Python interface to tell this is a nodal level set
    bool IsNodalLevelSet() const
    {
        return true;
    }

    /// inherit from BRep
    /// Check if a point is inside/outside of the BRep
    inline bool IsInside0(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        double phi = this->CalculateOnPoint(rGeometry, local_coords);
        return (phi < 0.0);
    }

    /// inherit from BRep
    /// Check if a point is inside/outside of the BRep
    bool IsInside1(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        return this->IsInside0(rGeometry, local_coords);
    }

    /// inherit from BRep
    /// Check if a geometry is cut by the BRep by sampling the geometry
    int CutStatusBySampling(const GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const final
    {
        std::vector<CoordinatesArrayType> SamplingLocalPoints;
        BRepMeshUtility::GenerateSamplingLocalPoints(SamplingLocalPoints, r_geom, nsampling);

        const std::vector<PointType> dummy;

        return this->CutStatus(r_geom, SamplingLocalPoints, dummy);
    }

    /// inherit from BRep
    /// Check if a geometry is cut by the level set
    int CutStatus(GeometryType& r_geom, const int& configuration) const final
    {
        std::vector<std::size_t> in_list, out_list, on_list;
        for (std::size_t v = 0; v < r_geom.size(); ++v)
        {
            double phi = this->GetValue(r_geom[v]);
            if (phi < -this->GetTolerance())
            {
                in_list.push_back(v);
            }
            else if (phi > this->GetTolerance())
            {
                out_list.push_back(v);
            }
            else
            {
                on_list.push_back(v);
            }
        }

        int stat;
        if (in_list.size() == 0 && out_list.size() == 0)
        {
            for (std::size_t v = 0; v < r_geom.size(); ++v)
            {
                KRATOS_WATCH(r_geom[v])
            }
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(on_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_ERROR << "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.";
        }
        else
        {
            if (in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if (out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }

    /// inherit from BRep
    /// Check if a set of points is cut by the BRep
    /// The geometry and the local information of the points are also provided.
    /// This function allows the use of BRep defined on grid (nodes)
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    int CutStatus(const GeometryType& r_geom,
                  const std::vector<CoordinatesArrayType>& r_local_points,
                  const std::vector<PointType>& r_points) const final
    {
        std::vector<std::size_t> in_list, out_list, on_list;
        for (std::size_t v = 0; v < r_local_points.size(); ++v)
        {
            double phi = this->CalculateOnPoint(r_geom, r_local_points[v]);
            if (phi < -this->GetTolerance())
            {
                in_list.push_back(v);
            }
            else if (phi > this->GetTolerance())
            {
                out_list.push_back(v);
            }
            else
            {
                on_list.push_back(v);
            }
        }

        int stat;
        if (in_list.size() == 0 && out_list.size() == 0)
        {
            for (std::size_t v = 0; v < r_local_points.size(); ++v)
            {
                KRATOS_WATCH(r_local_points[v])
            }
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(on_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_ERROR << "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.";
        }
        else
        {
            if (in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if (out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }

    ///@}
    ///@name Access
    ///@{

    /// Get level set value at a node
    virtual double GetValue(const NodeType& rNode) const
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get level set value at integration points of element
    void CalculateOnIntegrationPoints(Element::Pointer pElement, std::vector<double>& rValues) const
    {
        return this->CalculateOnIntegrationPoints(pElement->GetGeometry(), pElement->GetIntegrationMethod(), rValues);
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Interpolatory Nodal Level Set";
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

    /// Get level set value at integration points of geometry
    void CalculateOnIntegrationPoints(const GeometryType& rGeometry, const IntegrationMethod& ThisIntegrationMethod,
                                              std::vector<double>& rValues) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

        rValues.resize(integration_points.size());

        for (std::size_t j = 0; j < integration_points.size(); ++j)
        {
            rValues[j] = 0.0;

            for (std::size_t i = 0; i < rGeometry.size(); ++i)
            {
                rValues[j] += Ncontainer(j, i) * this->GetValue(rGeometry[i]);
            }
        }
    }

    /// Get level set value at point in geometry
    double CalculateOnPoint(const GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
    {
        Vector Ncontainer;
        rGeometry.ShapeFunctionsValues( Ncontainer, rLocalPoint );

        double result = 0.0;
        for (std::size_t i = 0; i < rGeometry.size(); ++i)
        {
            result += Ncontainer(i) * this->GetValue(rGeometry[i]);
        }

        return result;
    }

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
    InterpolatoryNodalLevelSet& operator=(InterpolatoryNodalLevelSet const& rOther);

    ///@}

}; // Class InterpolatoryNodalLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERPOLATORY_NODAL_LEVEL_SET_H_INCLUDED  defined
