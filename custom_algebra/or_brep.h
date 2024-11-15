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
//  Date:            12 Jul 2020
//

#if !defined(KRATOS_OR_BREP_H_INCLUDED )
#define  KRATOS_OR_BREP_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>

// External includes

// Project includes
#include "custom_algebra/brep.h"

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

/// BRep representing by OR operation of two breps
class OrBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OrBRep
    KRATOS_CLASS_POINTER_DEFINITION(OrBRep);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OrBRep(BRep::Pointer pBRep1, BRep::Pointer pBRep2)
        : mpBRep1(pBRep1), mpBRep2(pBRep2), BaseType()
    {}

    /// Copy constructor.
    OrBRep(OrBRep const& rOther)
        : BaseType(rOther)
        , mpBRep1(rOther.mpBRep1->CloneBRep())
        , mpBRep2(rOther.mpBRep2->CloneBRep())
    {}

    /// Destructor.
    virtual ~OrBRep() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new OrBRep(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        if (mpBRep1->WorkingSpaceDimension() != mpBRep2->WorkingSpaceDimension())
        {
            KRATOS_ERROR << "The working space dimension is not compatible";
        }
        return mpBRep1->WorkingSpaceDimension();
    }

    std::size_t LocalSpaceDimension() const final
    {
        if (mpBRep1->LocalSpaceDimension() != mpBRep2->LocalSpaceDimension())
        {
            KRATOS_ERROR << "The local space dimension is not compatible";
        }
        return mpBRep1->LocalSpaceDimension();
    }

    /// Check if a point is inside/outside of the BRep
    bool IsInside(const PointType& P) const final
    {
        if (mpBRep1->IsInside(P))
            return true;
        else
            return mpBRep2->IsInside(P);
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in reference configuration
    /// Since now C++ does not support virtual template function, this function must be separated to 2 functions
    bool IsInside0(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        if (mpBRep1->IsInside0(rGeometry, local_coords))
            return true;
        else
            return mpBRep2->IsInside0(rGeometry, local_coords);
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in current configuration
    bool IsInside1(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        if (mpBRep1->IsInside1(rGeometry, local_coords))
            return true;
        else
            return mpBRep2->IsInside1(rGeometry, local_coords);
    }

    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    int CutStatus(GeometryType& r_geom, const int& configuration) const final
    {
        int stat1 = mpBRep1->CutStatus(r_geom, configuration);
        int stat2 = mpBRep2->CutStatus(r_geom, configuration);
        return this->OrCutStatus(stat1, stat2);
    }

    /// Check if a set of points is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    int CutStatus(const std::vector<PointType>& r_points) const final
    {
        int stat1 = mpBRep1->CutStatus(r_points);
        int stat2 = mpBRep2->CutStatus(r_points);
        return this->OrCutStatus(stat1, stat2);
    }

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
        int stat1 = mpBRep1->CutStatus(r_geom, r_local_points, r_points);
        int stat2 = mpBRep2->CutStatus(r_geom, r_local_points, r_points);
        return this->OrCutStatus(stat1, stat2);
    }

    /// Check if a geometry is cut by the BRep by sampling the geometry
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    int CutStatusBySampling(const GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const final
    {
        int stat1 = mpBRep1->CutStatusBySampling(r_geom, nsampling, configuration);
        int stat2 = mpBRep2->CutStatusBySampling(r_geom, nsampling, configuration);
        return this->OrCutStatus(stat1, stat2);
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
    std::string Info() const final
    {
        std::stringstream ss;
        ss << "OR of (" << mpBRep1->Info() << " and " << mpBRep2->Info() << ")";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        mpBRep1->PrintData(rOStream);
        rOStream << std::endl;
        mpBRep2->PrintData(rOStream);
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

    BRep::Pointer mpBRep1;
    BRep::Pointer mpBRep2;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// OR the 2 cut statuses
    int OrCutStatus(const int& stat1, const int& stat2) const
    {
        if (stat1 == _OUT && stat2 == _OUT)
        {
            return _OUT;
        }
        else
        {
            if (stat1 == _IN || stat2 == _IN)
            {
                return _IN;
            }
            else
            {
                return _CUT;
            }
        }
    }

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
    OrBRep& operator=(OrBRep const& rOther);

    ///@}

}; // Class OrBRep

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OR_BREP_H_INCLUDED  defined
