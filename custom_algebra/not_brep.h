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
//  Date:            9 Sep 2017
//

#if !defined(KRATOS_NOT_BREP_H_INCLUDED )
#define  KRATOS_NOT_BREP_H_INCLUDED

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

/// BRep representing NOT operation
class NotBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NotBRep
    KRATOS_CLASS_POINTER_DEFINITION(NotBRep);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NotBRep(BRep::Pointer pBRep)
        : mpBRep(pBRep), BaseType()
    {}

    /// Copy constructor.
    NotBRep(NotBRep const& rOther)
        : BaseType(rOther)
        , mpBRep(rOther.mpBRep->CloneBRep())
    {}

    /// Destructor.
    virtual ~NotBRep() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new NotBRep(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return mpBRep->WorkingSpaceDimension();
    }

    std::size_t LocalSpaceDimension() const final
    {
        return mpBRep->LocalSpaceDimension();
    }

    /// Check if a point is inside/outside of the BRep
    bool IsInside(const PointType& P) const final
    {
        return !mpBRep->IsInside(P);
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in reference configuration
    /// Since now C++ does not support virtual template function, this function must be separated to 2 functions
    bool IsInside0(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        return !mpBRep->IsInside0(rGeometry, local_coords);
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in current configuration
    bool IsInside1(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        return !mpBRep->IsInside1(rGeometry, local_coords);
    }

    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    int CutStatus(GeometryType& r_geom, const int& configuration) const final
    {
        int stat = mpBRep->CutStatus(r_geom, configuration);
        return this->NotCutStatus(stat);
    }

    /// Check if a set of points is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    int CutStatus(const std::vector<PointType>& r_points) const final
    {
        int stat = mpBRep->CutStatus(r_points);
        return this->NotCutStatus(stat);
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
        int stat = mpBRep->CutStatus(r_geom, r_local_points, r_points);
        return this->NotCutStatus(stat);
    }

    /// Check if a geometry is cut by the BRep by sampling the geometry
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    int CutStatusBySampling(const GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const final
    {
        int stat = mpBRep->CutStatusBySampling(r_geom, nsampling, configuration);
        return this->NotCutStatus(stat);
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
        ss << "NOT operation of " << mpBRep->Info();
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
        mpBRep->PrintData(rOStream);
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

    BRep::Pointer mpBRep;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Reverse the cut status
    int NotCutStatus(const int& stat) const
    {
        if (stat == _OUT)
        {
            return _IN;
        }
        else if (stat == _IN)
        {
            return _OUT;
        }
        else
        {
            return _CUT;
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
    NotBRep& operator=(NotBRep const& rOther);

    ///@}

}; // Class NotBRep

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NotBRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NotBRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NOT_BREP_H_INCLUDED  defined
