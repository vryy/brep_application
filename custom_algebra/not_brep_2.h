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
//  Date:            31 Oct 2024
//

#if !defined(KRATOS_NOT_BREP_2_H_INCLUDED )
#define  KRATOS_NOT_BREP_2_H_INCLUDED

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
/// This BRep only provides an IsInside function and relies on the base funtions of BRep
/// to perform cut check operations.
class NotBRep2 : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NotBRep2
    KRATOS_CLASS_POINTER_DEFINITION(NotBRep2);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NotBRep2(BRep::Pointer pBRep)
        : BaseType(), mpBRep(pBRep)
    {}

    /// Copy constructor.
    NotBRep2(NotBRep2 const& rOther)
        : BaseType(rOther)
        , mpBRep(rOther.mpBRep->CloneBRep())
    {}

    /// Destructor.
    ~NotBRep2() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new NotBRep2(*this));
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
        ss << "NOT (2) of (" << mpBRep->Info() << ")";
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
    NotBRep2& operator=(NotBRep2 const& rOther);

    ///@}

}; // Class NotBRep2

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NOT_BREP_2_H_INCLUDED  defined
