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

#if !defined(KRATOS_AND_BREP_2_H_INCLUDED )
#define  KRATOS_AND_BREP_2_H_INCLUDED

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

/// BRep representing by AND operation of two breps
/// This BRep only provides an IsInside function and relies on the base funtions of BRep
/// to perform cut check operations.
class AndBRep2 : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AndBRep2
    KRATOS_CLASS_POINTER_DEFINITION(AndBRep2);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AndBRep2(BRep::Pointer pBRep1, BRep::Pointer pBRep2)
        : BaseType(), mpBRep1(pBRep1), mpBRep2(pBRep2)
    {}

    /// Copy constructor.
    AndBRep2(AndBRep2 const& rOther)
        : BaseType(rOther)
        , mpBRep1(rOther.mpBRep1->CloneBRep())
        , mpBRep2(rOther.mpBRep2->CloneBRep())
    {}

    /// Destructor.
    ~AndBRep2() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new AndBRep2(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        if (mpBRep1->WorkingSpaceDimension() != mpBRep2->WorkingSpaceDimension())
            KRATOS_ERROR << "The working space dimension is not compatible";
        return mpBRep1->WorkingSpaceDimension();
    }

    std::size_t LocalSpaceDimension() const final
    {
        if (mpBRep1->LocalSpaceDimension() != mpBRep2->LocalSpaceDimension())
            KRATOS_ERROR << "The local space dimension is not compatible";
        return mpBRep1->LocalSpaceDimension();
    }

    /// Check if a point is inside/outside of the BRep
    bool IsInside(const PointType& P) const final
    {
        const bool is_inside1 = mpBRep1->IsInside(P);
        if (is_inside1)
            return mpBRep2->IsInside(P);
        else
            return false;
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
        ss << "AND (2) of (" << mpBRep1->Info() << " and " << mpBRep2->Info() << ")";
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
    AndBRep2& operator=(AndBRep2 const& rOther);

    ///@}

}; // Class AndBRep2

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_AND_BREP_2_H_INCLUDED  defined
