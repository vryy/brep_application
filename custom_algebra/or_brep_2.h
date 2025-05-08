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

#if !defined(KRATOS_OR_BREP_2_H_INCLUDED )
#define  KRATOS_OR_BREP_2_H_INCLUDED

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
/// This BRep only provides an IsInside function and relies on the base funtions of BRep
/// to perform cut check operations.
class OrBRep2 : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OrBRep2
    KRATOS_CLASS_POINTER_DEFINITION(OrBRep2);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OrBRep2(BRep::Pointer pBRep1, BRep::Pointer pBRep2)
        : BaseType(), mpBRep1(pBRep1), mpBRep2(pBRep2)
    {}

    /// Copy constructor.
    OrBRep2(OrBRep2 const& rOther)
        : BaseType(rOther)
        , mpBRep1(rOther.mpBRep1->CloneBRep())
        , mpBRep2(rOther.mpBRep2->CloneBRep())
    {}

    /// Destructor.
    ~OrBRep2() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new OrBRep2(*this));
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
        ss << "OR (2) of (" << mpBRep1->Info() << " and " << mpBRep2->Info() << ")";
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
    OrBRep2& operator=(OrBRep2 const& rOther);

    ///@}

}; // Class OrBRep2

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OR_BREP_2_H_INCLUDED  defined
