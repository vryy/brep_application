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
//  Date:            5 Jan 2017
//

#if !defined(KRATOS_DIFFERENCE_LEVEL_SET_H_INCLUDED )
#define  KRATOS_DIFFERENCE_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_algebra/level_set/level_set.h"

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
/** Class for intersection of two level sets, i.e LS1 \ LS2
 * REF: Massing et al, CutFEM: Discretizing geometry and partial differential equations
*/
class DifferenceLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DifferenceLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(DifferenceLevelSet);

    typedef LevelSet BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DifferenceLevelSet(const BaseType::Pointer p_level_set_1, const BaseType::Pointer p_level_set_2)
        : BaseType(), mp_level_set_1(p_level_set_1), mp_level_set_2(p_level_set_2)
    {}

    /// Copy constructor.
    DifferenceLevelSet(DifferenceLevelSet const& rOther)
        : BaseType(rOther)
        , mp_level_set_1(rOther.mp_level_set_1->CloneLevelSet())
        , mp_level_set_2(rOther.mp_level_set_2->CloneLevelSet())
    {}

    /// Destructor.
    ~DifferenceLevelSet() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new DifferenceLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return mp_level_set_1->WorkingSpaceDimension();
    }

    double GetValue(const PointType& P) const final
    {
        return std::max(mp_level_set_1->GetValue(P), -mp_level_set_2->GetValue(P));
    }

    Vector GetGradient(const PointType& P) const final
    {
        if (mp_level_set_1->GetValue(P) > -mp_level_set_2->GetValue(P))
        {
            return mp_level_set_1->GetGradient(P);
        }
        else
        {
            return -mp_level_set_2->GetGradient(P);
        }
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
        return "Difference Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
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

    const BaseType::Pointer mp_level_set_1;
    const BaseType::Pointer mp_level_set_2;

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
    DifferenceLevelSet& operator=(DifferenceLevelSet const& rOther);

    ///@}

}; // Class DifferenceLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DIFFERENCE_LEVEL_SET_H_INCLUDED  defined
