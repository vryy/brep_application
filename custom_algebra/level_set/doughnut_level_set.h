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
//  Date:            5 Jan 2018
//

#if !defined(KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED )
#define  KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
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
/** Detail class definition.
 * REF: Massing et al, CutFEM: Discretizing geometry and partial differential equations
*/
class DoughnutLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DoughnutLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(DoughnutLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DoughnutLevelSet(const double& R, const double& r)
        : BaseType(), mR(R), mr(r)
    {
    }

    /// Copy constructor.
    DoughnutLevelSet(DoughnutLevelSet const& rOther)
        : BaseType(rOther), mR(rOther.mR), mr(rOther.mr)
    {}

    /// Destructor.
    ~DoughnutLevelSet() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new DoughnutLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }

    double GetValue(const PointType& P) const final
    {
        return pow(mR - sqrt(pow(P(0), 2) + pow(P(1), 2)), 2) + pow(P(2), 2) - pow(mr, 2);
    }

    Vector GetGradient(const PointType& P) const final
    {
        Vector grad(3);
        double aux = pow(mR - sqrt(pow(P(0), 2) + pow(P(1), 2)), 2);
        grad(0) = 2.0 * aux * ( -P(0) / sqrt(pow(P(0), 2) + pow(P(1), 2)) );
        grad(1) = 2.0 * aux * ( -P(1) / sqrt(pow(P(0), 2) + pow(P(1), 2)) );
        grad(2) = 2.0 * P(2);
        return grad;
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
        return "Doughnut Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << ", R: " << mR << ", r: " << mr;
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

    double mR, mr;

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
    DoughnutLevelSet& operator=(DoughnutLevelSet const& rOther);

    ///@}

}; // Class DoughnutLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DOUGHNUT_LEVEL_SET_H_INCLUDED  defined
