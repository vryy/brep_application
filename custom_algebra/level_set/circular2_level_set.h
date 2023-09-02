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
//  Date:            10 Feb 2017
//

#if !defined(KRATOS_CIRCULAR_2_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CIRCULAR_2_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/circular_level_set.h"

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

///
/**
 * Level set representing a circle
 * f = (X - cx)^2 + (Y - cy)^2 - R^2
 */
class Circular2LevelSet : public CircularLevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Circular2LevelSet
    KRATOS_CLASS_POINTER_DEFINITION(Circular2LevelSet);

    typedef CircularLevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Circular2LevelSet(const double& cX, const double& cY, const double& R)
        : BaseType(cX, cY, R)
    {}

    /// Copy constructor.
    Circular2LevelSet(Circular2LevelSet const& rOther)
        : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~Circular2LevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new Circular2LevelSet(*this));
    }

    double GetValue(const PointType& P) const final
    {
        return pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) - pow(mR, 2);
    }

    Vector GetGradient(const PointType& P) const final
    {
        Vector grad(3);
        grad(0) = 2.0 * (P(0) - mcX);
        grad(1) = 2.0 * (P(1) - mcY);
        grad(2) = 0.0;
        return grad;
    }

    Matrix GetGradientDerivatives(const PointType& P) const final
    {
        Matrix Jac(3, 3);
        noalias(Jac) = ZeroMatrix(3, 3);

        Jac(0, 0) = 2.0;
        Jac(0, 1) = 0.0;

        Jac(1, 0) = 0.0;
        Jac(1, 1) = 2.0;

        return Jac;
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
        return "Circular-2 Level Set";
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
    Circular2LevelSet& operator=(Circular2LevelSet const& rOther);

    ///@}

}; // Class Circular2LevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, Circular2LevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const Circular2LevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CIRCULAR_2_LEVEL_SET_H_INCLUDED  defined
