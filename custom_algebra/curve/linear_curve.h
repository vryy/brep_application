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
//  Date:            24 Oct 2024
//

#if !defined(KRATOS_LINEAR_CURVE_H_INCLUDED )
#define  KRATOS_LINEAR_CURVE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <iomanip>

// External includes

// Project includes
#include "custom_algebra/curve/curve.h"

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
/** Implementation for a linear curve in 3D.
 * A linear curve is mapping: t \in R^1 -> A + t*(B-A)
 */
class LinearCurve : public Curve
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearCurve
    KRATOS_CLASS_POINTER_DEFINITION(LinearCurve);

    typedef FunctionR1R3 SuperType;

    typedef Curve BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearCurve(const double xA, const double yA, const double zA,
                const double xB, const double yB, const double zB)
    : BaseType()
    , mxA(xA), myA(yA), mzA(zA)
    , mxB(xB), myB(yB), mzB(zB)
    {
    }

    /// Copy constructor.
    LinearCurve(LinearCurve const& rOther)
    : BaseType(rOther)
    , mxA(rOther.mxA), myA(rOther.myA), mzA(rOther.mzA)
    , mxB(rOther.mxB), myB(rOther.myB), mzB(rOther.mzB)
    {}

    /// Destructor.
    ~LinearCurve() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// inherit from Function
    SuperType::Pointer CloneFunction() const final
    {
        return SuperType::Pointer(new LinearCurve(*this));
    }

    /// inherit from Curve
    Curve::Pointer Clone() const final
    {
        return BaseType::Pointer(new LinearCurve(*this));
    }

    /// inherit from Function
    OutputType GetValue(const InputType& t) const final
    {
        OutputType P;

        P[0] = mxA + t * (mxB - mxA);
        P[1] = myA + t * (myB - myA);
        P[2] = mzA + t * (mzB - mzA);

        return P;
    }

    /// inherit from Function
    OutputType GetDerivative(const int& component, const InputType& t) const final
    {
        OutputType P;

        P[0] = mxB - mxA;
        P[1] = myB - myA;
        P[2] = mzB - mzA;

        return P;
    }

    /// inherit from Function
    OutputType GetSecondDerivative(const int& component_1, const int& component_2, const InputType& t) const final
    {
        OutputType P;

        P[0] = 0.0;
        P[1] = 0.0;
        P[2] = 0.0;

        return P;
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
        return "Linear Curve";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
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

    double mxA, myA, mzA;
    double mxB, myB, mzB;

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
    LinearCurve& operator=(LinearCurve const& rOther);

    ///@}

}; // Class LinearCurve

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LINEAR_CURVE_H_INCLUDED  defined
