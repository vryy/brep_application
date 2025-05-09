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
//  Date:            24 Feb 2017
//

#if !defined(KRATOS_PARAMETRIC_CURVE_H_INCLUDED )
#define  KRATOS_PARAMETRIC_CURVE_H_INCLUDED

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
/** Implementation for a parametric curve in 3D.
 * A parametric curve is mapping: t \in R^1 -> [x(t), y(t), z(t)]
 */
class ParametricCurve : public Curve
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParametricCurve
    KRATOS_CLASS_POINTER_DEFINITION(ParametricCurve);

    typedef FunctionR1R3 SuperType;

    typedef Curve BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParametricCurve(const FunctionR1R1::Pointer p_func_x,
                    const FunctionR1R1::Pointer p_func_y, const FunctionR1R1::Pointer p_func_z)
        : BaseType(), mp_func_x(p_func_x), mp_func_y(p_func_y), mp_func_z(p_func_z)
    {
        DataValueContainer::SetValue(CURVE_LOWER_BOUND, -1.0);
        DataValueContainer::SetValue(CURVE_UPPER_BOUND, 2.0);
        DataValueContainer::SetValue(CURVE_NUMBER_OF_SAMPLING, 10);
    }

    /// Copy constructor.
    ParametricCurve(ParametricCurve const& rOther)
        : BaseType(rOther)
        , mp_func_x(rOther.mp_func_x->CloneFunction())
        , mp_func_y(rOther.mp_func_y->CloneFunction())
        , mp_func_z(rOther.mp_func_z->CloneFunction())
    {}

    /// Destructor.
    ~ParametricCurve() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// inherit from Function
    SuperType::Pointer CloneFunction() const final
    {
        return SuperType::Pointer(new ParametricCurve(*this));
    }

    /// inherit from Curve
    Curve::Pointer Clone() const final
    {
        return BaseType::Pointer(new ParametricCurve(*this));
    }

    /// inherit from Function
    OutputType GetValue(const InputType& t) const final
    {
        OutputType P;

        P[0] = mp_func_x->GetValue(t);
        P[1] = mp_func_y->GetValue(t);
        P[2] = mp_func_z->GetValue(t);

        return P;
    }

    /// inherit from Function
    OutputType GetDerivative(const int& component, const InputType& t) const final
    {
        OutputType P;

        P[0] = mp_func_x->GetDerivative(component, t);
        P[1] = mp_func_y->GetDerivative(component, t);
        P[2] = mp_func_z->GetDerivative(component, t);

        return P;
    }

    /// inherit from Function
    OutputType GetSecondDerivative(const int& component_1, const int& component_2, const InputType& t) const final
    {
        OutputType P;

        P[0] = mp_func_x->GetSecondDerivative(component_1, component_2, t);
        P[1] = mp_func_y->GetSecondDerivative(component_1, component_2, t);
        P[2] = mp_func_z->GetSecondDerivative(component_1, component_2, t);

        return P;
    }

    /// inherit from Function
    SuperType::Pointer GetDiffFunction(const int& component) const final
    {
        return SuperType::Pointer(
                   new ParametricCurve(
                       mp_func_x->GetDiffFunction(component),
                       mp_func_y->GetDiffFunction(component),
                       mp_func_z->GetDiffFunction(component)
                   )
               );
    }

    void Export(const std::string& filename, const double& tmin, const double& tmax, const std::size_t& nsampling, const int& deriv) const
    {
        std::ofstream fid(filename.c_str());
        fid << std::setprecision(6) << std::scientific;

        OutputType P, DP, D2P;

        fid << "t\t\t\t\tx\t\t\t\ty\t\t\t\tz";
        if (deriv > 0) { fid << "\t\t\t\tdx\t\t\t\tdy\t\t\t\tdz"; }
        if (deriv > 1) { fid << "\t\t\t\td2x\t\t\t\td2y\t\t\t\td2z"; }
        fid << "\n";

        for (std::size_t i = 0; i < nsampling; ++i)
        {
            double t = tmin + i * (tmax - tmin) / (nsampling - 1);
            P = this->GetValue(t);
            fid << t << "\t" << P[0] << "\t" << P[1] << "\t" << P[2];

            if (deriv > 0)
            {
                DP = this->GetDerivative(0, t);
                fid << "\t" << DP[0] << "\t" << DP[1] << "\t" << DP[2];
            }

            if (deriv > 1)
            {
                D2P = this->GetSecondDerivative(0, 0, t);
                fid << "\t" << D2P[0] << "\t" << D2P[1] << "\t" << D2P[2];
            }

            fid << std::endl;
        }

        fid.close();
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
        return "Parametric Curve";
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

    FunctionR1R1::Pointer mp_func_x;
    FunctionR1R1::Pointer mp_func_y;
    FunctionR1R1::Pointer mp_func_z;

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
    ParametricCurve& operator=(ParametricCurve const& rOther);

    ///@}

}; // Class ParametricCurve

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PARAMETRIC_CURVE_H_INCLUDED  defined
