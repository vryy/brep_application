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
//  Date:            22 Feb 2017
//

#if !defined(KRATOS_NEGATE_FUNCTION_H_INCLUDED )
#define  KRATOS_NEGATE_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"

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
///@name  NegateFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general NegateFunction
*/
template<class TFunction>
class NegateFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NegateFunction
    KRATOS_CLASS_POINTER_DEFINITION(NegateFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NegateFunction(const typename BaseType::Pointer p_func)
        : BaseType(), mp_func(p_func)
    {}

    /// Copy constructor.
    NegateFunction(NegateFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~NegateFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new NegateFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return -mp_func->GetValue(P);
    }

    double GetDerivative(const int& component, const InputType& P) const final
    {
        return -mp_func->GetDerivative(component, P);
    }

    double GetSecondDerivative(const int& component_1, const int& component_2, const InputType& P) const final
    {
        return -mp_func->GetSecondDerivative(component_1, component_2, P);
    }

    std::string GetFormula(const std::string& Format) const final
    {
        std::stringstream ss;
        ss << "-" << mp_func->GetFormula(Format);
        return ss.str();
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(new NegateFunction(mp_func->GetDiffFunction(component)));
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
        return "Negate Function of " + mp_func->Info();
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

    const typename BaseType::Pointer mp_func;

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
    NegateFunction& operator=(NegateFunction const& rOther);

    ///@}

}; // Class NegateFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream NegateFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, NegateFunction<TFunction>& rThis)
{}

/// output stream NegateFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const NegateFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEGATE_FUNCTION_H_INCLUDED  defined
