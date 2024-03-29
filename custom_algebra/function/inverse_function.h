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
//  Date:            23 Feb 2017
//

#if !defined(KRATOS_INVERSE_FUNCTION_H_INCLUDED )
#define  KRATOS_INVERSE_FUNCTION_H_INCLUDED

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
///@name  InverseFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general InverseFunction
*/
template<class TFunction>
class InverseFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InverseFunction
    KRATOS_CLASS_POINTER_DEFINITION(InverseFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InverseFunction(const typename BaseType::Pointer p_func)
        : BaseType(), mp_func(p_func)
    {}

    /// Copy constructor.
    InverseFunction(InverseFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~InverseFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new InverseFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return 1.0 / mp_func->GetValue(P);
    }

    std::string GetFormula(const std::string& Format) const final
    {
        std::stringstream ss;
        ss << "1.0/" << mp_func->GetFormula(Format);
        return ss.str();
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(
                   new NegateFunction<TFunction>(
                       typename BaseType::Pointer(
                           new ProductFunction<TFunction>(
                               mp_func->GetDiffFunction(component),
                               typename BaseType::Pointer(
                                   new PowFunction<TFunction>(2, typename BaseType::Pointer( new InverseFunction(mp_func) ) )
                               )
                           )
                       )
                   )
               );
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
        return "Inverse Function of " + mp_func->Info();
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
    InverseFunction& operator=(InverseFunction const& rOther);

    ///@}

}; // Class InverseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream InverseFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, InverseFunction<TFunction>& rThis)
{}

/// output stream InverseFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const InverseFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INVERSE_FUNCTION_H_INCLUDED  defined
