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
//  Date:            15 Feb 2017
//

#if !defined(KRATOS_SUM_FUNCTION_H_INCLUDED )
#define  KRATOS_SUM_FUNCTION_H_INCLUDED

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
///@name  SumFunction
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general SumFunction
*/
template<class TFunction>
class SumFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SumFunction
    KRATOS_CLASS_POINTER_DEFINITION(SumFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SumFunction(const typename BaseType::Pointer p_func_1, const typename BaseType::Pointer p_func_2)
        : BaseType(), mp_func_1(p_func_1), mp_func_2(p_func_2)
    {}

    /// Copy constructor.
    SumFunction(SumFunction const& rOther)
        : BaseType(rOther)
        , mp_func_1(rOther.mp_func_1->CloneFunction())
        , mp_func_2(rOther.mp_func_2->CloneFunction())
    {}

    /// Destructor.
    virtual ~SumFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new SumFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return mp_func_1->GetValue(P) + mp_func_2->GetValue(P);
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return mp_func_1->GetFormula(Format) + "+" + mp_func_2->GetFormula(Format);
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(
                   new SumFunction(
                       mp_func_1->GetDiffFunction(component),
                       mp_func_2->GetDiffFunction(component)
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
        return "Sum Function of " + mp_func_1->Info() + " and " + mp_func_2->Info();
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

    const typename BaseType::Pointer mp_func_1;
    const typename BaseType::Pointer mp_func_2;

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
    SumFunction& operator=(SumFunction const& rOther);

    ///@}

}; // Class SumFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream SumFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, SumFunction<TFunction>& rThis)
{}

/// output stream SumFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const SumFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SUM_FUNCTION_H_INCLUDED  defined
