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

#if !defined(KRATOS_POW_FUNCTION_H_INCLUDED )
#define  KRATOS_POW_FUNCTION_H_INCLUDED

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
///@name  PowFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general PowFunction
*/
template<class TFunction>
class PowFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowFunction
    KRATOS_CLASS_POINTER_DEFINITION(PowFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PowFunction(const double a, const typename BaseType::Pointer p_func)
        : BaseType(), ma(a), mp_func(p_func)
    {}

    PowFunction(const typename BaseType::Pointer p_func, const double a)
        : BaseType(), ma(a), mp_func(p_func)
    {}

    /// Copy constructor.
    PowFunction(PowFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction()), ma(rOther.ma)
    {}

    /// Destructor.
    virtual ~PowFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new PowFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return pow(mp_func->GetValue(P), ma);
    }

    std::string GetFormula(const std::string& Format) const final
    {
        std::stringstream ss;
        if (Format == "matlab")
        {
            if (ma == 1.0)
            {
                ss << mp_func->GetFormula(Format);
            }
            else if (ma == 0.0)
            {
                ss << "1.0";
            }
            else
            {
                ss << "(" << mp_func->GetFormula(Format) << ")^" << ma;
            }
        }
        return ss.str();
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        if (ma == 1.0)
        {
            return mp_func->GetDiffFunction(component);
        }
        else if (ma == 0.0)
        {
            return typename BaseType::Pointer(new ZeroFunction<TFunction>());
        }
        else
            return typename BaseType::Pointer(
                       new ProductFunction<TFunction>(
                           typename BaseType::Pointer(
                               new ScaleFunction<TFunction>(ma,
                                       typename BaseType::Pointer(new PowFunction(ma - 1, mp_func))
                                                           )
                           ),
                           mp_func->GetDiffFunction(component)
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
        return "Pow Function of " + mp_func->Info();
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

    double ma;
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
    PowFunction& operator=(PowFunction const& rOther);

    ///@}

}; // Class PowFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream PowFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, PowFunction<TFunction>& rThis)
{}

/// output stream PowFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const PowFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_POW_FUNCTION_H_INCLUDED  defined
