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

#if !defined(KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED )
#define  KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/product_function.h"
#include "custom_algebra/function/negate_function.h"
#include "custom_algebra/function/pow_function.h"
#include "custom_algebra/function/scalar_function.h"

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
///@name  SinFunction
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general CosFunction
*/
template<class TFunction>
class CosFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CosFunction
    KRATOS_CLASS_POINTER_DEFINITION(CosFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CosFunction(const typename BaseType::Pointer p_func)
        : BaseType(), mp_func(p_func)
    {}

    /// Copy constructor.
    CosFunction(CosFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~CosFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new CosFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return cos(mp_func->GetValue(P));
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return "cos(" + mp_func->GetFormula(Format) + ")";
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final;

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
        return "Cos Function";
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
    CosFunction& operator=(CosFunction const& rOther);

    ///@}

}; // Class CosFunction

/** Class for a general SinFunction
*/
template<class TFunction>
class SinFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SinFunction
    KRATOS_CLASS_POINTER_DEFINITION(SinFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SinFunction(const typename BaseType::Pointer p_func)
        : BaseType(), mp_func(p_func)
    {}

    /// Copy constructor.
    SinFunction(SinFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~SinFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new SinFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return sin(mp_func->GetValue(P));
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return "sin(" + mp_func->GetFormula(Format) + ")";
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final;

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
        return "Sin Function";
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
    SinFunction& operator=(SinFunction const& rOther);

    ///@}

}; // Class SinFunction

/** Class for a general AcosFunction
*/
template<class TFunction>
class AcosFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AcosFunction
    KRATOS_CLASS_POINTER_DEFINITION(AcosFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AcosFunction(const typename BaseType::Pointer p_func)
        : BaseType(), mp_func(p_func)
    {}

    /// Copy constructor.
    AcosFunction(AcosFunction const& rOther)
        : BaseType(rOther), mp_func(rOther.mp_func->CloneFunction())
    {}

    /// Destructor.
    virtual ~AcosFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new AcosFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return acos(mp_func->GetValue(P));
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return "acos(" + mp_func->GetFormula(Format) + ")";
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final;

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
        return "Acos Function";
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
    AcosFunction& operator=(AcosFunction const& rOther);

    ///@}

}; // Class AcosFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream SinFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, SinFunction<TFunction>& rThis)
{}

/// output stream SinFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const SinFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/// input stream CosFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, CosFunction<TFunction>& rThis)
{}

/// output stream CosFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const CosFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

/// input stream AcosFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, AcosFunction<TFunction>& rThis)
{}

/// output stream CosFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const AcosFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

template<class TFunction>
typename CosFunction<TFunction>::BaseType::Pointer CosFunction<TFunction>::GetDiffFunction(const int& component) const
{
    return typename BaseType::Pointer(
               new NegateFunction<TFunction>(
                   typename BaseType::Pointer(
                       new ProductFunction<TFunction>(
                           typename BaseType::Pointer(new SinFunction<TFunction>(mp_func)),
                           mp_func->GetDiffFunction(component)
                       )
                   )
               )
           );
}

template<class TFunction>
typename SinFunction<TFunction>::BaseType::Pointer SinFunction<TFunction>::GetDiffFunction(const int& component) const
{
    return typename BaseType::Pointer(
               new ProductFunction<TFunction>(
                   typename BaseType::Pointer(new CosFunction<TFunction>(mp_func)),
                   mp_func->GetDiffFunction(component)
               )
           );
}

template<class TFunction>
typename AcosFunction<TFunction>::BaseType::Pointer AcosFunction<TFunction>::GetDiffFunction(const int& component) const
{
    return typename BaseType::Pointer(
               new NegateFunction<TFunction>(
                   typename BaseType::Pointer(
                       new ProductFunction<TFunction>(
                           mp_func->GetDiffFunction(component),
                           typename BaseType::Pointer(
                               new PowFunction<TFunction>(
                                   -0.5,
                                   typename BaseType::Pointer(
                                       new SumFunction<TFunction>(
                                           typename BaseType::Pointer( new ScalarFunction<TFunction>(1.0) ),
                                           typename BaseType::Pointer(
                                                   new NegateFunction<TFunction>(
                                                           typename BaseType::Pointer( new PowFunction<TFunction>(2.0, mp_func) )
                                                   )
                                           )
                                       )
                                   )
                               )
                           )
                       )
                   )
               )
           );
}

}  // namespace Kratos.

#endif // KRATOS_TRIGONOMETRIC_FUNCTION_H_INCLUDED  defined
