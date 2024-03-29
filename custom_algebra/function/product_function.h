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
//  Date:            14 Feb 2017
//

#if !defined(KRATOS_PRODUCT_FUNCTION_H_INCLUDED )
#define  KRATOS_PRODUCT_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/sum_function.h"
#include "custom_algebra/function/product_function.h"

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
///@name  ProductFunction
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ProductFunction
*/
template<class TFunction>
class ProductFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ProductFunction
    KRATOS_CLASS_POINTER_DEFINITION(ProductFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ProductFunction(const typename BaseType::Pointer p_func_1, const typename BaseType::Pointer p_func_2)
        : BaseType(), mp_func_1(p_func_1), mp_func_2(p_func_2)
    {}

    /// Copy constructor.
    ProductFunction(ProductFunction const& rOther)
        : BaseType(rOther)
        , mp_func_1(rOther.mp_func_1->CloneFunction())
        , mp_func_2(rOther.mp_func_2->CloneFunction())
    {}

    /// Destructor.
    virtual ~ProductFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new ProductFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return mp_func_1->GetValue(P) * mp_func_2->GetValue(P);
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return mp_func_1->GetFormula(Format) + "*" + mp_func_2->GetFormula(Format);
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(
                   new SumFunction<BaseType>(
                       typename BaseType::Pointer(
                           new ProductFunction(
                               mp_func_1->GetDiffFunction(component),
                               mp_func_2
                           )
                       ),
                       typename BaseType::Pointer(
                           new ProductFunction(
                               mp_func_1,
                               mp_func_2->GetDiffFunction(component)
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
        return "Product Function of " + mp_func_1->Info() + " and " + mp_func_2->Info();
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
    ProductFunction& operator=(ProductFunction const& rOther);

    ///@}

}; // Class ProductFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream ProductFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, ProductFunction<TFunction>& rThis)
{}

/// output stream ProductFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const ProductFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PRODUCT_FUNCTION_H_INCLUDED  defined
