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

#if !defined(KRATOS_SCALAR_FUNCTION_H_INCLUDED )
#define  KRATOS_SCALAR_FUNCTION_H_INCLUDED

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
///@name  ScalarFunction
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ScalarFunction
*/
template<class TFunction>
class ScalarFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ScalarFunction
    KRATOS_CLASS_POINTER_DEFINITION(ScalarFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScalarFunction(const double& S)
        : BaseType(), mS(S)
    {}

    /// Copy constructor.
    ScalarFunction(ScalarFunction const& rOther) : BaseType(), mS(rOther.mS)
    {}

    /// Destructor.
    virtual ~ScalarFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static typename BaseType::Pointer Create(const double& S)
    {
        return typename BaseType::Pointer(new ScalarFunction(S));
    }

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new ScalarFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return mS;
    }

    std::string GetFormula(const std::string& Format) const final
    {
        std::stringstream ss;
        ss << mS;
        return ss.str();
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(new ZeroFunction<BaseType>());
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
        return "Scalar Function";
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

    const double mS;

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
    ScalarFunction& operator=(ScalarFunction const& rOther);

    ///@}

}; // Class ScalarFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream ScalarFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, ScalarFunction<TFunction>& rThis)
{}

/// output stream ScalarFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const ScalarFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SCALAR_FUNCTION_H_INCLUDED  defined
