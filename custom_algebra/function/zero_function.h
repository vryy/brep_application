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

#if !defined(KRATOS_ZERO_FUNCTION_H_INCLUDED )
#define  KRATOS_ZERO_FUNCTION_H_INCLUDED

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
///@name  ZeroFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general ZeroFunction
*/
template<class TFunction>
class ZeroFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ZeroFunction
    KRATOS_CLASS_POINTER_DEFINITION(ZeroFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ZeroFunction() : BaseType()
    {}

    /// Copy constructor.
    ZeroFunction(ZeroFunction const& rOther) : BaseType(rOther)
    {}

    /// Destructor.
    ~ZeroFunction() override
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static typename BaseType::Pointer Create()
    {
        return typename BaseType::Pointer(new ZeroFunction());
    }

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new ZeroFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        return 0.0;
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return "0.0";
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(new ZeroFunction());
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
        return "Zero Function";
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
    ZeroFunction& operator=(ZeroFunction const& rOther);

    ///@}

}; // Class ZeroFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream ZeroFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, ZeroFunction<TFunction>& rThis)
{}

/// output stream ZeroFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const ZeroFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ZERO_FUNCTION_H_INCLUDED  defined
