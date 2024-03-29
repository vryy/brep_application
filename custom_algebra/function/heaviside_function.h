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

#if !defined(KRATOS_HEAVISIDE_FUNCTION_H_INCLUDED )
#define  KRATOS_HEAVISIDE_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/zero_function.h"
#include "custom_algebra/brep.h"

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
///@name  HeavisideFunctions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for a general HeavisideFunction
*/
template<class TFunction>
class HeavisideFunction : public TFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HeavisideFunction
    KRATOS_CLASS_POINTER_DEFINITION(HeavisideFunction);

    typedef TFunction BaseType;

    typedef typename BaseType::InputType InputType;

    typedef typename BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HeavisideFunction(const BRep& r_brep)
        : BaseType(), mr_brep(r_brep)
    {}

    HeavisideFunction(const BRep::Pointer p_brep)
        : BaseType(), mr_brep(*p_brep)
    {}

    /// Copy constructor.
    HeavisideFunction(HeavisideFunction const& rOther)
        : BaseType(rOther), mr_brep(rOther.mr_brep)
    {}

    /// Destructor.
    virtual ~HeavisideFunction()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer CloneFunction() const final
    {
        return typename BaseType::Pointer(new HeavisideFunction(*this));
    }

    double GetValue(const InputType& P) const final
    {
        if (mr_brep.IsInside(static_cast<BRep::PointType>(P)))
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }

    std::string GetFormula(const std::string& Format) const final
    {
        return "H(L)";
    }

    typename BaseType::Pointer GetDiffFunction(const int& component) const final
    {
        return typename BaseType::Pointer(new ZeroFunction<TFunction>());
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
        return "Heaviside Function";
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

    const BRep& mr_brep;

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
    HeavisideFunction& operator=(HeavisideFunction const& rOther);

    ///@}

}; // Class HeavisideFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream HeavisideFunction
template<class TFunction>
inline std::istream& operator >> (std::istream& rIStream, HeavisideFunction<TFunction>& rThis)
{}

/// output stream HeavisideFunction
template<class TFunction>
inline std::ostream& operator << (std::ostream& rOStream, const HeavisideFunction<TFunction>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HEAVISIDE_FUNCTION_H_INCLUDED  defined
