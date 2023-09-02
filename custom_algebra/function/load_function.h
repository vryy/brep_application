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
//  Date:            26 Aug 2019
//

#if !defined(KRATOS_LOAD_FUNCTION_H_INCLUDED )
#define  KRATOS_LOAD_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <sstream>
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
///@name  LoadFunctionR3Rns
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Class for the load apply on line/surface
*/
class LoadFunctionR3Rn : public FunctionR3Rn
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LoadFunctionR3Rn
    KRATOS_CLASS_POINTER_DEFINITION(LoadFunctionR3Rn);

    typedef FunctionR3Rn BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LoadFunctionR3Rn()
        : BaseType()
    {}

    /// Copy constructor.
    LoadFunctionR3Rn(LoadFunctionR3Rn const& rOther)
        : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~LoadFunctionR3Rn()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BaseType::Pointer CloneFunction() const final
    {
        return BaseType::Pointer(new LoadFunctionR3Rn(*this));
    }

    void AddComponent(FunctionR3R1::Pointer pComp)
    {
        mpLoadComponents.push_back(pComp);
    }

    OutputType GetValue(const InputType& P) const final
    {
        std::size_t ncomponent = mpLoadComponents.size();
        Vector Load(ncomponent);

        for (std::size_t i = 0; i < ncomponent; ++i)
        {
            Load(i) = mpLoadComponents[i]->GetValue(P);
        }

        return Load;
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
        std::stringstream ss;
        ss << "Load Function R^3->R^" << mpLoadComponents.size();
        return ss.str();
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

    std::vector<FunctionR3R1::Pointer> mpLoadComponents;

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
    LoadFunctionR3Rn& operator=(LoadFunctionR3Rn const& rOther);

    ///@}

}; // Class LoadFunctionR3Rn

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream LoadFunctionR3Rn
inline std::istream& operator >> (std::istream& rIStream, LoadFunctionR3Rn& rThis)
{}

/// output stream LoadFunctionR3Rn
inline std::ostream& operator << (std::ostream& rOStream, const LoadFunctionR3Rn& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOAD_FUNCTION_H_INCLUDED  defined
