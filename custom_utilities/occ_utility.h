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
//  Date:            23 Aug 2019
//

#if !defined(KRATOS_OPENCASCADE_UTILITY_H_INCLUDED )
#define  KRATOS_OPENCASCADE_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include <TopoDS_Shape.hxx>

// Project includes
#include "includes/define.h"
#include "custom_utilities/occ_define.h"

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
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * class for auxilliary routines for OpenCasCade
 */
class OCCUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OCCUtility
    KRATOS_CLASS_POINTER_DEFINITION(OCCUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OCCUtility() {}

    /// Destructor.
    virtual ~OCCUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Make the OCC bottle
    // REF: https://www.opencascade.com/doc/occt-6.7.0/overview/html/tutorial.html#sec5
    OCC::shared_ptr<TopoDS_Shape> MakeBottle(double width, double height, double thickness)
    {
        return OCC::MakeBottle(width, height, thickness);
    }

    /// Make the OCC sphere
    // REF: https://algotopia.com/contents/opencascade/opencascade_basic
    //      https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_prim_a_p_i___make_sphere.html
    OCC::shared_ptr<TopoDS_Shape> MakeSphere(double cx, double cy, double cz, double r)
    {
        return OCC::MakeSphere(cx, cy, cz, r);
    }

    /// Read the STEP file
    // REF: https://www.opencascade.com/doc/occt-7.0.0/overview/html/occt_user_guides__step.html
    OCC::shared_ptr<TopoDS_Shape> ReadSTEP(std::string filename)
    {
        return OCC::ReadSTEP(Standard_CString(filename.c_str()));
    }

    /// Write the STEP file
    // REF: https://www.opencascade.com/doc/occt-7.0.0/overview/html/occt_user_guides__step.html
    void WriteSTEP(OCC::shared_ptr<TopoDS_Shape> pShape, std::string filename)
    {
        return OCC::WriteSTEP(pShape, Standard_CString(filename.c_str()));
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
    virtual std::string Info() const
    {
        return "OpenCasCade Utility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    OCCUtility& operator=(OCCUtility const& rOther);

    /// Copy constructor.
    OCCUtility(OCCUtility const& rOther);

    ///@}

}; // Class OCCUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, OCCUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const OCCUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OPENCASCADE_UTILITY_H_INCLUDED  defined
