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


#if !defined(KRATOS_CYLINDER_2_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CYLINDER_2_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/cylinder_level_set.h"
#include "custom_utilities/brep_mesh_utility.h"

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

///
/**
 * Level set representing the cylinder
 */
class Cylinder2LevelSet : public CylinderLevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Cylinder2LevelSet
    KRATOS_CLASS_POINTER_DEFINITION(Cylinder2LevelSet);

    typedef CylinderLevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Cylinder2LevelSet(const double& cX, const double& cY, const double& cZ, const double& dX, const double& dY, const double& dZ, const double& R)
    : BaseType(cX, cY, cZ, dX, dY, dZ, R)
    {}

    /// Copy constructor.
    Cylinder2LevelSet(Cylinder2LevelSet const& rOther)
    : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~Cylinder2LevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new Cylinder2LevelSet(*this));
    }


    double GetValue(const PointType& P) const final
    {
        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
//        double pX = (P(0) - mcX) * mdX;
//        double pY = (P(1) - mcY) * mdY;
//        double pZ = (P(2) - mcZ) * mdZ;
        return pow(P(0) - pX, 2) + pow(P(1) - pY, 2) + pow(P(2) - pZ, 2) - pow(mR, 2);
    }


    Vector GetGradient(const PointType& P) const final
    {
//        double pX = (P(0) - mcX) * mdX;
//        double pY = (P(1) - mcY) * mdY;
//        double pZ = (P(2) - mcZ) * mdZ;
//        Vector grad(3);
//        grad(0) = 2.0 * (P(0) - pX) * (1.0 - mdX);
//        grad(1) = 2.0 * (P(1) - pY) * (1.0 - mdY);
//        grad(2) = 2.0 * (P(2) - pZ) * (1.0 - mdZ);

        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
        Vector grad(3);
        grad(0) = 2.0 * (P(0) - pX) * (1.0 - mdX*mdX);
        grad(1) = 2.0 * (P(1) - pY) * (1.0 - mdY*mdY);
        grad(2) = 2.0 * (P(2) - pZ) * (1.0 - mdZ*mdZ);
        return grad;
    }


    Matrix GetGradientDerivatives(const PointType& P) const final
    {
        Matrix Jac(3, 3);

        Jac(0, 0) = 2.0 * (1.0 - mdX*mdX) * (1.0 - mdX*mdX);
        Jac(0, 1) = 2.0 * (1.0 - mdX*mdX) * (-mdX*mdY);
        Jac(0, 2) = 2.0 * (1.0 - mdX*mdX) * (-mdX*mdZ);

        Jac(1, 0) = 2.0 * (1.0 - mdY*mdY) * (-mdY*mdX);
        Jac(1, 1) = 2.0 * (1.0 - mdY*mdY) * (1.0 - mdY*mdY);
        Jac(1, 2) = 2.0 * (1.0 - mdY*mdY) * (-mdY*mdZ);

        Jac(2, 0) = 2.0 * (1.0 - mdZ*mdZ) * (-mdZ*mdX);
        Jac(2, 1) = 2.0 * (1.0 - mdZ*mdZ) * (-mdZ*mdY);
        Jac(2, 2) = 2.0 * (1.0 - mdZ*mdZ) * (1.0 - mdZ*mdZ);

        return Jac;
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
        return "Cylinder-2 Level Set";
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
    Cylinder2LevelSet& operator=(Cylinder2LevelSet const& rOther);

    ///@}

}; // Class Cylinder2LevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                Cylinder2LevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const Cylinder2LevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef PI

#endif // KRATOS_CYLINDER_2_LEVEL_SET_H_INCLUDED  defined
