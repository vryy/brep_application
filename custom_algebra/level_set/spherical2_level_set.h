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


#if !defined(KRATOS_SPHERICAL_2_LEVEL_SET_H_INCLUDED )
#define  KRATOS_SPHERICAL_2_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/spherical_level_set.h"


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
/** Detail class definition.
*/
class Spherical2LevelSet : public SphericalLevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Spherical2LevelSet
    KRATOS_CLASS_POINTER_DEFINITION(Spherical2LevelSet);

    typedef SphericalLevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Spherical2LevelSet(const double& cX, const double& cY, const double& cZ, const double& R)
    : BaseType(cX, cY, cZ, R)
    {}

    /// Copy constructor.
    Spherical2LevelSet(Spherical2LevelSet const& rOther)
    : BaseType(rOther)
    {}

    /// Destructor.
    virtual ~Spherical2LevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new Spherical2LevelSet(*this));
    }


    double GetValue(const PointType& P) const final
    {
        return pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) + pow(P(2) - mcZ, 2) - pow(mR, 2);
    }


    Vector GetGradient(const PointType& P) const final
    {
        Vector grad(3);
        grad(0) = 2.0 * (P(0) - mcX);
        grad(1) = 2.0 * (P(1) - mcY);
        grad(2) = 2.0 * (P(2) - mcZ);
        return grad;
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
        return "Spherical-2 Level Set";
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
    Spherical2LevelSet& operator=(Spherical2LevelSet const& rOther);

    ///@}

}; // Class Spherical2LevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                Spherical2LevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const Spherical2LevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERICAL_2_LEVEL_SET_H_INCLUDED  defined
