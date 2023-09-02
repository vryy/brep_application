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

#if !defined(KRATOS_SPHERICAL_LEVEL_SET_H_INCLUDED )
#define  KRATOS_SPHERICAL_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"

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
 * Level set representing the signed-distance function to a sphere
 */
class SphericalLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SphericalLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(SphericalLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SphericalLevelSet(const double& cX, const double& cY, const double& cZ, const double& R)
        : BaseType(), mcX(cX), mcY(cY), mcZ(cZ), mR(R)
    {}

    /// Copy constructor.
    SphericalLevelSet(SphericalLevelSet const& rOther)
        : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mcZ(rOther.mcZ), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~SphericalLevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const override
    {
        return LevelSet::Pointer(new SphericalLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }

    double GetValue(const PointType& P) const override
    {
        return sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) + pow(P(2) - mcZ, 2)) - mR;
    }

    Vector GetGradient(const PointType& P) const override
    {
        Vector grad(3);
        double aux = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) + pow(P(2) - mcZ, 2));
        grad(0) = (P(0) - mcX) / aux;
        grad(1) = (P(1) - mcY) / aux;
        grad(2) = (P(2) - mcZ) / aux;
        return grad;
    }

    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double vector_length = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2) + pow(P(2) - mcZ, 2));
        if (vector_length == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project point that's in the center of Brep sphere  ", "");
        }

        Proj(0) = (P(0) - mcX) * mR / vector_length + mcX;
        Proj(1) = (P(1) - mcY) * mR / vector_length + mcY;
        Proj(2) = (P(2) - mcZ) * mR / vector_length + mcZ;

        return 0;
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
    std::string Info() const override
    {
        return "Spherical Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", cZ: " << mcZ << ", R: " << mR;
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

    double mcX, mcY, mcZ, mR;

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
    SphericalLevelSet& operator=(SphericalLevelSet const& rOther);

    ///@}

}; // Class SphericalLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SphericalLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SphericalLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERICAL_LEVEL_SET_H_INCLUDED  defined
