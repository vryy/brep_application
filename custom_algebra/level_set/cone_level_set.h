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


#if !defined(KRATOS_CONE_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CONE_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"
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

/// Short class definition.
/** Detail class definition.
*/
class ConeLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConeLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(ConeLevelSet);

    typedef LevelSet BaseType;

    #if defined(__clang__) || defined(__INTEL_COMPILER)
    static constexpr double PI = 3.1415926535897932384626433832795028841971693;
    #elif defined(__GNUC__) || defined(__GNUG__)
    static constexpr double PI = std::atan(1.0)*4;
    #else
    static constexpr double PI = std::atan(1.0)*4;
    #endif

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConeLevelSet(const double& cX, const double& cY, const double& cZ, const double& dX, const double& dY, const double& dZ, const double& phi)
    : BaseType(), mcX(cX), mcY(cY), mcZ(cZ), mphi(phi)
    {
        mLength = sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));

        if(mLength == 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The director vector can't be null", "")

        mdX = dX / mLength;
        mdY = dY / mLength;
        mdZ = dZ / mLength;
    }

    /// Copy constructor.
    ConeLevelSet(ConeLevelSet const& rOther)
    : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mcZ(rOther.mcZ)
    , mdX(rOther.mdX), mdY(rOther.mdY), mdZ(rOther.mdZ), mphi(rOther.mphi)
    {}

    /// Destructor.
    virtual ~ConeLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new ConeLevelSet(*this));
    }


    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }


    double GetValue(const PointType& P) const final
    {
        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;

        return pow(P(0) - pX, 2) + pow(P(1) - pY, 2) + pow(P(2) - pZ, 2) - pow(t * std::tan(mphi*PI/180), 2);
    }


    // virtual Vector GetGradient(const PointType& P) const
    // {
    //     double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
    //     double pX = mcX + t*mdX;
    //     double pY = mcY + t*mdY;
    //     double pZ = mcZ + t*mdZ;
    //     Vector grad(3);
    //
    //     return grad;
    // }


    // /// Generate the sampling points on the level set surface
    // std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
    //     const double& start_angle, const double& end_angle) const
    // {
    //
    // }


    // /// Generate the sampling points on the level set surface
    // std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial) const
    // {
    //     return GeneratePoints(nsampling_axial, nsampling_radial, 0.0, 2*PI);
    // }

    /// projects a point on the surface of level_set (this does not project perpendicularly instead it projects in the radial direction)
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
        double vector_length = sqrt(pow(P(0)-pX, 2) + pow(P(1)-pY, 2) + pow(P(2)-pZ, 2));
        if (vector_length == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project point that's on the center line  ", "");

        Proj(0) = (P(0) - pX) * (t * std::tan(mphi*PI/180)) / vector_length + pX;
        Proj(1) = (P(1) - pY) * (t * std::tan(mphi*PI/180)) / vector_length + pY;
        Proj(2) = (P(2) - pZ) * (t * std::tan(mphi*PI/180)) / vector_length + pZ;

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
    std::string Info() const final
    {
        return "Cone Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", cZ: " << mcZ
                 << "dX: " << mdX << ", dY: " << mdY << ", dZ: " << mdZ
                 << ", phi: " << mphi;
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


    double mcX, mcY, mcZ; // point on center line
    double mdX, mdY, mdZ; // director vector
    double mLength;
    double mphi;


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
    ConeLevelSet& operator=(ConeLevelSet const& rOther);

    ///@}

}; // Class ConeLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ConeLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ConeLevelSet& rThis)
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

#endif // KRATOS_CONE_LEVEL_SET_H_INCLUDED  defined
