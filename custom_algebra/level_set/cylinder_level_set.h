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


#if !defined(KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/brep_mesh_utility.h"

/*#define PI 3.1415926535897932384626433832795028841971693*/

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
class CylinderLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CylinderLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(CylinderLevelSet);

    typedef LevelSet BaseType;

    static constexpr double PI = std::atan(1.0)*4;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CylinderLevelSet(const double& cX, const double& cY, const double& cZ, const double& dX, const double& dY, const double& dZ, const double& R)
    : BaseType(), mcX(cX), mcY(cY), mcZ(cZ), mR(R)
    {
        mLength = sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));

        if(mLength == 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The director vector can't be null", "")

        mdX = dX / mLength;
        mdY = dY / mLength;
        mdZ = dZ / mLength;
    }

    /// Copy constructor.
    CylinderLevelSet(CylinderLevelSet const& rOther)
    : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mcZ(rOther.mcZ)
    , mdX(rOther.mdX), mdY(rOther.mdY), mdZ(rOther.mdZ), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~CylinderLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new CylinderLevelSet(*this));
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


    /// Generate the sampling points on the level set surface
    std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
        const double& start_angle, const double& end_angle, const double& tmin, const double& tmax) const
    {
        // KRATOS_WATCH(nsampling_axial)
        // KRATOS_WATCH(nsampling_radial)

        std::vector<std::vector<PointType> > results;
        double small_angle = (end_angle - start_angle) / nsampling_radial;

        double t, d;
        PointType P, T, B, V, Up;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        for (std::size_t i = 0; i < nsampling_axial; ++i)
        {
            t = tmin + i*(tmax-tmin)/(nsampling_axial-1);

            P[0] = mcX + t*mdX*mLength;
            P[1] = mcY + t*mdY*mLength;
            P[2] = mcZ + t*mdZ*mLength;
            // KRATOS_WATCH(P)

            T[0] = mdX;
            T[1] = mdY;
            T[2] = mdZ;
            // KRATOS_WATCH(T)

            noalias(B) = MathUtils<double>::CrossProduct(Up, T);
            // KRATOS_WATCH(B)

            std::vector<PointType> radial_points(nsampling_radial);
            for (std::size_t j = 0; j < nsampling_radial; ++j)
            {
                d = start_angle + j*small_angle;
                noalias(V) = std::cos(d)*Up + std::sin(d)*B;

                noalias(radial_points[j]) = P + mR*V;
            }

            results.push_back(radial_points);
        }

        return results;
    }


    /// Generate the sampling points on the level set surface
    std::vector<std::vector<PointType> > GeneratePoints(const std::size_t& nsampling_axial, const std::size_t& nsampling_radial) const
    {
        return GeneratePoints(nsampling_axial, nsampling_radial, 0.0, 2*PI, 0.0, 1.0);
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4ElementsClosedLoop(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial) const
    {
        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points = this->GeneratePoints(nsampling_axial, nsampling_radial);
        int order = 1;
        int close_dir = 2;
        int activation_dir = 1;
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4ElementsClosedLoop(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial,
        const double& tmin,
        const double& tmax) const
    {
        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points = this->GeneratePoints(nsampling_axial, nsampling_radial, 0.0, 2*PI, tmin, tmax);
        int order = 1;
        int close_dir = 2;
        int activation_dir = 1;
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4Elements(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial,
        const double& start_radial_angle,
        const double& end_radial_angle) const
    {
        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points = this->GeneratePoints(nsampling_axial, nsampling_radial, start_radial_angle, end_radial_angle, 0.0, 1.0);
        int order = 1;
        int close_dir = 0;
        int activation_dir = 1;
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }

    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
        double vector_length = sqrt(pow(P(0)-pX, 2) + pow(P(1)-pY, 2) + pow(P(2)-pZ, 2));
        if (vector_length == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project point that's on the center line  ", "");

        Proj(0) = (P(0) - pX) * mR / vector_length + pX;
        Proj(1) = (P(1) - pY) * mR / vector_length + pY;
        Proj(2) = (P(2) - pZ) * mR / vector_length + pZ;

        return 0;
    }

    /// compute the derivatives of the projection point w.r.t to the original point.
    /// The derivatives are organized as;
    ///     [d Proj[0] / d P[0], d Proj[0] / d P[1], d Proj[0] / d P[2]]
    ///     [d Proj[1] / d P[0], d Proj[1] / d P[1], d Proj[1] / d P[2]]
    ///     [d Proj[2] / d P[0], d Proj[2] / d P[1], d Proj[2] / d P[2]]
    void ProjectionDerivatives(const PointType& P, Matrix& Derivatives) const final
    {
        if (Derivatives.size1() != 3 || Derivatives.size2() != 3)
            Derivatives.resize(3, 3, false);

        double t = (P(0) - mcX) * mdX + (P(1) - mcY) * mdY + (P(2) - mcZ) * mdZ;
        double pX = mcX + t*mdX;
        double pY = mcY + t*mdY;
        double pZ = mcZ + t*mdZ;
        double vector_length = sqrt(pow(P(0)-pX, 2) + pow(P(1)-pY, 2) + pow(P(2)-pZ, 2));
        if (vector_length == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project point that's on the center line  ", "");

        Vector dt(3);
        dt(0) = mdX; dt(1) = mdY; dt(2) = mdZ;

        Vector dpX(3), dpY(3), dpZ(3);

        noalias(dpX) = mdX*dt;
        noalias(dpY) = mdY*dt;
        noalias(dpZ) = mdZ*dt;

        Vector dvector_length(3);
        dvector_length(0) = ((P(0) - pX)*(1.0 - dpX(0)) + (P(1)-pY)*(-dpY(0)) + (P(2)-pZ)*(-dpZ(0))) / vector_length;
        dvector_length(1) = ((P(0) - pX)*(-dpX(1)) + (P(1)-pY)*(1.0 - dpY(1)) + (P(2)-pZ)*(-dpZ(1))) / vector_length;
        dvector_length(2) = ((P(0) - pX)*(-dpX(2)) + (P(1)-pY)*(-dpY(2)) + (P(2)-pZ)*(1.0 - dpZ(2))) / vector_length;

        Derivatives(0, 0) = (1.0-dpX(0))*mR/vector_length + dpX(0) - (P(0)-pX)*mR*dvector_length(0)/pow(vector_length, 2);
        Derivatives(0, 1) = (-dpX(1))*mR/vector_length + dpX(1) - (P(0)-pX)*mR*dvector_length(1)/pow(vector_length, 2);
        Derivatives(0, 2) = (-dpX(2))*mR/vector_length + dpX(2) - (P(0)-pX)*mR*dvector_length(2)/pow(vector_length, 2);

        Derivatives(1, 0) = (-dpY(0))*mR/vector_length + dpY(0) - (P(1)-pY)*mR*dvector_length(0)/pow(vector_length, 2);
        Derivatives(1, 1) = (1.0-dpY(1))*mR/vector_length + dpY(1) - (P(1)-pY)*mR*dvector_length(1)/pow(vector_length, 2);
        Derivatives(1, 2) = (-dpY(2))*mR/vector_length + dpY(2) - (P(1)-pY)*mR*dvector_length(2)/pow(vector_length, 2);

        Derivatives(2, 0) = (-dpZ(0))*mR/vector_length + dpZ(0) - (P(2)-pZ)*mR*dvector_length(0)/pow(vector_length, 2);
        Derivatives(2, 1) = (-dpZ(1))*mR/vector_length + dpZ(1) - (P(2)-pZ)*mR*dvector_length(1)/pow(vector_length, 2);
        Derivatives(2, 2) = (1.0-dpZ(2))*mR/vector_length + dpZ(2) - (P(2)-pZ)*mR*dvector_length(2)/pow(vector_length, 2);
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
        return "Cylinder Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", cZ: " << mcZ
                 << "dX: " << mdX << ", dY: " << mdY << ", dZ: " << mdZ
                 << ", R: " << mR;
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
    double mR;


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
    CylinderLevelSet& operator=(CylinderLevelSet const& rOther);

    ///@}

}; // Class CylinderLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                CylinderLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CylinderLevelSet& rThis)
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

#endif // KRATOS_CYLINDER_LEVEL_SET_H_INCLUDED  defined
