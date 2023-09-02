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
//  Date:            10 Feb 2017
//

#if !defined(KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED

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

///
/**
 * Level set representing the signed distance function to a circle
 * f = sqrt((X - cx)^2 + (Y - cy)^2) - R
 */
class CircularLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CircularLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(CircularLevelSet);

    typedef LevelSet BaseType;

#if defined(__clang__) || defined(__INTEL_COMPILER)
    static constexpr double PI = 3.1415926535897932384626433832795028841971693;
#elif defined(__GNUC__) || defined(__GNUG__)
    static constexpr double PI = std::atan(1.0) * 4;
#else
    static constexpr double PI = std::atan(1.0) * 4;
#endif

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircularLevelSet(const double& cX, const double& cY, const double& R)
        : BaseType(), mcX(cX), mcY(cY), mR(R)
    {}

    /// Copy constructor.
    CircularLevelSet(CircularLevelSet const& rOther)
        : BaseType(rOther), mcX(rOther.mcX), mcY(rOther.mcY), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~CircularLevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const override
    {
        return LevelSet::Pointer(new CircularLevelSet(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return 2;
    }

    double GetValue(const PointType& P) const override
    {
        return sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2)) - mR;
    }

    Vector GetGradient(const PointType& P) const override
    {
        Vector grad(3);
        double aux = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2));
        grad(0) = (P(0) - mcX) / aux;
        grad(1) = (P(1) - mcY) / aux;
        grad(2) = 0.0;
        return grad;
    }

    Matrix GetGradientDerivatives(const PointType& P) const override
    {
        Matrix Jac(3, 3);
        noalias(Jac) = ZeroMatrix(3, 3);

        double aux = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2));

        Jac(0, 0) = 1.0 / aux - pow(P(0) - mcX, 2) / pow(aux, 3);
        Jac(0, 1) = -(P(0) - mcX) * (P(1) - mcY) / pow(aux, 3);

        Jac(1, 0) = Jac(0, 1);
        Jac(1, 1) = 1.0 / aux - pow(P(1) - mcY, 2) / pow(aux, 3);

        return Jac;
    }

    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double vector_length = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2));
        if (vector_length == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project node that's in the center of Brep circle", "");
        }

        Proj(0) = (P(0) - mcX) * mR / vector_length + mcX;
        Proj(1) = (P(1) - mcY) * mR / vector_length + mcY;
        Proj(2) = 0.0;

        return 0;
    }

    /// compute the derivatives of the projection point w.r.t to the original point.
    void ProjectionDerivatives(const PointType& P, Matrix& Derivatives) const final
    {
        if (Derivatives.size1() != 3 || Derivatives.size2() != 3)
        {
            Derivatives.resize(3, 3, false);
        }

        double vector_length = sqrt(pow(P(0) - mcX, 2) + pow(P(1) - mcY, 2));
        if (vector_length == 0)
        {
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project node that's in the center of Brep circle", "");
        }

        noalias(Derivatives) = ZeroMatrix(3, 3);

        Vector dvector_length(2);
        dvector_length(0) = 2.0 * (P(0) - mcX);
        dvector_length(1) = 2.0 * (P(1) - mcY);

        Derivatives(0, 0) = mR / vector_length - (P(0) - mcX) * mR * dvector_length(0) / pow(vector_length, 2);
        Derivatives(0, 1) = -(P(0) - mcX) * mR * dvector_length(1) / pow(vector_length, 2);

        Derivatives(1, 0) = -(P(1) - mcY) * mR * dvector_length(0) / pow(vector_length, 2);
        Derivatives(1, 1) = mR / vector_length - (P(1) - mcY) * mR * dvector_length(1) / pow(vector_length, 2);
    }

    /***********EXCLUSIVE INTERFACE****************/

    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<PointType>& radial_points,
                        const double& start_angle, const double& end_angle,
                        const std::size_t& nsampling_radial) const
    {
        radial_points.resize(nsampling_radial);
        double small_angle = (end_angle - start_angle) / nsampling_radial;

        PointType V;
        V[2] = 0.0;
        double d;
        for (std::size_t j = 0; j < nsampling_radial; ++j)
        {
            d = start_angle + j * small_angle;
            V[0] = mcX + mR * std::cos(d);
            V[1] = mcY + mR * std::sin(d);
            noalias(radial_points[j]) = V;
        }
    }

    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<PointType>& radial_points,
                        const std::size_t& nsampling_radial) const
    {
        this->GeneratePoints(radial_points, 0.0, 2 * PI, nsampling_radial);
    }

    /// Create the elements based on sampling points on the line
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateLineElements(ModelPart& r_model_part,
            const std::string& sample_element_name,
            Properties::Pointer pProperties,
            const double& start_angle,
            const double& end_angle,
            const std::size_t& nsampling_radial,
            const bool close = false) const
    {
        // firstly create the sampling points on surface
        std::vector<PointType> sampling_points;
        this->GeneratePoints(sampling_points, start_angle, end_angle, nsampling_radial);
        int type = 1; // linear element
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateLineElements(r_model_part, sampling_points, sample_element_name, type, close, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
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
        return "Circular Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "cX: " << mcX << ", cY: " << mcY << ", R: " << mR;
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

    double mcX, mcY, mR;

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
    CircularLevelSet& operator=(CircularLevelSet const& rOther);

    ///@}

}; // Class CircularLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, CircularLevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const CircularLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CIRCULAR_LEVEL_SET_H_INCLUDED  defined
