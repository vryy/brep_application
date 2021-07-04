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
//  Date:            1 Feb 2018
//


#if !defined(KRATOS_DISTANCE_TO_CURVE_LEVEL_SET_H_INCLUDED )
#define  KRATOS_DISTANCE_TO_CURVE_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_algebra/curve/curve.h"
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
 * Level Set computes the distance to a curve
 */
class DistanceToCurveLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceToCurveLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(DistanceToCurveLevelSet);

    typedef LevelSet BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DistanceToCurveLevelSet(const Curve::Pointer pAlignCurve, const double& R)
    : BaseType(), mpCurve(pAlignCurve), mR(R)
    {}

    /// Copy constructor.
    DistanceToCurveLevelSet(DistanceToCurveLevelSet const& rOther)
    : BaseType(rOther), mpCurve(rOther.mpCurve->Clone()), mR(rOther.mR)
    {}

    /// Destructor.
    virtual ~DistanceToCurveLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new DistanceToCurveLevelSet(*this));
    }


    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }


    const Curve& AlignmentCurve() const
    {
        return *mpCurve;
    }


    double Radius() const
    {
        return mR;
    }


    double GetValue(const PointType& P) const final
    {
        return mpCurve->ComputeDistance(P) - mR;
    }


    Vector GetGradient(const PointType& P) const final
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "Not yet implemented")
    }


    /// projects a point on the surface of level_set using Bisection
    /// inherit from LevelSet
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        mpCurve->ProjectOnCurve(P, Proj);

        if (P(0) == Proj(0) && P(1) == Proj(1) && P(2) == Proj(2))
            KRATOS_THROW_ERROR(std::invalid_argument, "trying to project point that's on the curve of Brep distance_to_curve  ", "");

        noalias(Proj) = (P - Proj) * mR / norm_2(P - Proj) + Proj;

        return 0;
    }


    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<std::vector<PointType> >& results, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
        const double& start_angle, const double& end_angle, const double& tmin, const double& tmax) const
    {
        // KRATOS_WATCH(nsampling_axial)
        // KRATOS_WATCH(nsampling_radial)
        // KRATOS_WATCH(tmin)
        // KRATOS_WATCH(tmax)

        double small_angle = (end_angle - start_angle) / nsampling_radial;

        double t, d;
        PointType P, T, N, B, V, Up, Aux;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        results.resize(nsampling_axial);
        for (std::size_t i = 0; i < nsampling_axial; ++i)
        {
            t = tmin + i*(tmax-tmin)/(nsampling_axial-1);
            // KRATOS_WATCH(t)

            noalias(P) = mpCurve->GetValue(t);
            // KRATOS_WATCH(P)

            noalias(T) = mpCurve->GetDerivative(0, t);
            T *= 1.0/norm_2(T);
            // KRATOS_WATCH(T)

            /****************
            Using Frenet frame can cause twist
            // noalias(N) = mpCurve->GetSecondDerivative(0, 0, t);
            // N *= 1.0/norm_2(N);
            // KRATOS_WATCH(N)

            // noalias(B) = MathUtils<double>::CrossProduct(T, N);
            // KRATOS_WATCH(B)
            *****************/

            /****************/
            // this method is better, although the up vector must be defined, see nrbsweep.m
            if (i == 0)
            {
                noalias(Aux) = MathUtils<double>::CrossProduct(Up, T);
                noalias(B) = Aux / norm_2(Aux);
            }
            else
            {
                noalias(Aux) = B - inner_prod(B, T)*T;
                noalias(B) = Aux / norm_2(Aux);
            }
            // KRATOS_WATCH(B)

            noalias(N) = MathUtils<double>::CrossProduct(B, T);
            // KRATOS_WATCH(N)
            /****************/

            results[i].resize(nsampling_radial);
            for (std::size_t j = 0; j < nsampling_radial; ++j)
            {
                d = start_angle + j*small_angle;
                noalias(V) = std::cos(d)*N + std::sin(d)*B;

                noalias(results[i][j]) = P + mR*V;
            }
        }
    }


    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<std::vector<PointType> >& results, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
        const double& start_angle, const double& end_angle) const
    {
        this->GeneratePoints(results, nsampling_axial, nsampling_radial, start_angle, end_angle, 0.0, 1.0);
    }


    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<std::vector<PointType> >& results, const std::size_t& nsampling_axial, const std::size_t& nsampling_radial) const
    {
        const double Pi = 3.1415926535897932384626433832795028841971693;
        this->GeneratePoints(results, nsampling_axial, nsampling_radial, 0.0, 2*Pi);
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
        const double Pi = 3.1415926535897932384626433832795028841971693;

        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points;
        this->GeneratePoints(sampling_points, nsampling_axial, nsampling_radial, 0.0, 2*Pi, tmin, tmax);
        int order = 1;
        int close_dir = 2;
        int activation_dir = 1;
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points,
            sample_element_name, order, close_dir, activation_dir, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }


    /// Create the conditions based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> CreateQ4ConditionsClosedLoop(ModelPart& r_model_part,
        const std::string& sample_condition_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial,
        const double& tmin,
        const double& tmax,
        const bool& reverse) const
    {
        const double Pi = 3.1415926535897932384626433832795028841971693;

        // firstly create the sampling points on surface
        std::vector<std::vector<PointType> > sampling_points;
        this->GeneratePoints(sampling_points, nsampling_axial, nsampling_radial, 0.0, 2*Pi, tmin, tmax);
        int order = 1;
        int close_dir = 2;
        int activation_dir = 1;
        int initial_activation_level = 0;
        BRepMeshUtility::ConditionMeshInfoType Info = BRepMeshUtility::CreateQuadConditions(r_model_part, sampling_points,
            sample_condition_name, order, close_dir, activation_dir, initial_activation_level, reverse, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
    }


    /// Create the elements based on sampling points on the surface
    std::pair<ModelPart::NodesContainerType, ModelPart::ElementsContainerType> CreateQ4ElementsClosedLoop(ModelPart& r_model_part,
        const std::string& sample_element_name,
        Properties::Pointer pProperties,
        const std::size_t& nsampling_axial,
        const std::size_t& nsampling_radial) const
    {
        return CreateQ4ElementsClosedLoop(r_model_part, sample_element_name, pProperties, nsampling_axial, nsampling_radial, 0.0, 1.0);
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
        std::vector<std::vector<PointType> > sampling_points;
        this->GeneratePoints(sampling_points, nsampling_axial, nsampling_radial, start_radial_angle, end_radial_angle);
        int order = 1;
        int close_dir = 0;
        int activation_dir = 1;
        BRepMeshUtility::ElementMeshInfoType Info = BRepMeshUtility::CreateQuadElements(r_model_part, sampling_points, sample_element_name, order, close_dir, activation_dir, pProperties);
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
    std::string Info() const final
    {
        return "Distance to Curve Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "curve: " << *mpCurve;
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


    const Curve::Pointer mpCurve;
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
    DistanceToCurveLevelSet& operator=(DistanceToCurveLevelSet const& rOther);

    ///@}

}; // Class DistanceToCurveLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                DistanceToCurveLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const DistanceToCurveLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SWEEP_LEVEL_SET_H_INCLUDED  defined
