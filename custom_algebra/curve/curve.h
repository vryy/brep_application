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
//  Date:            10 Feb 2020
//

#if !defined(KRATOS_CURVE_H_INCLUDED )
#define  KRATOS_CURVE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/level_set/level_set.h"
#include "brep_application_variables.h"

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
/** Abstract class for a curve in 3D
 * A curve is mapping: t \in R^1 -> P \in R^3
 */
class Curve : public FunctionR1R3, public DataValueContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Curve
    KRATOS_CLASS_POINTER_DEFINITION(Curve);

    typedef FunctionR1R3 BaseType;

    typedef BaseType::InputType InputType; // double

    typedef BaseType::OutputType OutputType; // array_1d<double, 3>

    typedef LevelSet::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Curve() : BaseType(), DataValueContainer()
    {
        DataValueContainer::SetValue(CURVE_SEARCH_TOLERANCE, 1.0e-10);
        DataValueContainer::SetValue(CURVE_MAX_ITERATIONS, 300);
    }

    /// Copy constructor.
    Curve(Curve const& rOther) : BaseType(rOther), DataValueContainer(rOther)
    {}

    /// Destructor.
    virtual ~Curve() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// inherit from Function
    BaseType::Pointer CloneFunction() const override
    {
        return BaseType::Pointer(new Curve(*this));
    }

    /// Curve cloning
    virtual Curve::Pointer Clone() const
    {
        return Curve::Pointer(new Curve(*this));
    }

    /// inherit from Function
    OutputType GetValue(const InputType& t) const override
    {
        KRATOS_ERROR << "Error calling abstract function";
    }

    /// inherit from Function
    OutputType GetDerivative(const int& component, const InputType& t) const override
    {
        KRATOS_ERROR << "Error calling abstract function";
    }

    /// inherit from Function
    OutputType GetSecondDerivative(const int& component_1, const int& component_2, const InputType& t) const override
    {
        KRATOS_ERROR << "Error calling abstract function";
    }

    /// inherit from Function
    BaseType::Pointer GetDiffFunction(const int& component) const override
    {
        KRATOS_ERROR << "Error calling abstract function";
    }

    /******************** SPECIFIC CURVE OPERATIONS **********************/

    /// Compute the length of the curve
    double ComputeLength() const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        // return this->ComputeLengthByTrapezoidalQuadrature(tmin, tmax, nsampling);
        return this->ComputeLengthByGaussQuadrature(tmin, tmax, 2, nsampling);
    }

    /// Compute the length of the curve on a curve section
    double ComputeLength(const double tmin, const double tmax) const
    {
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        // return this->ComputeLengthByTrapezoidalQuadrature(tmin, tmax, nsampling);
        return this->ComputeLengthByGaussQuadrature(tmin, tmax, 2, nsampling);
    }

    /// Compute the distance from a point to a curve
    /// Because the bisection algorithm is used, user must notice some additional parameters may be necessary to give successful return
    /// On output, the distance and local coordinates of the projection point are returned
    double ComputeDistance(const PointType& P, double& t) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);
        double d = ComputeDistanceByBisection(P, t, tmin, tmax, nsampling, tol);
        return d;
    }

    /// Compute the distance from a point to a curve
    /// Because the bisection algorithm is used, user must notice some additional parameters may be necessary to give successful return
    double ComputeDistance(const PointType& P) const
    {
        double t;
        return this->ComputeDistance(P, t);
    }

    /// Compute the projection on the the curve
    /// Because the bisection algorithm is used, user must notice some additional parameters may be necessary to give successful return
    int ProjectOnCurve(const PointType& P, PointType& Proj, double& t) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);
        return this->ProjectOnCurveUsingBisection(P, Proj, t, tmin, tmax, nsampling, tol);
    }

    /// Compute the projection on the the curve
    /// Because the bisection algorithm is used, user must notice some additional parameters may be necessary to give successful return
    int ProjectOnCurve(const PointType& P, PointType& Proj) const
    {
        double t;
        return this->ProjectOnCurve(P, Proj, t);
    }

    /// Compute the projection on the the curve
    /// In case the projection does not exist, either point on two side will be returned
    PointType ComputeProjection(const PointType& P) const
    {
        PointType Proj;
        this->ProjectOnCurve(P, Proj);
        return Proj;
    }

    /// Compute the intersection of the curve with the level set
    int ComputeIntersection(const LevelSet& rLevelSet, PointType& P, PointType& dP) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);
        double t;
        int error_code = this->ComputeIntersectionByBisection(rLevelSet, P, t, tmin, tmax, nsampling, tol);
        noalias(dP) = this->GetDerivative(0, t);
        return error_code;
    }

    /// Compute the intersection of the curve with the level set
    /// In case the point is not found, it will be set to locate in infinity
    PointType ComputeIntersection(const LevelSet& rLevelSet) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        int nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);
        PointType IntPoint;
        double t;
        int error_code = this->ComputeIntersectionByBisection(rLevelSet, IntPoint, t, tmin, tmax, nsampling, tol);
        if (error_code != 0)
            std::cout << "WARNING: Error computing the intersection with level set " << rLevelSet.Info()
                      << ", error code = " << error_code << std::endl;
        return IntPoint;
    }

    /// Compute a uniform distribution on the curve. The distance criteria is the arc length.
    void ComputeUniformDivision(const double tmin, const double tmax,
                const unsigned int num_division, std::vector<double>& tvec, const double tolerance) const;

    /// Compute a division of parameters to have equal distance sampling points on curve
    /// This algorithm is recursive and is very slow when num_division > 6
    void ComputeEquallyDistanceDivisionByRecursiveBisection(const double tmin, const double tmax,
                const unsigned int num_division, std::vector<double>& tvec, const double tolerance) const;

    /******************** END OF SPECIFIC CURVE OPERATIONS **********************/

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
        return "Curve";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    /// Compute the length of the curve using trapezoidal integration rule
    double ComputeLengthByTrapezoidalQuadrature(const double tmin, const double tmax, const std::size_t nsampling) const;

    /// Compute the length of the curve using Gauss integration rule
    double ComputeLengthByGaussQuadrature(const double tmin, const double tmax, const int order, const std::size_t nsampling) const;

    /// Compute the length of the curve using specific integration rule
    template<class TQuadratureType>
    double ComputeLengthByQuadrature(const double tmin, const double tmax, const std::size_t nsampling) const;

    /// Compute the distance from a point to the curve based on Newton Raphson, quick but not always give a solution :(
    /// On output, the distance and local coordinates of the projection point are returned
    double ComputeDistanceByNewtonRaphson(const PointType& P, double& t,
            const double tol, const unsigned int max_iters) const;

    /// Compute the distance from a point to the curve based on bisection algorithm, slow but stable :)
    /// The sampling parameter is for computing the projection, the more it is more stable but take more time
    /// On return, the distance and local coordinates of the projection point are returned
    /// This subroutine require a bound to search for the local projection point.
    /// A sampling is required so that at least one segment exists by which the projection function <dP, P'-P> has reversed sign
    double ComputeDistanceByBisection(const PointType& P, double& t, const double tmin, const double tmax, const int nsampling, const double tol) const;

    /// Compute the projection from a point to the curve using Bisection
    /// The sampling parameter is to provide a hint to predict the good span, where the sign on ending of the span are not same
    /// On return:
    ///     + 0: projection point is found inside the parametric domain of the curve
    ///     + 1: the projection cannot be detected in the given parametric domain and is reset to the left tip
    ///     + 2: the projection cannot be detected in the given parametric domain and is reset to the right tip
    int ProjectOnCurveUsingBisection(const PointType& P, PointType& Proj, double& t,
            const double tmin, const double tmax, const int nsampling, const double tol) const;

    /// Compute the intersection of the curve with a level set using Bisection algorithm
    /// The sampling parameter is to provide a hint to predict the good span, where the sign on ending of the span are not same
    /// On return:
    ///     + 0: the intersection point is found inside the parametric domain of the curve
    ///     + 1: the intersection point cannot be detected in the given parametric domain and is reset to the left tip
    ///     + 2: the intersection point cannot be detected in the given parametric domain and is reset to the right tip
    int ComputeIntersectionByBisection(const LevelSet& rLevelSet, PointType& IntPoint, double& t,
            const double tmin, const double tmax, const std::size_t nsampling, const double tol) const;

    /// Compute the parameter of the "middle" point on curve which has the same distance to left and right
    void ProjectMiddleUsingBisection(const double tmin, const double tmax, double& tmid, const double tolerance) const;

    /// Compute the point tmid which makes the curve section from [tmin,tmid] equal to arc length
    void FindPointOnCurveWithArcLengthUsingBisection(const double tmin, const double tmax, double& t,
            const double arc_length, const double tolerance) const;

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
    Curve& operator=(Curve const& rOther);

    ///@}

}; // Class Curve

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, Curve& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const Curve& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CURVE_H_INCLUDED  defined
