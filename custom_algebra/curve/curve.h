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
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"
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
        KRATOS_THROW_ERROR(std::logic_error, "Error calling abstract function", __FUNCTION__)
    }


    /// inherit from Function
    OutputType GetDerivative(const int& component, const InputType& t) const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling abstract function", __FUNCTION__)
    }


    /// inherit from Function
    OutputType GetSecondDerivative(const int& component_1, const int& component_2, const InputType& t) const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling abstract function", __FUNCTION__)
    }


    /// inherit from Function
    BaseType::Pointer GetDiffFunction(const int& component) const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling abstract function", __FUNCTION__)
    }

    /******************** SPECIFIC CURVE OPERATIONS **********************/

    /// Compute the distance from a point to a curve
    /// Because the bisection algorithm is used, user must notice some additional parameters may be necessary to give successful return
    /// On output, the distance and local coordinates of the projection point are returned
    double ComputeDistance(const PointType& P, double& t) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        double nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double d = ComputeDistanceByBisection(P, t, tmin, tmax, nsampling);
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
        double nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        return this->ProjectOnCurveUsingBisection(P, Proj, t, tmin, tmax, nsampling);
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
        double nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        double t;
        int error_code = this->ComputeIntersectionByBisection(rLevelSet, P, t, tmin, tmax, nsampling);
        noalias(dP) = this->GetDerivative(0, t);
        return error_code;
    }

    /// Compute the intersection of the curve with the level set
    /// In case the point is not found, it will be set to locate in infinity
    PointType ComputeIntersection(const LevelSet& rLevelSet) const
    {
        double tmin = DataValueContainer::GetValue(CURVE_LOWER_BOUND);
        double tmax = DataValueContainer::GetValue(CURVE_UPPER_BOUND);
        double nsampling = DataValueContainer::GetValue(CURVE_NUMBER_OF_SAMPLING);
        PointType IntPoint;
        double t;
        int error_code = this->ComputeIntersectionByBisection(rLevelSet, IntPoint, t, tmin, tmax, nsampling);
        if (error_code != 0)
            std::cout << "WARNING: Error computing the intersection with level set " << rLevelSet.Info()
                      << ", error code = " << error_code << std::endl;
        return IntPoint;
    }

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

    /// Compute the distance from a point to the curve based on Newton Raphson, quick but not always give a solution :(
    /// On output, the distance and local coordinates of the projection point are returned
    double ComputeDistanceByNewtonRaphson(const PointType& P, double& t) const
    {
        const double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);
        const int max_iters = DataValueContainer::GetValue(CURVE_MAX_ITERATIONS);

        // firstly compute the projection of point P to the curve
        int iter = 0;
        PointType Proj, dProj, ddProj;
        double rhs, drhs;

        do
        {
            noalias(Proj) = this->GetValue(t);
            noalias(dProj) = this->GetDerivative(0, t);
            noalias(ddProj) = this->GetSecondDerivative(0, 0, t);
            rhs = inner_prod(dProj, P - Proj);
            if (fabs(rhs) < tol) break;
            drhs = inner_prod(ddProj, P - Proj) - inner_prod(dProj, dProj);
            t -= rhs/drhs;
            ++iter;
            if (iter > max_iters) break;
        }
        while (fabs(rhs) > tol);

        if (iter > max_iters)
        {
            KRATOS_WATCH(P)
            KRATOS_WATCH(this->GetValue(0.0))
            KRATOS_WATCH(this->GetValue(1.0))
            KRATOS_WATCH(t)
            KRATOS_WATCH(rhs)
            KRATOS_THROW_ERROR(std::logic_error, "The local iteration does not converge", "")
        }
        // KRATOS_WATCH(t)
        // KRATOS_WATCH(P)
        // KRATOS_WATCH(Proj)

        // compute the distance
        return norm_2(P - Proj);
    }

    /// Compute the distance from a point to the curve based on bisection algorithm, slow but stable :)
    /// The sampling parameter is for computing the projection, the more it is more stable but take more time
    /// On return, the distance and local coordinates of the projection point are returned
    /// This subroutine require a bound to search for the local projection point. A sampling is required so that at least one segment exists by which the projection function <dP, P'-P> has reversed sign
    double ComputeDistanceByBisection(const PointType& P, double& t, const double& tmin, const double& tmax, const int& nsampling) const
    {
        PointType Proj;
        int stat = ProjectOnCurveUsingBisection(P, Proj, t, tmin, tmax, nsampling);

        // if the bisection can't find the point in [tmin, tmax] region. The projection point will be computed before tmin or after tmax using the tangent information on that point.
        if (stat != 0)
        {
            PointType First = static_cast<PointType>(this->GetValue(tmin));
            PointType dFirst = static_cast<PointType>(this->GetDerivative(0, tmin));

            PointType Last = static_cast<PointType>(this->GetValue(tmax));
            PointType dLast = static_cast<PointType>(this->GetDerivative(0, tmax));

            double test1 = inner_prod(P - First, dFirst);
            double test2 = inner_prod(P - Last, dLast);

            if ((test1 < 0.0) &&  (test2 < 0.0)) // the point is before tmin
            {
                // project on to the line created by First and dFirst
                double a = inner_prod(P - First, dFirst) / inner_prod(dFirst, dFirst);
                noalias(Proj) = First + a*dFirst;
            }
            else if ((test1 > 0.0) &&  (test2 < 0.0)) // the point is in between
            {
                KRATOS_WATCH(P)
                KRATOS_WATCH(stat)
                KRATOS_WATCH(tmin)
                KRATOS_WATCH(tmax)
                KRATOS_WATCH(nsampling)
                KRATOS_WATCH(test1)
                KRATOS_WATCH(test2)
                KRATOS_THROW_ERROR(std::logic_error, "There is error in computing the bisection. ProjectOnCurveUsingBisection shall detect it at the beginning of this function.", "")
            }
            else if ((test1 > 0.0) &&  (test2 > 0.0)) // the point is after tmax
            {
                // project on to the line created by Last and dLast
                double a = inner_prod(P - Last, dLast) / inner_prod(dLast, dLast);
                noalias(Proj) = Last + a*dLast;
            }
            else if ((test1 < 0.0) &&  (test2 > 0.0)) // this can happen if the curve has U-shape
            {
                // project to closer line
                double a = inner_prod(P - First, dFirst) / inner_prod(dFirst, dFirst);
                noalias(Proj) = First + a*dFirst;
                double d1 = norm_2(P - Proj);

                double b = inner_prod(P - Last, dLast) / inner_prod(dLast, dLast);
                noalias(Proj) = Last + b*dLast;
                double d2 = norm_2(P - Proj);

                if (d1 < d2)
                    noalias(Proj) = First + a*dFirst;
                else
                    noalias(Proj) = Last + b*dLast;
            }
        }

        return norm_2(P - Proj);
    }

    /// Compute the projection from a point to the curve using Bisection
    /// The sampling parameter is to provide a hint to predict the good span, where the sign on ending of the span are not same
    /// On return:
    ///     + 0: projection point is found inside the parametric domain of the curve
    ///     + 1: the projection cannot be detected in the given parametric domain and is reset to the left tip
    ///     + 2: the projection cannot be detected in the given parametric domain and is reset to the right tip
    int ProjectOnCurveUsingBisection(const PointType& P, PointType& Proj, double& t, const double& tmin, const double& tmax, const int& nsampling) const
    {
        const double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);

        // firstly do the sampling
        std::vector<double> f(nsampling);
        PointType dProj;
        for (std::size_t i = 0; i < nsampling+1; ++i)
        {
            t = tmin + i*(tmax-tmin)/nsampling;
            noalias(Proj) = this->GetValue(t);
            noalias(dProj) = this->GetDerivative(0, t);
            f[i] = inner_prod(dProj, P - Proj);
        }

        // secondly find the good span
        int stat = -1;
        for (std::size_t i = 0; i < nsampling; ++i)
        {
            if (fabs(f[i]) < tol)
            {
                stat = 0;
                t = tmin + i*(tmax-tmin)/nsampling;
                noalias(Proj) = this->GetValue(t);
                break;
            }

            if (f[i]*f[i+1] < 0.0)
            {
                // found the segment, do the bisection
                double left = tmin + i*(tmax-tmin)/nsampling;
                double right = tmin + (i+1)*(tmax-tmin)/nsampling;
                double mid;
                double fleft = f[i], fright = f[i+1], fmid;
                while ((right - left) > tol)
                {
                    mid = 0.5*(left + right);

                    noalias(Proj) = this->GetValue(mid);
                    noalias(dProj) = this->GetDerivative(0, mid);
                    fmid = inner_prod(dProj, P - Proj);

                    if (fabs(fmid) < tol)
                        break;

                    if (fmid * fleft > 0.0)
                    {
                        left = mid;
                        noalias(Proj) = this->GetValue(left);
                        noalias(dProj) = this->GetDerivative(0, left);
                        fleft = inner_prod(dProj, P - Proj);
                    }

                    if (fmid * fright > 0.0)
                    {
                        right = mid;
                        noalias(Proj) = this->GetValue(right);
                        noalias(dProj) = this->GetDerivative(0, right);
                        fright = inner_prod(dProj, P - Proj);
                    }
                }

                stat = 0;
                t = mid;
                break;
            }
        }

        if (stat != 0)
        {
            if (f[0] < -tol)
            {
                t = tmin; // set default value to t
                noalias(Proj) = this->GetValue(t);
                stat = 1;
            }
            if (f[nsampling-1] > tol)
            {
                t = tmax; // set default value to t
                noalias(Proj) = this->GetValue(t);
                stat = 2;
            }
        }

        // if (!found)
        // {
        //     KRATOS_WATCH(P)
        //     KRATOS_WATCH(this->GetValue(0.0))
        //     KRATOS_WATCH(this->GetValue(1.0))
        //     std::cout << "f:";
        //     for (std::size_t i = 0; i < f.size(); ++i)
        //         std::cout << " " << f[i];
        //     std::cout << std::endl;
        //     KRATOS_THROW_ERROR(std::logic_error, "Bisection error: there are no valid segment", "")
        // }

        return 0;
    }

    /// Compute the intersection of the curve with a level set using Bisection algorithm
    /// The sampling parameter is to provide a hint to predict the good span, where the sign on ending of the span are not same
    /// On return:
    ///     + 0: the intersection point is found inside the parametric domain of the curve
    ///     + 1: the intersection point cannot be detected in the given parametric domain and is reset to the left tip
    ///     + 2: the intersection point cannot be detected in the given parametric domain and is reset to the right tip
    int ComputeIntersectionByBisection(const LevelSet& rLevelSet, PointType& IntPoint, double& t, const double& tmin, const double& tmax, const std::size_t& nsampling) const
    {
        const double tol = DataValueContainer::GetValue(CURVE_SEARCH_TOLERANCE);

        // firstly do the sampling
        std::vector<double> f(nsampling);
        for (std::size_t i = 0; i < nsampling+1; ++i)
        {
            t = tmin + i*(tmax-tmin)/nsampling;
            noalias(IntPoint) = this->GetValue(t);
            f[i] = rLevelSet.GetValue(IntPoint);
        }

        // secondly find the good span
        int stat = -1;
        for (std::size_t i = 0; i < nsampling; ++i)
        {
            if (fabs(f[i]) < tol)
            {
                stat = 0;
                t = tmin + i*(tmax-tmin)/nsampling;
                noalias(IntPoint) = this->GetValue(t);
                break;
            }

            if (f[i]*f[i+1] < 0.0)
            {
                // found the segment, do the bisection
                double left = tmin + i*(tmax-tmin)/nsampling;
                double right = tmin + (i+1)*(tmax-tmin)/nsampling;
                double mid;
                double fleft = f[i], fright = f[i+1], fmid;
                while ((right - left) > tol)
                {
                    mid = 0.5*(left + right);

                    noalias(IntPoint) = this->GetValue(mid);
                    fmid = rLevelSet.GetValue(IntPoint);

                    if (fabs(fmid) < tol)
                        break;

                    if (fmid * fleft > 0.0)
                    {
                        left = mid;
                        noalias(IntPoint) = this->GetValue(left);
                        fleft = rLevelSet.GetValue(IntPoint);
                    }

                    if (fmid * fright > 0.0)
                    {
                        right = mid;
                        noalias(IntPoint) = this->GetValue(right);
                        fright = rLevelSet.GetValue(IntPoint);
                    }
                }

                stat = 0;
                t = mid;
                break;
            }
        }

        if (stat != 0)
        {
            if (f[0] < -tol)
            {
                t = tmin; // set default value to t
                noalias(IntPoint) = this->GetValue(t);
                stat = 1;
            }
            if (f[nsampling-1] > tol)
            {
                t = tmax; // set default value to t
                noalias(IntPoint) = this->GetValue(t);
                stat = 2;
            }
        }

        // if (!found)
        // {
        //     KRATOS_WATCH(P)
        //     KRATOS_WATCH(this->GetValue(0.0))
        //     KRATOS_WATCH(this->GetValue(1.0))
        //     std::cout << "f:";
        //     for (std::size_t i = 0; i < f.size(); ++i)
        //         std::cout << " " << f[i];
        //     std::cout << std::endl;
        //     KRATOS_THROW_ERROR(std::logic_error, "Bisection error: there are no valid segment", "")
        // }

        return 0;

    }

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
