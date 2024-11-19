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

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "custom_algebra/curve/curve.h"

namespace Kratos
{

double Curve::ComputeLengthByTrapezoidalQuadrature(const double tmin, const double tmax, const std::size_t nsampling) const
{
    const double h = (tmax - tmin) / nsampling;

    std::vector<double> f(nsampling + 1);
    for (std::size_t i = 0; i < nsampling+1; ++i)
    {
        const double t = tmin + i*h;

        const OutputType d = this->GetDerivative(0, t);

        f[i] = norm_2(d);
    }

    double length = 0.0;
    for (std::size_t i = 0; i < nsampling; ++i)
    {
        length += 0.5*h*(f[i] + f[i+1]);
    }

    return length;
}

double Curve::ComputeLengthByGaussQuadrature(const double tmin, const double tmax, const int order, const std::size_t nsampling) const
{
    if (order == 1)
    {
        return ComputeLengthByQuadrature<LineGaussLegendreIntegrationPoints1>(tmin, tmax, nsampling);
    }
    else if (order == 2)
    {
        return ComputeLengthByQuadrature<LineGaussLegendreIntegrationPoints2>(tmin, tmax, nsampling);
    }
    else if (order == 3)
    {
        return ComputeLengthByQuadrature<LineGaussLegendreIntegrationPoints3>(tmin, tmax, nsampling);
    }
    else
        KRATOS_ERROR << "Unsupported order " << order;
}

double Curve::ComputeLengthRecursively(const double tmin, const double tmax, const double tolerance) const
{
    // start value for number of sampling
    std::size_t nsampling = 5;
    const int order = 3;
    double prev_length = this->ComputeLengthByGaussQuadrature(tmin, tmax, order, nsampling);

    bool converged = false;
    while (!converged)
    {
        nsampling += 5;
        double length = this->ComputeLengthByGaussQuadrature(tmin, tmax, order, nsampling);
        if (std::abs(length - prev_length) < tolerance) converged = true;
        prev_length = length;
    }

    return prev_length;
}

template<class TQuadratureType>
double Curve::ComputeLengthByQuadrature(const double tmin, const double tmax, const std::size_t nsampling) const
{
    const double h = (tmax - tmin) / nsampling;

    const auto& integration_points = TQuadratureType::IntegrationPoints();

    double length = 0.0;
    for (std::size_t i = 0; i < nsampling; ++i)
    {
        const double t1 = tmin + i*h;
        const double t2 = tmin + (i+1)*h;

        for (std::size_t j = 0; j < integration_points.size(); ++j)
        {
            const double a = 0.5*(integration_points[j].X() + 1.0);
            const double w = integration_points[j].Weight();
            const double t = (1-a)*t1 + a*t2;
            const OutputType d = this->GetDerivative(0, t);
            const double f = norm_2(d);
            length += 0.5*h*w*f;
        }
    }

    return length;
}

double Curve::ComputeDistanceByNewtonRaphson(const PointType& P, double& t,
        const double tol, const unsigned int max_iters) const
{
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
        if (fabs(rhs) < tol) { break; }
        drhs = inner_prod(ddProj, P - Proj) - inner_prod(dProj, dProj);
        t -= rhs / drhs;
        ++iter;
        if (iter > max_iters) { break; }
    }
    while (fabs(rhs) > tol);

    if (iter > max_iters)
    {
        KRATOS_WATCH(P)
        KRATOS_WATCH(this->GetValue(0.0))
        KRATOS_WATCH(this->GetValue(1.0))
        KRATOS_WATCH(t)
        KRATOS_WATCH(rhs)
        KRATOS_ERROR << "The local iteration does not converge";
    }
    // KRATOS_WATCH(t)
    // KRATOS_WATCH(P)
    // KRATOS_WATCH(Proj)

    // compute the distance
    return norm_2(P - Proj);
}

double Curve::ComputeDistanceByBisection(const PointType& P, double& t, const double tmin, const double tmax, const int nsampling, const double tol) const
{
    PointType Proj;
    int stat = ProjectOnCurveUsingBisection(P, Proj, t, tmin, tmax, nsampling, tol);

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
            noalias(Proj) = First + a * dFirst;
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
            KRATOS_ERROR << "There is error in computing the bisection. ProjectOnCurveUsingBisection shall detect it at the beginning of this function.";
        }
        else if ((test1 > 0.0) &&  (test2 > 0.0)) // the point is after tmax
        {
            // project on to the line created by Last and dLast
            double a = inner_prod(P - Last, dLast) / inner_prod(dLast, dLast);
            noalias(Proj) = Last + a * dLast;
        }
        else if ((test1 < 0.0) &&  (test2 > 0.0)) // this can happen if the curve has U-shape
        {
            // project to closer line
            double a = inner_prod(P - First, dFirst) / inner_prod(dFirst, dFirst);
            noalias(Proj) = First + a * dFirst;
            double d1 = norm_2(P - Proj);

            double b = inner_prod(P - Last, dLast) / inner_prod(dLast, dLast);
            noalias(Proj) = Last + b * dLast;
            double d2 = norm_2(P - Proj);

            if (d1 < d2)
            {
                noalias(Proj) = First + a * dFirst;
            }
            else
            {
                noalias(Proj) = Last + b * dLast;
            }
        }
    }

    return norm_2(P - Proj);
}

int Curve::ProjectOnCurveUsingBisection(const PointType& P, PointType& Proj, double& t,
            const double tmin, const double tmax, const int nsampling, const double tol) const
{
    // firstly do the sampling
    std::vector<double> f(nsampling);
    PointType dProj;
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        t = tmin + i * (tmax - tmin) / nsampling;
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
            t = tmin + i * (tmax - tmin) / nsampling;
            noalias(Proj) = this->GetValue(t);
            break;
        }

        if (f[i]*f[i + 1] < 0.0)
        {
            // found the segment, do the bisection
            double left = tmin + i * (tmax - tmin) / nsampling;
            double right = tmin + (i + 1) * (tmax - tmin) / nsampling;
            double mid;
            double fleft = f[i], fright = f[i + 1], fmid;
            while ((right - left) > tol)
            {
                mid = 0.5 * (left + right);

                noalias(Proj) = this->GetValue(mid);
                noalias(dProj) = this->GetDerivative(0, mid);
                fmid = inner_prod(dProj, P - Proj);

                if (fabs(fmid) < tol)
                {
                    break;
                }

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
        if (f[nsampling - 1] > tol)
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
    //     KRATOS_ERROR << "Bisection error: there are no valid segment";
    // }

    return 0;
}

int Curve::ComputeIntersectionByBisection(const LevelSet& rLevelSet, PointType& IntPoint, double& t,
            const double tmin, const double tmax, const std::size_t nsampling, const double tol) const
{
    // firstly do the sampling
    std::vector<double> f(nsampling);
    for (std::size_t i = 0; i < nsampling + 1; ++i)
    {
        t = tmin + i * (tmax - tmin) / nsampling;
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
            t = tmin + i * (tmax - tmin) / nsampling;
            noalias(IntPoint) = this->GetValue(t);
            break;
        }

        if (f[i]*f[i + 1] < 0.0)
        {
            // found the segment, do the bisection
            double left = tmin + i * (tmax - tmin) / nsampling;
            double right = tmin + (i + 1) * (tmax - tmin) / nsampling;
            double mid;
            double fleft = f[i], fright = f[i + 1], fmid;
            while ((right - left) > tol)
            {
                mid = 0.5 * (left + right);

                noalias(IntPoint) = this->GetValue(mid);
                fmid = rLevelSet.GetValue(IntPoint);

                if (fabs(fmid) < tol)
                {
                    break;
                }

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
        if (f[nsampling - 1] > tol)
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
    //     KRATOS_ERROR << "Bisection error: there are no valid segment";
    // }

    return 0;
}

void Curve::ComputeEquallyDistanceDivisionByRecursiveBisection(const double tmin, const double tmax,
                const unsigned int num_division, std::vector<double>& tvec, const double tolerance) const
{
    if (num_division < 2)
        KRATOS_ERROR << "num_division must be > 2. The value is " << num_division;

    if (num_division == 2)
    {
        if (tvec.size() != 3)
            tvec.resize(3);
        tvec[0] = tmin;
        tvec[2] = tmax;
        this->ProjectMiddleUsingBisection(tmin, tmax, tvec[1], tolerance);
    }
    else
    {
        double tleft = tmin, tright = tmax;

        const auto Pmin = this->GetValue(tmin);
        const auto Pmax = this->GetValue(tmax);

        while (true)
        {
            const double tmid = 0.5*(tleft + tright);
            const auto Pmid = this->GetValue(tmid);

            const double d1 = norm_2(Pmin - Pmid);

            std::vector<double> tvecr;
            this->ComputeEquallyDistanceDivisionByRecursiveBisection(tmid, tmax, num_division-1, tvecr, tolerance);
            const auto Pr1 = this->GetValue(tvecr[0]);
            const auto Pr2 = this->GetValue(tvecr[1]);
            const double d2 = norm_2(Pr1 - Pr2);

            if (std::abs(d1 - d2) < tolerance)
            {
                if (tvec.size() != num_division+1)
                    tvec.resize(num_division+1);

                tvec[0] = tmin;
                for (std::size_t i = 0; i < num_division; ++i)
                    tvec[i+1] = tvecr[i];

                break;
            }

            if (d1 > d2)
            {
                tright = tmid;
            }
            else
            {
                tleft = tmid;
            }
        }
    }
}

void Curve::ProjectMiddleUsingBisection(const double tmin, const double tmax, double& tmid, const double tolerance) const
{
    double tleft = tmin, tright = tmax;

    const auto Pmin = this->GetValue(tmin);
    const auto Pmax = this->GetValue(tmax);

    while (true)
    {
        tmid = 0.5*(tleft + tright);
        const auto Pmid = this->GetValue(tmid);

        const double d1 = norm_2(Pmin - Pmid);
        const double d2 = norm_2(Pmax - Pmid);

        if (std::abs(d1 - d2) < tolerance)
            return;

        if (d1 > d2)
        {
            tright = tmid;
        }
        else
        {
            tleft = tmid;
        }
    }
}

void Curve::ComputeUniformDivision(const double tmin, const double tmax,
            const unsigned int num_division, std::vector<double>& tvec, const double tolerance) const
{
    // compute the length
    const double length = this->ComputeLength(tmin, tmax);

    // compute the estimated arc length
    const double arc_length = length / num_division;

    // find the uniform distributed point sequentially
    if (tvec.size() != num_division + 1)
        tvec.resize(num_division + 1);

    tvec[0] = tmin;
    for (unsigned int i = 0; i < num_division-1; ++i)
    {
        this->FindPointOnCurveWithArcLengthUsingBisection(tvec[i], tmax, tvec[i+1], arc_length, tolerance);
    }
    tvec[num_division] = tmax;
}

void Curve::FindPointOnCurveWithArcLengthUsingBisection(const double tmin, const double tmax, double& tmid,
            const double arc_length, const double tolerance) const
{
    double tleft = tmin, tright = tmax;

    while (true)
    {
        tmid = 0.5*(tleft + tright);

        const double length = this->ComputeLength(tmin, tmid);

        if (std::abs(length - arc_length) < tolerance)
            return;

        if (length > arc_length)
        {
            tright = tmid;
        }
        else
        {
            tleft = tmid;
        }
    }
}

} // end namespace Kratos
