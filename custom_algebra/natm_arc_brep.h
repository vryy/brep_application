//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         NATMArcBRep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            22 Jun 2020
//


#if !defined(KRATOS_NATM_ARC_BREP_H_INCLUDED )
#define  KRATOS_NATM_ARC_BREP_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/brep_math_utility.h"


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

struct NATMArc
{
    double X, Y, R;
    double Alpha1, Alpha2;
    double Angle1, Angle2, Angle3;
};

inline std::ostream& operator << (std::ostream& rOStream, const NATMArc& rThis)
{
    rOStream << "(";
    rOStream << "X: " << rThis.X;
    rOStream << ", Y: " << rThis.Y;
    rOStream << ", R: " << rThis.R;
    rOStream << ", Alpha1: " << rThis.Alpha1;
    rOStream << ", Alpha2: " << rThis.Alpha2;
    rOStream << ", Angle1: " << rThis.Angle1;
    rOStream << ", Angle2: " << rThis.Angle2;
    rOStream << ", Angle3: " << rThis.Angle3;
    rOStream << ")";
    return rOStream;
}

/// Short class definition.
/**
 * BRep represented for NATM tunnel constructed from piecewise arc curves
 */
class NATMArcBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NATMArcBRep
    KRATOS_CLASS_POINTER_DEFINITION(NATMArcBRep);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::PointType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NATMArcBRep() : BaseType()
    {}

    /// Copy constructor.
    NATMArcBRep(NATMArcBRep const& rOther) : BaseType(rOther)
    , mArcs(rOther.mArcs), mRefAngle(rOther.mRefAngle), mRefCenter(rOther.mRefCenter)
    {}

    /// Destructor.
    virtual ~NATMArcBRep()
    {}


    ///@}
    ///@name Self Operations
    ///@{

    ///
    void SetReferenceCenter(const double& X, const double& Y)
    {
        mRefCenter[0] = X;
        mRefCenter[1] = Y;
        mRefCenter[2] = 0.0;
    }

    /// The reference point shall be the first point of the profile
    void SetReferencePoint(const double& X, const double& Y)
    {
        mRefAngle = std::atan2(Y - mRefCenter[1], X - mRefCenter[0]);
        KRATOS_WATCH(mRefAngle)
    }

    /// Add the arc to the BRep. The arc must be in counter-clockwise direction.
    void AddArc(const double& X1, const double& Y1,
                const double& X2, const double& Y2,
                const double& X3, const double& Y3)
    {
        double x1 = X2 - X1, y1 = Y2 - Y1;
        double x2 = X3 - X1, y2 = Y3 - Y1;
        double xc, yc;

        BRepMathUtility::Solve(2*x1, 2*y1, x1*x1 + y1*y1,
                2*x2, 2*y2, x2*x2 + y2*y2, xc, yc);

        NATMArc a;
        a.X = X1 + xc;
        a.Y = Y1 + yc;
        a.R = sqrt(pow(x1 - xc, 2) + pow(y1 - yc, 2));

        const double pi = std::atan(1)*4;
        // KRATOS_WATCH(pi)

        a.Alpha1 = std::atan2(Y1 - mRefCenter[1], X1 - mRefCenter[0]) - mRefAngle;
        // KRATOS_WATCH(a.Alpha1)
        if (a.Alpha1 < 0.0) a.Alpha1 += 2*pi;
        a.Alpha2 = std::atan2(Y3 - mRefCenter[1], X3 - mRefCenter[0]) - mRefAngle;
        // KRATOS_WATCH(a.Alpha2)
        if (a.Alpha2 <= 0.0) a.Alpha2 += 2*pi;

        Vector v1(2); v1[0] = X1 - a.X; v1[1] = Y1 - a.Y;
        Vector v2(2); v2[0] = X2 - a.X; v2[1] = Y2 - a.Y;
        Vector v3(2); v3[0] = X3 - a.X; v3[1] = Y3 - a.Y;
        a.Angle1 = std::atan2(v1[1], v1[0]);
        if (a.Angle1 < 0.0) a.Angle1 += 2*pi;
        a.Angle2 = a.Angle1 + std::acos( inner_prod(v1, v2) / (norm_2(v1) * norm_2(v2)) );
        a.Angle3 = a.Angle2 + std::acos( inner_prod(v3, v2) / (norm_2(v3) * norm_2(v2)) );

        // KRATOS_WATCH(a)

        mArcs.push_back(a);
    }

    ///@}
    ///@name Derived Operations
    ///@{

    /// Clone this BRep
    virtual BRep::Pointer CloneBRep() const
    {
        return BRep::Pointer(new NATMArcBRep(*this));
    }

    /// Get working space dimension
    virtual std::size_t WorkingSpaceDimension() const
    {
        return 2;
    }

    /// Get local space dimension
    virtual std::size_t LocalSpaceDimension() const
    {
        return 1;
    }

    /// Check if a point is inside/outside of the NATMArcBRep
    virtual bool IsInside(const PointType& P) const
    {
        // check which arc the point is in
        const double pi = std::atan(1)*4;
        double alpha = std::atan2(P[1] - mRefCenter[1], P[0] - mRefCenter[0]) - mRefAngle;
        if (alpha <= 0.0) alpha += 2*pi;

        bool is_in = false;
        for (std::size_t i = 0; i < mArcs.size(); ++i)
        {
            if (mArcs[i].Alpha2 > mArcs[i].Alpha1)
            {
                if (mArcs[i].Alpha2 >= alpha && alpha >= mArcs[i].Alpha1)
                    is_in = true;
            }

            if (mArcs[i].Alpha1 > mArcs[i].Alpha2)
            {
                if (mArcs[i].Alpha1 >= alpha && alpha >= mArcs[i].Alpha2)
                    is_in = true;
            }

            if (is_in)
            {
                // KRATOS_WATCH(P)
                // KRATOS_WATCH(i)
                // KRATOS_WATCH(mArcs[i])
                // KRATOS_WATCH(alpha)

                // compute the intersection point
                double a = P[1] - mRefCenter[1];
                double b = mRefCenter[0] - P[0];
                double c = -a*P[0] - b*P[1];
                std::vector<double> points = BRepMathUtility::Intersect(a, b, c, mArcs[i].X, mArcs[i].Y, mArcs[i].R);

                // KRATOS_WATCH(points.size())

                if (points.size() == 0 || points.size() == 2)
                {
                    KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the profile. The line from center to point must cut the arc in 2 points", "")
                }

                // KRATOS_WATCH(points[0])
                // KRATOS_WATCH(points[1])
                // KRATOS_WATCH(points[2])
                // KRATOS_WATCH(points[3])

                // find the right intersection point
                double angle1 = std::atan2(points[1] - mArcs[i].Y, points[0] - mArcs[i].X);
                if (angle1 < 0.0) angle1 += 2*pi;
                double angle2 = std::atan2(points[3] - mArcs[i].Y, points[2] - mArcs[i].X);
                if (angle2 < 0.0) angle2 += 2*pi;

                // KRATOS_WATCH(angle1)
                // KRATOS_WATCH(angle2)

                double xi, yi;
                if ( (angle1 >= mArcs[i].Angle1 && angle1 <= mArcs[i].Angle3)
                  || (angle1+2*pi >= mArcs[i].Angle1 && angle1+2*pi <= mArcs[i].Angle3) )
                {
                    xi = points[0];
                    yi = points[1];
                }
                else if ( (angle2 >= mArcs[i].Angle1 && angle2 <= mArcs[i].Angle3)
                       || (angle2+2*pi >= mArcs[i].Angle1 && angle2+2*pi <= mArcs[i].Angle3) )
                {
                    xi = points[2];
                    yi = points[3];
                }
                else
                {
                    // KRATOS_WATCH("outside")
                    return false;
                }

                double d = sqrt(pow(P[0] - mRefCenter[0], 2) + pow(P[1] - mRefCenter[1], 2));
                double di = sqrt(pow(xi - mRefCenter[0], 2) + pow(yi - mRefCenter[1], 2));

                // KRATOS_WATCH(xi)
                // KRATOS_WATCH(yi)

                // KRATOS_WATCH(d)
                // KRATOS_WATCH(di)

                if (d < di)
                {
                    // KRATOS_WATCH("inside")
                    return true;
                }
                else
                {
                    // KRATOS_WATCH("outside")
                    return false;
                }
            }
        }

        KRATOS_WATCH(P)
        KRATOS_WATCH(alpha)
        KRATOS_THROW_ERROR(std::logic_error, "That is not possible, the point does not lie within any domain covered by the arc", "")

        return false;
    }

    /// Check if a point is on the boundary within a tolerance
    virtual bool IsOnBoundary(const PointType& P, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
    }

    /// Compute the intersection of the NATMArcBRep with a line connect by 2 points.
    virtual PointType Bisect(const PointType& P1, const PointType& P2, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Not yet implemented", "")
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
    virtual std::string Info() const
    {
        return "NATM Arc BRep";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    double mRefAngle;
    array_1d<double, 3> mRefCenter;
    std::vector<NATMArc> mArcs;

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
    NATMArcBRep& operator=(NATMArcBRep const& rOther);

    ///@}

}; // Class NATMArcBRep

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NATMArcBRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NATMArcBRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NATM_ARC_BREP_H_INCLUDED  defined
