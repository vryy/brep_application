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
//  Date:            20 Mar 2021
//


#if !defined(KRATOS_BREP_INTERSECTION_UTILITY_H_INCLUDED )
#define  KRATOS_BREP_INTERSECTION_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/brep.h"
#include "custom_algebra/level_set/distance_to_curve_level_set.h"
#include "custom_algebra/level_set/planar_level_set.h"


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
/** class to compute the intersection between BRep, Level Set, etc
*/
class BRepIntersectionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRepIntersectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(BRepIntersectionUtility);

    typedef BRep::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRepIntersectionUtility() {}

    /// Destructor.
    virtual ~BRepIntersectionUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Sampling the intersection of two BReps
    static void Intersect(std::vector<PointType>& rPoints, const BRep& rBRep1, const BRep& rBRep2, const std::size_t& nsampling)
    {
        if (dynamic_cast<const DistanceToCurveLevelSet*>(&rBRep1) && dynamic_cast<const PlanarLevelSet*>(&rBRep2))
        {
            const DistanceToCurveLevelSet& rLS1 = dynamic_cast<const DistanceToCurveLevelSet&>(rBRep1);
            const PlanarLevelSet& rLS2 = dynamic_cast<const PlanarLevelSet&>(rBRep2);
            Intersect(rPoints, rLS1, rLS2, nsampling);
        }
        else
        {
            std::stringstream ss;
            ss << "Intersection between " << rBRep1.Info() << " and " << rBRep2.Info() << " is not implemented";
            KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
        }
    }

    /// Sampling the intersection between DistanceToCurveLevelSet and PlanarLevelSet
    static void Intersect(std::vector<PointType>& rPoints, const DistanceToCurveLevelSet& rLS1, const PlanarLevelSet& rLS2, const std::size_t& nsampling)
    {
        // compute the intersection of the alignment curve and the plane
        PointType C, T;
        int error_code = rLS1.AlignmentCurve().ComputeIntersection(rLS2, C, T);
        if (error_code != 0)
            KRATOS_THROW_ERROR(std::logic_error, "The alignment curve and the plane do not intersect", "")
        // KRATOS_WATCH(C)
        // KRATOS_WATCH(T)

        // compute the angle between the plane and the axial (tangent) vector
        PointType N;
        rLS2.NormalVector(N);
        // KRATOS_WATCH(N)

        double sina = inner_prod(N, T) / (norm_2(N) * norm_2(T));

        // compute the effective projected radius
        double eradius = rLS1.Radius() / sina;
        // KRATOS_WATCH(eradius)

        // compute a local Frenet frame at the intersection point
        PointType T1, T2;
        T2[0] = 0.0; T2[1] = 0.0; T2[2] = 1.0;
        noalias(T1) = MathUtils<double>::CrossProduct(N, T2);
        if (norm_2(T1) < 1.0e-10)
        {
            T2[0] = 0.0; T2[1] = 1.0; T2[2] = 0.0;
            noalias(T1) = MathUtils<double>::CrossProduct(N, T2);
        }
        T1 /= norm_2(T1);
        noalias(T2) = MathUtils<double>::CrossProduct(N, T1);
        T2 /= norm_2(T2);
        // KRATOS_WATCH(T1)
        // KRATOS_WATCH(T2)

        // perform the sampling on the effective circle
        double a, c, s;
        constexpr double pi = atan(1.0)*4;
        constexpr double rfactor = 2.0; // we take the factor 2 to make sure the point is outside the distance-to-curve level set
        PointType P;
        rPoints.resize(nsampling);
        for (std::size_t i = 0; i < nsampling; ++i)
        {
            a = 2*pi*i/nsampling;
            c = cos(a);
            s = sin(a);

            noalias(P) = C + rfactor*eradius*(c*T1 + s*T2);

            // compute the intersection between CP and the distance-to-curve level set
            error_code = rLS1.Bisect(rPoints[i], C, P, rLS1.Tolerance());
            if (error_code != 0)
                std::cout << "Error computing sampling points on intersection plane at i = " << i
                          << ", P = " << P
                          << std::endl;
        }
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
        return "BRep Intersection Utility";
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
    BRepIntersectionUtility& operator=(BRepIntersectionUtility const& rOther);

    /// Copy constructor.
    BRepIntersectionUtility(BRepIntersectionUtility const& rOther);

    ///@}

}; // Class BRepIntersectionUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRepIntersectionUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRepIntersectionUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_BREP_INTERSECTION_UTILITY_H_INCLUDED  defined
