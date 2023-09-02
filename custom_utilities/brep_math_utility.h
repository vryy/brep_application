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
//  Date:            22 Jun 2020
//

#if !defined(KRATOS_BREP_MATH_UTILITY_H_INCLUDED )
#define  KRATOS_BREP_MATH_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

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
/** class for auxilliary routines
*/
class BRepMathUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRepMathUtility
    KRATOS_CLASS_POINTER_DEFINITION(BRepMathUtility);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRepMathUtility() {}

    /// Destructor.
    virtual ~BRepMathUtility() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Solve a1 * x + b1 * y = c1
     *       a2 * x + b2 * y = c2
     */
    static void Solve(const double& a1, const double& b1, const double& c1,
                      const double& a2, const double& b2, const double& c2, double& x, double& y)
    {
        x = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
        y = (c2 * a1 - c1 * a2) / (a1 * b2 - a2 * b1);
    }

    /*
     * Solve the quadratic equation a*x^2 + b*x + c = 0
     * In return:   0: no solution
     *              1: one solution
     *              2: two solutions
     */
    static std::vector<double> SolveQuadratic(const double& a, const double& b, const double& c)
    {
        double b2 = b / 2;
        double delta = b2 * b2 - a * c;
        std::vector<double> x;

        if (delta < 0.0)
        {
            return x; // no solution
        }
        else if (delta == 0.00)
        {
            if (a == 0.0)
            {
                return x; // no solution
            }
            else
            {
                x.push_back(-b2 / a); // one solution
                return x;
            }
        }
        else
        {
            if (a == 0.0)
            {
                x.push_back(-c / b); // one solution
                return x;
            }
            else
            {
                x.push_back((-b2 - sqrt(delta)) / a);
                x.push_back((-b2 + sqrt(delta)) / a);
                // two solution
                return x;
            }
        }
    }

    /**
     * Compute the intersection between line ax+by+c=0 and circle (x-xc)^2 + (y-yc)^2 = R^2
     * Depending on the cut situation, the output may have 0, 2 or 4 results
     * The input shall exclude the special sitation, e.g. a=b=0 or R=0
     */
    static std::vector<double> Intersect(const double& a, const double& b, const double& c,
                                         const double& xc, const double& yc, const double& R)
    {
        std::vector<double> Output;
        if (a == 0.0)
        {
            double y = -c / b;

            if (R >= abs(y - yc))
            {
                Output.push_back(xc - sqrt(pow(R, 2) - pow(y - yc, 2)));
                Output.push_back(y);
                Output.push_back(xc + sqrt(pow(R, 2) - pow(y - yc, 2)));
                Output.push_back(y);
            }
        }
        else if (b == 0.0)
        {
            double x = -c / a;

            if (R >= abs(x - xc))
            {
                Output.push_back(x);
                Output.push_back(yc - sqrt(pow(R, 2) - pow(x - xc, 2)));
                Output.push_back(x);
                Output.push_back(yc + sqrt(pow(R, 2) - pow(x - xc, 2)));
            }
        }
        else
        {
            double A = 1.0 + pow(a / b, 2);
            double B = 2.0 * (-xc + a / b * (c / b + yc));
            double C = pow(xc, 2) + pow(c / b + yc, 2) - pow(R, 2);

            std::vector<double> x = SolveQuadratic(A, B, C);
            if (x.size() == 1)
            {
                Output.push_back(x[0]);
                Output.push_back(-(a * x[0] + c) / b);
            }
            else if (x.size() == 2)
            {
                Output.push_back(x[0]);
                Output.push_back(-(a * x[0] + c) / b);
                Output.push_back(x[1]);
                Output.push_back(-(a * x[1] + c) / b);
            }
        }

        return Output;
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
        return "BRep Utility";
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
    BRepMathUtility& operator=(BRepMathUtility const& rOther);

    /// Copy constructor.
    BRepMathUtility(BRepMathUtility const& rOther);

    ///@}

}; // Class BRepMathUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRepMathUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRepMathUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BREP_MATH_UTILITY_H_INCLUDED  defined
