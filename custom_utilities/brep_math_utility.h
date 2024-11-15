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
#include "geometries/geometry_data.h"

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
template<typename TDataType = double>
class BRepMathUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRepMathUtility
    KRATOS_CLASS_POINTER_DEFINITION(BRepMathUtility);

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

    static inline constexpr int GetMaxIntegrationOrder()
    {
        return 5;
    }

    static inline GeometryData::IntegrationMethod GetIntegrationMethod(const int integration_order)
    {
        if (integration_order > GetMaxIntegrationOrder())
        {
            KRATOS_ERROR << "Does not support for integration rule with order > " << GetMaxIntegrationOrder();
        }

        GeometryData::IntegrationMethod ThisIntegrationMethod;

        if (integration_order == 1)
        {
            ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
        }
        else if (integration_order == 2)
        {
            ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
        else if (integration_order == 3)
        {
            ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
        }
        else if (integration_order == 4)
        {
            ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
        }
        else if (integration_order == 5)
        {
            ThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
        }
        else
        {
            KRATOS_ERROR << "Unsupported integration order " << integration_order;
        }

        return ThisIntegrationMethod;
    }

    /**
     * Solve a1 * x + b1 * y = c1
     *       a2 * x + b2 * y = c2
     */
    static void Solve(const TDataType a1, const TDataType b1, const TDataType c1,
                      const TDataType a2, const TDataType b2, const TDataType c2, TDataType& x, TDataType& y)
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
    static std::vector<TDataType> SolveQuadratic(const TDataType a, const TDataType b, const TDataType c)
    {
        TDataType b2 = b / 2;
        TDataType delta = b2 * b2 - a * c;
        std::vector<TDataType> x;

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
    static std::vector<TDataType> Intersect(const TDataType a, const TDataType b, const TDataType c,
                                         const TDataType xc, const TDataType yc, const TDataType R)
    {
        std::vector<TDataType> Output;
        if (a == 0.0)
        {
            TDataType y = -c / b;

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
            TDataType x = -c / a;

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
            TDataType A = 1.0 + pow(a / b, 2);
            TDataType B = 2.0 * (-xc + a / b * (c / b + yc));
            TDataType C = pow(xc, 2) + pow(c / b + yc, 2) - pow(R, 2);

            std::vector<TDataType> x = SolveQuadratic(A, B, C);
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
template<typename TDataType>
inline std::istream& operator >> (std::istream& rIStream, BRepMathUtility<TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template<typename TDataType>
inline std::ostream& operator << (std::ostream& rOStream, const BRepMathUtility<TDataType>& rThis)
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
