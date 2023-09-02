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

#if !defined(KRATOS_LEVEL_SET_H_INCLUDED )
#define  KRATOS_LEVEL_SET_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "geometries/geometry_data.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/brep.h"

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
/** Abstract class for a level set in space, both for implicit level set or nodal interpolated level set
 * A level set is characterized by mapping phi: R^3 -> R, and partitions the space to 3 domains by criteria:
 *   phi(x) < 0  <=>   x is inside of \Omega
 *   phi(x) = 0  <=>   x is on \Gamma
 *   phi(x) > 0  <=>   x is outside of \Gamma
 * REF: Massing et al, CutFEM: Discretizing geometry and partial differential equations
 */
class LevelSet : public FunctionR3R1, public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LevelSet
    KRATOS_CLASS_POINTER_DEFINITION(LevelSet);

    typedef FunctionR3R1 BaseType;

    typedef BRep::GeometryType GeometryType;

    typedef BRep::NodeType NodeType;

    typedef BRep::PointType PointType;

    typedef BRep::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LevelSet() : BRep(), BaseType() {}

    /// Copy constructor.
    LevelSet(LevelSet const& rOther) : BRep(rOther), BaseType(rOther) {}

    /// Destructor.
    virtual ~LevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// inherit from BRep
    BRep::Pointer CloneBRep() const override
    {
        return BRep::Pointer(this->CloneLevelSet());
    }

    /// inherit from Function
    FunctionR3R1::Pointer CloneFunction() const override
    {
        return FunctionR3R1::Pointer(this->CloneLevelSet());
    }

    /// Clone this level set
    virtual LevelSet::Pointer CloneLevelSet() const
    {
        return LevelSet::Pointer(new LevelSet());
    }

    ///@}
    ///@name Access
    ///@{

    /// inherit from BRep
    std::size_t WorkingSpaceDimension() const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// inherit from BRep
    std::size_t LocalSpaceDimension() const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// Get level set value at a point
    double GetValue(const double& X, const double& Y, const double& Z) const
    {
        PointType P;
        P[0] = X;
        P[1] = Y;
        P[2] = Z;
        return this->GetValue(P);
    }

    /// inherit from Function
    double GetValue(const array_1d<double, 3>& P) const final
    {
        PointType PP(P);
        return this->GetValue(PP);
    }

    /// Get level set value at a point
    virtual double GetValue(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

//    virtual double GetValue(GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
//    }

    /// Compute the gradient at a point, i.e. [d Phi / dX, d Phi / dY, d Phi / dZ]
    Vector GetGradient(const double& X, const double& Y, const double& Z) const
    {
        PointType P;
        P[0] = X;
        P[1] = Y;
        P[2] = Z;
        return this->GetGradient(P);
    }

    /// inherit from Function
    Vector GetGradient(const array_1d<double, 3>& P) const final
    {
        PointType PP(P);
        return this->GetGradient(PP);
    }

    /// Get the gradient at point
    virtual Vector GetGradient(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// compute the derivatives of the gradient w.r.t the global point
    /// [d^2 Phi / d X^2, d^2 Phi / dX dY, d^2 Phi / dX dZ
    ///  d^2 Phi / dX dY, d^2 Phi / d^2 Y, d^2 Phi / dY dZ
    ///  d^2 Phi / dX dZ, d^2 Phi / dY dZ, d^2 Phi / d Z^2]
    virtual Matrix GetGradientDerivatives(const PointType& P) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

//    virtual Vector GetGradient(GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
//    {
//        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
//    }

    /// inherit from BRep
    bool IsInside(const PointType& P) const override
    {
        return (this->GetValue(P) < 0.0);
    }

    /// inherit from BRep
    bool IsOnBoundary(const PointType& P) const override
    {
        return (fabs(this->GetValue(P)) < this->Tolerance());
    }

    /// inherit from BRep
    /// Check if a geometry is cut by the level set
    int CutStatus(GeometryType& r_geom, const int& configuration) const override
    {
        if (configuration == 0)
        {
            std::vector<PointType> points(r_geom.size());
            for (std::size_t i = 0; i < r_geom.size(); ++i)
            {
                noalias(points[i]) = r_geom[i].GetInitialPosition();
            }
            return this->CutStatusOfPoints(points, this->GetTolerance());
        }
        else if (configuration == 1)
        {
            return this->CutStatusOfPoints(r_geom, this->GetTolerance());
            // REMARK: this will use the current position of node, e.g. in dynamics
        }
        else
        {
            return BRep::_INVALID;
        }
    }

    /// inherit from BRep
    /// Check if a set of points is cut by the level set
    int CutStatus(const std::vector<PointType>& r_points) const override
    {
        return this->CutStatusOfPoints(r_points, this->GetTolerance());
    }

    /// inherit from BRep
    /// Check if a geometry is cut by the BRep by sampling the geometry
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    int CutStatusBySampling(GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const override
    {
        return BRep::CutStatusBySampling(r_geom, nsampling, configuration);
    }

    /// inherit from BRep
    /// Compute the intersection of the level set with a line connect by 2 points.
    /// Note that, the checking of the intersection of the level set with the line is not performed. Hence one should ensure that before calling this function.
    int Bisect(PointType& P, const PointType& P1, const PointType& P2, const double& Tol) const override
    {
        double f1 = this->GetValue(P1);
        double f2 = this->GetValue(P2);
        if (f1 * f2 > 0.0)
        {
            // KRATOS_THROW_ERROR(std::logic_error, "Bisect does not work with two ends at the same side", "")
            return -1;
        }

        double left = 0.0;
        double right = 1.0;

        bool converged = false;
        while (!converged)
        {
            double mid = (left + right) / 2;
            noalias(P) = (1.0 - mid) * P1 + mid * P2;
            double fm = this->GetValue(P);

            if (fabs(fm) < Tol)
            {
                converged = true;
            }
            else
            {
                if (fm * f1 < 0.0)
                {
                    right = mid;
                    f2 = fm;
                }
                else
                {
                    left = mid;
                    f1 = fm;
                }

                // if(right-left < Tol)
                if (fabs(f2 - f1) < Tol)
                {
                    converged = true;
                }
            }
        }

        return 0;
    }

    /// inherit from BRep
    /// Compute the intersection of the level set with a line connect by 2 points.
    /// Note that, the checking of the intersection of the level set with the line is not performed. Hence one should ensure that before calling this function.
    int Bisect(PointType& P, const PointType& P1, const PointType& P2, const PointType& P3, const double& Tol) const override
    {
        double f1 = this->GetValue(P1);
        double f2 = this->GetValue(P2);
        if (f1 * f2 > 0.0)
        {
            // KRATOS_THROW_ERROR(std::logic_error, "Bisect does not work with two ends at the same side", "")
            return -1;
        }

        double left = -1.0;
        double right = 1.0;

        bool converged = false;
        while (!converged)
        {
            double mid = (left + right) / 2;
            noalias(P) = P1 * L3_N1(mid) + P2 * L3_N2(mid) + P3 * L3_N3(mid);
            double fm = this->GetValue(P);

            if (fabs(fm) < Tol)
            {
                converged = true;
            }
            else
            {
                if (fm * f1 < 0.0)
                {
                    right = mid;
                    f2 = fm;
                }
                else
                {
                    left = mid;
                    f1 = fm;
                }

                // if(right-left < Tol)
                if (fabs(f2 - f1) < Tol)
                {
                    converged = true;
                }
            }
        }

        return 0;
    }

    /// inherit from BRep
    void GetNormal(const PointType& P, PointType& rNormal) const override
    {
        Vector G = this->GetGradient(P);
        noalias(rNormal) = ZeroVector(3);
        for (std::size_t i = 0; i < G.size(); ++i)
        {
            rNormal[i] = G(i);
        }
    }

    /// inherit from BRep
    void GetNormalDerivatives(const PointType& P, Matrix& Derivatives) const override
    {
        Derivatives = this->GetGradientDerivatives(P);
    }

    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// compute the derivatives of the projection point w.r.t to the original point.
    /// The derivatives are organized as;
    ///     [d Proj[0] / d P[0], d Proj[0] / d P[1], d Proj[0] / d P[2]]
    ///     [d Proj[1] / d P[0], d Proj[1] / d P[1], d Proj[1] / d P[2]]
    ///     [d Proj[2] / d P[0], d Proj[2] / d P[1], d Proj[2] / d P[2]]
    void ProjectionDerivatives(const PointType& P, Matrix& Derivatives) const override
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Level Set";
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

    template<class TPointsContainerType>
    int CutStatusOfPoints(const TPointsContainerType& r_points, const double& tolerance) const
    {
        std::vector<std::size_t> in_list, out_list, on_list;
        for (std::size_t v = 0; v < r_points.size(); ++v)
        {
            double phi = this->GetValue(r_points[v]);
            if (phi < -tolerance)
            {
                in_list.push_back(v);
            }
            else if (phi > tolerance)
            {
                out_list.push_back(v);
            }
            else
            {
                on_list.push_back(v);
            }
        }

        int stat;
        if (in_list.size() == 0 && out_list.size() == 0)
        {
            if (on_list.size() > 0)
            {
                return BRep::_ON; // every points are on
            }
            else
            {
                std::cout << "----------------------------------------" << std::endl;
                for (std::size_t v = 0; v < r_points.size(); ++v)
                {
                    KRATOS_WATCH(r_points[v])
                }
                KRATOS_WATCH(in_list.size())
                KRATOS_WATCH(out_list.size())
                KRATOS_WATCH(on_list.size())
                KRATOS_WATCH(tolerance)
                KRATOS_THROW_ERROR(std::logic_error, "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.", "")
                // std::cout << "!!!WARNING!!!The geometry is degenerated." << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                return BRep::_INVALID; // geometry is degenerated
            }
        }
        else
        {
            if (on_list.size() > 0)
            {
                stat = BRep::_CUT;
                return stat;
            }

            if (in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if (out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }

    /// Shape function of 3-node line
    static inline double L3_N1(const double& xi) {return 0.5 * (xi - 1.0) * xi;}
    static inline double L3_N2(const double& xi) {return 0.5 * (xi + 1.0) * xi;}
    static inline double L3_N3(const double& xi) {return 1.0 - xi * xi;}

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
    LevelSet& operator=(LevelSet const& rOther);

    ///@}

}; // Class LevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, LevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const LevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LEVEL_SET_H_INCLUDED  defined
