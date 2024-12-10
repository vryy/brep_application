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
//  Date:            13 Mar 2017
//

#if !defined(KRATOS_BREP_H_INCLUDED )
#define  KRATOS_BREP_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "section/section.h"

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
/** Abstract class for a boundary representation
*/
class KRATOS_API(BREP_APPLICATION) BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BRep
    KRATOS_CLASS_POINTER_DEFINITION(BRep);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    static const int _CUT = -1;
    static const int _IN  = 0;
    static const int _OUT = 1;
    static const int _ON = -2;
    static const int _INVALID = -3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BRep();

    /// Copy constructor.
    BRep(BRep const& rOther) : mTOL(rOther.mTOL), mName(rOther.mName), mEchoLevel(0)
    {}

    /// Destructor.
    virtual ~BRep();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Clone this BRep
    virtual BRep::Pointer CloneBRep() const
    {
        return BRep::Pointer(new BRep(*this));
    }

    /// Expose this BRep as a shared_ptr. Host BRep can also use this to expose the internal BRep.
    virtual BRep::ConstPointer pBRep() const
    {
        return BRep::ConstPointer(this);
    }

    /// Set name for this BRep
    void SetName(const std::string& Name) {mName = Name;}

    /// Get name of this BRep
    std::string GetName() const {return mName;}
    std::string Name() const {return mName;}

    /// Set for geometric tolerance
    void SetTolerance(const double& TOL) {mTOL = TOL;}

    /// Get for geometric tolerance
    double GetTolerance() const {return mTOL;}
    double Tolerance() const {return mTOL;}

    /// Set the echo level
    void SetEchoLevel(const int echo_level) {mEchoLevel = echo_level;}

    /// Get the echo level
    int GetEchoLevel() const {return mEchoLevel;}

    /// Get working space dimension
    virtual std::size_t WorkingSpaceDimension() const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Get local space dimension
    virtual std::size_t LocalSpaceDimension() const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Check if a point is inside/outside of the BRep
    virtual bool IsInside(const PointType& P) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in reference configuration
    /// Since now C++ does not support virtual template function, this function must be separated to 2 functions
    bool IsInside(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords, const int& configuration) const
    {
        if (configuration == 0)
        {
            return IsInside0(rGeometry, local_coords);
        }
        else if (configuration == 1)
        {
            return IsInside1(rGeometry, local_coords);
        }
        else
        {
            KRATOS_ERROR << "Invalid configuration " << configuration;
        }
    }

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in reference configuration
    /// Since now C++ does not support virtual template function, this function must be separated to 2 functions
    virtual bool IsInside0(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const;

    /// Check if a point is inside/outside of the BRep
    /// The point will be interpolated in current configuration
    virtual bool IsInside1(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const;

    /// Check if a point is on the boundary
    virtual bool IsOnBoundary(const PointType& P) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Check if an element is cut by the brep
    int CutStatus(Element::Pointer p_elem, const int& configuration) const;

    /// Check if a geometry is cut by the brep
    int CutStatus(GeometryType::Pointer p_geom, const int& configuration) const;

    /// Check if a geometry is cut by the BRep
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    virtual int CutStatus(GeometryType& r_geom, const int& configuration) const;

    /// Check if a set of points is cut by the BRep
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    virtual int CutStatus(const std::vector<PointType>& r_points) const;

    /// Check if a set of points is cut by the BRep
    /// The geometry and the local information of the points are also provided.
    /// This function allows the use of BRep defined on grid (nodes)
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    virtual int CutStatus(const GeometryType& r_geom,
                          const std::vector<CoordinatesArrayType>& r_local_points,
                          const std::vector<PointType>& r_points) const;

    /// return the string of the cut status
    static std::string CutStatusStr(const int& stat);

    /// Check if an element is cut by the brep by sampling the elemental geometry
    int CutStatusBySampling(Element::Pointer p_elem, const std::size_t& nsampling, const int& configuration) const;

    /// Check if a geometry is cut by the brep by sampling the geometry
    int CutStatusBySampling(GeometryType::Pointer p_geom, const std::size_t& nsampling, const int& configuration) const;

    /// Check if a geometry is cut by the BRep by sampling the geometry
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    virtual int CutStatusBySampling(const GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const;

    /// Compute the intersection of the BRep with a line
    virtual int Bisect(PointType& P, const std::vector<PointType>& Points, const double& Tol) const
    {
        if (Points.size() == 2)
        {
            return Bisect(P, Points[0], Points[1], Tol);
        }
        else if (Points.size() == 3)
        {
            return Bisect(P, Points[0], Points[1], Points[2], Tol);
        }
        else
        {
            KRATOS_ERROR << __FUNCTION__ << ": Invalid number of points";
        }
        return -2;
    }

    /// Compute the intersection of the BRep with a line connected by 2 points
    virtual int Bisect(PointType& P, const PointType& P1, const PointType& P2, const double& Tol) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Compute the intersection of the BRep with a line connected by 3 points
    /// Remark: P3 is the middle point, the line is interpolated by using Lagrange shape function
    virtual int Bisect(PointType& P, const PointType& P1, const PointType& P2, const PointType& P3, const double& Tol) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Get/compute the normal vector at a point on the BRep
    virtual void GetNormal(const PointType& P, PointType& rNormal) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Get/compute the normal vector derivatives w.r.t the global point
    /// The derivatives are organized as:
    ///     [d N[0] / d P[0], d N[0] / d P[1], d N[0] / d P[2]]
    ///     [d N[1] / d P[0], d N[1] / d P[1], d N[1] / d P[2]]
    ///     [d N[2] / d P[0], d N[2] / d P[1], d N[2] / d P[2]]
    virtual void GetNormalDerivatives(const PointType& P, Matrix& Derivatives) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Get/compute the tangential vectors at a point on the BRep
    virtual void GetTangent(const PointType& P, std::vector<PointType>& rTangentialVectors) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Get/compute the tangential vectors derivatives w.r.t the global point
    /// The derivatives are organized as:
    /// [
    ///     [d T1[0] / d P[0], d T1[0] / d P[1], d T1[0] / d P[2]]
    ///     [d T1[1] / d P[0], d T1[1] / d P[1], d T1[1] / d P[2]]
    ///     [d T1[2] / d P[0], d T1[2] / d P[1], d T1[2] / d P[2]]
    /// ]
    /// [
    ///     [d T2[0] / d P[0], d T2[0] / d P[1], d T2[0] / d P[2]]
    ///     [d T2[1] / d P[0], d T2[1] / d P[1], d T2[1] / d P[2]]
    ///     [d T2[2] / d P[0], d T2[2] / d P[1], d T2[2] / d P[2]]
    /// ]
    virtual void GetTangentDerivatives(const PointType& P, std::vector<Matrix>& Derivatives) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// projects a point on the surface of brep
    virtual int ProjectOnSurface(const PointType& P, PointType& Proj) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// compute the derivatives of the projection point w.r.t to the original point.
    /// The derivatives are organized as;
    ///     [d Proj[0] / d P[0], d Proj[0] / d P[1], d Proj[0] / d P[2]]
    ///     [d Proj[1] / d P[0], d Proj[1] / d P[1], d Proj[1] / d P[2]]
    ///     [d Proj[2] / d P[0], d Proj[2] / d P[1], d Proj[2] / d P[2]]
    virtual void ProjectionDerivatives(const PointType& P, Matrix& Derivatives) const
    {
        KRATOS_ERROR << "Calling the base class";
    }

    /// Intersect the BRep with a geometry
    virtual Section::Pointer Intersect(GeometryType& r_geom) const;

    ///@}
    ///@name Access
    ///@{

    /// Set value of internal variable
    virtual void SetValue(const Variable<bool>& rVariable, const bool& rValue);

    /// Set value of internal variable
    virtual void SetValue(const Variable<int>& rVariable, const int& rValue);

    /// Set value of internal variable
    virtual void SetValue(const Variable<double>& rVariable, const double& rValue);

    /// Get value of internal variable
    virtual bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) const;

    /// Get value of internal variable
    virtual int& GetValue(const Variable<int>& rThisVariable, int& rValue) const;

    /// Get value of internal variable
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue) const;

    ///@}
    ///@name Inquiry
    ///@{

    template<class TPointsContainerType>
    int CutStatusOfPoints(const TPointsContainerType& r_points) const
    {
        std::vector<std::size_t> in_list, out_list;
        bool check;
        for (std::size_t v = 0; v < r_points.size(); ++v)
        {
            check = this->IsInside(r_points[v]);
            if (check)
            {
                in_list.push_back(v);
            }
            else
            {
                out_list.push_back(v);
            }
        }

        int stat;
        if (in_list.size() == 0 && out_list.size() == 0)
        {
            for (std::size_t v = 0; v < r_points.size(); ++v)
            {
                KRATOS_WATCH(r_points[v])
            }
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_ERROR << "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.";
        }
        else
        {
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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BRep";
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

    double mTOL;
    std::string mName;
    int mEchoLevel;

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
    BRep& operator=(BRep const& rOther);

    ///@}

}; // Class BRep

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, BRep& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const BRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BREP_H_INCLUDED  defined
