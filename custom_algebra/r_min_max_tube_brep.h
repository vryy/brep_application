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
//  Date:            24 Oct 2024
//

#if !defined(KRATOS_TUBE_BREP_RMIN_RMAX_H_INCLUDED )
#define  KRATOS_TUBE_BREP_RMIN_RMAX_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>

// External includes

// Project includes
#include "custom_algebra/brep.h"
#include "custom_algebra/curve/curve.h"

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

/// A specialized BRep defined by an alignment curve and a tube BRep
/// If the distance to the alignment curve is smaller than Rmin -> In
/// If the distance to the alignment curve is larger than Rmax -> Out
/// The remaining case (Rmin < R < Rmax) is governed by the tube Brep
class RminmaxTubeBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RminmaxTubeBRep
    KRATOS_CLASS_POINTER_DEFINITION(RminmaxTubeBRep);

    typedef BRep BaseType;

    typedef BaseType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RminmaxTubeBRep(const BRep::Pointer pBRep, const Curve::Pointer pAlignCurve, const double Rmin, const double Rmax)
        : mpBRep(pBRep), mpCurve(pAlignCurve), mRmin(Rmin), mRmax(Rmax), BaseType()
    {}

    /// Copy constructor.
    RminmaxTubeBRep(RminmaxTubeBRep const& rOther)
        : BaseType(rOther)
        , mpBRep(rOther.mpBRep->CloneBRep())
        , mpCurve(rOther.mpCurve->Clone())
        , mRmin(rOther.mRmin)
        , mRmax(rOther.mRmax)
    {}

    /// Destructor.
    ~RminmaxTubeBRep() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    BRep::Pointer CloneBRep() const final
    {
        return BRep::Pointer(new RminmaxTubeBRep(*this));
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return mpBRep->WorkingSpaceDimension();
    }

    std::size_t LocalSpaceDimension() const final
    {
        return mpBRep->LocalSpaceDimension();
    }

    /// Check if a point is inside/outside of the BRep
    bool IsInside(const PointType& P) const final
    {
        const double R = mpCurve->ComputeDistance(P);
        // KRATOS_WATCH(R)
        // KRATOS_WATCH(mRmin)
        // KRATOS_WATCH(mRmax)

        // PointType proj;
        // double t;
        // mpCurve->ProjectOnCurve(P, proj, t);
        // KRATOS_WATCH(P)
        // KRATOS_WATCH(proj)
        // KRATOS_WATCH(t)
        // KRATOS_WATCH(norm_2(P - proj))

        // const auto v = mpCurve->GetDerivative(0, t);
        // KRATOS_WATCH(v)
        // KRATOS_WATCH(inner_prod(P - proj, v))

        // std::cout << "curve sampling:" << std::endl;
        // const unsigned int nsampling = 300;
        // for (unsigned int i = 0; i < nsampling; ++i)
        // {
        //     const double t = ((double)i / nsampling);
        //     const auto p = mpCurve->GetValue(t);
        //     std::cout << " " << p[0] << "," << p[1] << "," << p[2] << std::endl;
        // }

        if (R < mRmin) return true;
        if (R > mRmax) return false;

        return mpBRep->IsInside(P);
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
        std::stringstream ss;
        ss << "RminmaxTubeBRep of (" << mpBRep->Info() << " and " << mpCurve->Info() << ")";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        mpBRep->PrintData(rOStream);
        rOStream << std::endl;
        mpCurve->PrintData(rOStream);
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

    BRep::Pointer mpBRep;
    Curve::Pointer mpCurve;
    double mRmin, mRmax;

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
    RminmaxTubeBRep& operator=(RminmaxTubeBRep const& rOther);

    ///@}

}; // Class RminmaxTubeBRep

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TUBE_BREP_RMIN_RMAX_H_INCLUDED  defined
