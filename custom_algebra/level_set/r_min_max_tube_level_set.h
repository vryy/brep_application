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
//  Date:            8 Nov 2024
//

#if !defined(KRATOS_TUBE_LEVEL_SET_RMIN_RMAX_H_INCLUDED )
#define  KRATOS_TUBE_LEVEL_SET_RMIN_RMAX_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>

// External includes

// Project includes
#include "custom_algebra/level_set/level_set.h"
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

/// A specialized Level set defined by an alignment curve and a tube LevelSet
/// If the distance to the alignment curve is smaller than Rmin -> Vmin
/// If the distance to the alignment curve is larger than Rmax -> Vmax
/// The remaining case (Rmin < R < Rmax) is governed by the tube LevelSet
class RminmaxTubeLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RminmaxTubeLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(RminmaxTubeLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RminmaxTubeLevelSet(const LevelSet::Pointer pLevelSet, const Curve::Pointer pAlignCurve,
        const double Rmin, const double Rmax, const double Vmin, const double Vmax)
        : mpLevelSet(pLevelSet), mpCurve(pAlignCurve)
        , mRmin(Rmin), mRmax(Rmax), mVmin(Vmin), mVmax(Vmax)
        , BaseType()
    {}

    /// Copy constructor.
    RminmaxTubeLevelSet(RminmaxTubeLevelSet const& rOther)
        : BaseType(rOther)
        , mpLevelSet(rOther.mpLevelSet->CloneLevelSet())
        , mpCurve(rOther.mpCurve->Clone())
        , mRmin(rOther.mRmin)
        , mRmax(rOther.mRmax)
        , mVmin(rOther.mVmin)
        , mVmax(rOther.mVmax)
    {}

    /// Destructor.
    ~RminmaxTubeLevelSet() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new RminmaxTubeLevelSet(*this));
    }

    BRep::ConstPointer pBRep() const final
    {
        return mpLevelSet;
    }

    std::size_t WorkingSpaceDimension() const final
    {
        return mpLevelSet->WorkingSpaceDimension();
    }

    std::size_t LocalSpaceDimension() const final
    {
        return mpLevelSet->LocalSpaceDimension();
    }

    double GetValue(const PointType& P) const final
    {
        const double R = mpCurve->ComputeDistance(P);

        if (R < mRmin) return mVmin;
        if (R > mRmax) return mVmax;

        return mpLevelSet->GetValue(P);
    }

    /// Check if a point is inside/outside of the LevelSet
    bool IsInside(const PointType& P) const final
    {
        const double R = mpCurve->ComputeDistance(P);

        if (R < mRmin) return true;
        if (R > mRmax) return false;

        return mpLevelSet->IsInside(P);
    }

    ///@}
    ///@name Access
    ///@{

    void SetValue(const Variable<bool>& rVariable, const bool& rValue) override
    {
        mpLevelSet->SetValue(rVariable, rValue);
    }

    void SetValue(const Variable<int>& rVariable, const int& rValue) override
    {
        mpLevelSet->SetValue(rVariable, rValue);
    }

    void SetValue(const Variable<double>& rVariable, const double& rValue) override
    {
        mpLevelSet->SetValue(rVariable, rValue);
    }

    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) const override
    {
        return mpLevelSet->GetValue(rThisVariable, rValue);
    }

    int& GetValue(const Variable<int>& rThisVariable, int& rValue) const override
    {
        return mpLevelSet->GetValue(rThisVariable, rValue);
    }

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) const override
    {
        return mpLevelSet->GetValue(rThisVariable, rValue);
    }

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
        ss << "RminmaxTubeLevelSet of (" << mpLevelSet->Info() << " and " << mpCurve->Info() << ")";
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
        mpLevelSet->PrintData(rOStream);
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

    LevelSet::Pointer mpLevelSet;
    Curve::Pointer mpCurve;
    double mRmin, mRmax, mVmin, mVmax;

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
    RminmaxTubeLevelSet& operator=(RminmaxTubeLevelSet const& rOther);

    ///@}

}; // Class RminmaxTubeLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TUBE_LEVEL_SET_RMIN_RMAX_H_INCLUDED  defined
