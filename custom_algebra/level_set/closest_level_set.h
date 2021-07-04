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
//  Date:            12 May 2021
//


#if !defined(KRATOS_CLOSEST_LEVEL_SET_H_INCLUDED )
#define  KRATOS_CLOSEST_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "custom_algebra/level_set/level_set.h"


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
/** Select the closest level set between two
*/
class ClosestLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ClosestLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(ClosestLevelSet);

    typedef LevelSet BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ClosestLevelSet(const BaseType::Pointer p_level_set_1, const BaseType::Pointer p_level_set_2)
    : BaseType(), mp_level_set_1(p_level_set_1), mp_level_set_2(p_level_set_2)
    {}

    /// Copy constructor.
    ClosestLevelSet(ClosestLevelSet const& rOther)
    : BaseType(rOther)
    , mp_level_set_1(rOther.mp_level_set_1->CloneLevelSet())
    , mp_level_set_2(rOther.mp_level_set_2->CloneLevelSet())
    {}

    /// Destructor.
    virtual ~ClosestLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new ClosestLevelSet(*this));
    }


    std::size_t WorkingSpaceDimension() const final
    {
        return mp_level_set_1->WorkingSpaceDimension();
    }


    double GetValue(const PointType& P) const final
    {
        if (std::abs(mp_level_set_1->GetValue(P)) < std::abs(mp_level_set_2->GetValue(P)))
            return mp_level_set_1->GetValue(P);
        else
            return mp_level_set_2->GetValue(P);
    }


    Vector GetGradient(const PointType& P) const final
    {
        if (std::abs(mp_level_set_1->GetValue(P)) < std::abs(mp_level_set_2->GetValue(P)))
            return mp_level_set_1->GetGradient(P);
        else
            return mp_level_set_2->GetGradient(P);
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
        ss << "(Closest Level Set of " << mp_level_set_1->Info() << " and " << mp_level_set_2->Info() << ")";
        return ss.str();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        // rOStream << "(" << Info() << " of ";
        // mp_level_set_1->PrintData(rOStream);
        // rOStream << " and ";
        // mp_level_set_1->PrintData(rOStream);
        // rOStream << ")";
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


    const BaseType::Pointer mp_level_set_1;
    const BaseType::Pointer mp_level_set_2;


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
    ClosestLevelSet& operator=(ClosestLevelSet const& rOther);

    ///@}

}; // Class ClosestLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ClosestLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ClosestLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CLOSEST_LEVEL_SET_H_INCLUDED  defined
