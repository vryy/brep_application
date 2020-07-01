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
//  Date:            9 Sep 2017
//


#if !defined(KRATOS_NOT_BREP_H_INCLUDED )
#define  KRATOS_NOT_BREP_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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

/// BRep representing NOT operation
class NotBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NotBRep
    KRATOS_CLASS_POINTER_DEFINITION(NotBRep);

    typedef BRep BaseType;

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;

    typedef typename NodeType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NotBRep(BRep::Pointer pBRep)
    : mpBRep(pBRep), BaseType()
    {}

    /// Copy constructor.
    NotBRep(NotBRep const& rOther)
    : BaseType(rOther)
    , mpBRep(rOther.mpBRep->CloneBRep())
    {}

    /// Destructor.
    virtual ~NotBRep() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual BRep::Pointer CloneBRep() const
    {
        return BRep::Pointer(new NotBRep(*this));
    }


    virtual std::size_t WorkingSpaceDimension() const
    {
        return mpBRep->WorkingSpaceDimension();
    }


    virtual std::size_t LocalSpaceDimension() const
    {
        return mpBRep->LocalSpaceDimension();
    }

    virtual bool IsInside(const PointType& P) const
    {
        return !mpBRep->IsInside(P);
    }

    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(GeometryType& r_geom, const int& configuration) const
    {
        if(mpBRep->CutStatus(r_geom, configuration) == _OUT)
        {
            return _IN;
        }
        else if(mpBRep->CutStatus(r_geom, configuration) == _IN)
        {
            return _OUT;
        }
        else
        {
            return _CUT;
        }
    }

    /// Check if a set of points is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(const std::vector<PointType>& r_points) const
    {
        if(mpBRep->CutStatus(r_points) == _OUT)
        {
            return _IN;
        }
        else if(mpBRep->CutStatus(r_points) == _IN)
        {
            return _OUT;
        }
        else
        {
            return _CUT;
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
        return "NOT operation of a BRep";
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

    BRep::Pointer mpBRep;

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
    NotBRep& operator=(NotBRep const& rOther);

    ///@}

}; // Class NotBRep

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NotBRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NotBRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NOT_BREP_H_INCLUDED  defined
