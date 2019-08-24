//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         OCCBRep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            23 Aug 2019
//


#if !defined(KRATOS_OPENCASCADE_BREP_H_INCLUDED )
#define  KRATOS_OPENCASCADE_BREP_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <BRepClass3d_SolidClassifier.hxx>


// Project includes
#include "includes/define.h"
#include "custom_algebra/brep.h"
#include "custom_utilities/occ_define.h"


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
/**
 * BRep represented by OpenCasCade BRep object
 */
class OCCBRep : public BRep
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OCCBRep
    KRATOS_CLASS_POINTER_DEFINITION(OCCBRep);

    typedef BRep BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::PointType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OCCBRep() : BaseType()
    {}

    /// Copy constructor.
    OCCBRep(OCCBRep const& rOther) : BaseType(rOther), mpShape(rOther.mpShape)
    {}

    /// Destructor.
    virtual ~OCCBRep()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Set the underlying OCC geometry
    void SetShape(OCC::shared_ptr<TopoDS_Shape> pShape)
    {
        mpShape = pShape;
    }

    /// Clone this OCCBRep
    virtual BRep::Pointer CloneOCCBRep() const
    {
        return BRep::Pointer(new OCCBRep());
    }

    /// Get working space dimension
    virtual std::size_t WorkingSpaceDimension() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// Get local space dimension
    virtual std::size_t LocalSpaceDimension() const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
    }

    /// Check if a point is inside/outside of the OCCBRep
    // REF: https://www.opencascade.com/content/point-inside-solid-or-not
    virtual bool IsInside(const PointType& P) const
    {
        BRepClass3d_SolidClassifier solidClassifier(*mpShape);
        gp_Pnt occPoint(P[0], P[1], P[2]);
        solidClassifier.Perform(occPoint, this->GetTolerance());
        return (solidClassifier.State() == TopAbs_State::TopAbs_IN);
    }

    /// Check if a point is on the boundary within a tolerance
    virtual bool IsOnBoundary(const PointType& P, const double& tol) const
    {
        BRepClass3d_SolidClassifier solidClassifier(*mpShape);
        gp_Pnt occPoint(P[0], P[1], P[2]);
        solidClassifier.Perform(occPoint, tol);
        return (solidClassifier.State() == TopAbs_State::TopAbs_ON);
    }

    /// Check if a geometry is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(GeometryType& r_geom, const int& configuration) const
    {
        if (configuration == 0)
        {
            std::vector<PointType> points(r_geom.size());
            for (std::size_t i = 0; i < r_geom.size(); ++i)
                noalias(points[i]) = r_geom[i].GetInitialPosition();
            return CutStatusOfPoints(points);
        }
        else if (configuration == 1)
        {
            return CutStatusOfPoints(r_geom);
            // REMARK: this will use the current position of node, e.g. in dynamics
        }
    }

    /// Check if a set of points is cut by the level set
    /// 0: the cell is completely inside the domain bounded by level set
    /// 1: completely outside
    /// -1: the cell is cut by level set
    virtual int CutStatus(const std::vector<PointType>& r_points) const
    {
        return CutStatusOfPoints(r_points);
    }

    /// Compute the intersection of the OCCBRep with a line connect by 2 points.
    virtual PointType Bisect(const PointType& P1, const PointType& P2, const double& tol) const
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
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
        return "OpenCasCade BRep";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        if (mpShape != NULL)
            OCC::DumpObject(rOStream, *mpShape);
        else
            rOStream << "OCCBRep contains no object";
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

    OCC::shared_ptr<TopoDS_Shape> mpShape;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template<class TPointsContainerType>
    int CutStatusOfPoints(const TPointsContainerType& r_points) const
    {
        std::vector<std::size_t> in_list, out_list;
        bool check;
        for(std::size_t v = 0; v < r_points.size(); ++v)
        {
            check = this->IsInside(r_points[v]);
            if(check)
                in_list.push_back(v);
            else
                out_list.push_back(v);
        }

        int stat;
        if(in_list.size() == 0 && out_list.size() == 0)
        {
            for(std::size_t v = 0; v < r_points.size(); ++v)
                KRATOS_WATCH(r_points[v])
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_THROW_ERROR(std::logic_error, "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.", "")
        }
        else
        {
            if(in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if(out_list.size() == 0)
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
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    OCCBRep& operator=(OCCBRep const& rOther);

    ///@}

}; // Class OCCBRep

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, OCCBRep& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const OCCBRep& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_OPENCASCADE_BREP_H_INCLUDED  defined
