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

    typedef BaseType::NodeType NodeType;

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
    ~OCCBRep() override
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
    BRep::Pointer CloneOCCBRep() const
    {
        return BRep::Pointer(new OCCBRep());
    }

    /// Get working space dimension
    std::size_t WorkingSpaceDimension() const final
    {
        // TODO
        KRATOS_ERROR << "Calling the base " << __FUNCTION__;
    }

    /// Get local space dimension
    std::size_t LocalSpaceDimension() const final
    {
        // TODO
        KRATOS_ERROR << "Calling the base " << __FUNCTION__;
    }

    /// Check if a point is inside/outside of the OCCBRep
    // REF: https://www.opencascade.com/content/point-inside-solid-or-not
    bool IsInside(const PointType& P) const final
    {
        BRepClass3d_SolidClassifier solidClassifier(*mpShape);
        gp_Pnt occPoint(P[0], P[1], P[2]);
        solidClassifier.Perform(occPoint, this->GetTolerance());
        return (solidClassifier.State() == TopAbs_State::TopAbs_IN);
    }

    /// Check if a point is on the boundary within a tolerance
    bool IsOnBoundary(const PointType& P) const final
    {
        BRepClass3d_SolidClassifier solidClassifier(*mpShape);
        gp_Pnt occPoint(P[0], P[1], P[2]);
        solidClassifier.Perform(occPoint, this->Tolerance());
        return (solidClassifier.State() == TopAbs_State::TopAbs_ON);
    }

    /// Compute the intersection of the OCCBRep with a line connect by 2 points.
    int Bisect(PointType& P, const PointType& P1, const PointType& P2, const double& tol) const final
    {
        KRATOS_ERROR << "Calling the base " << __FUNCTION__;
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
        return "OpenCasCade BRep";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        if (mpShape != NULL)
        {
            OCC::DumpObject(rOStream, *mpShape);
        }
        else
        {
            rOStream << "OCCBRep contains no object";
        }
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
