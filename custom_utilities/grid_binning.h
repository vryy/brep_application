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
//  Date:            25 Feb 2021
//

#if !defined(KRATOS_GRID_BINNING_H_INCLUDED )
#define  KRATOS_GRID_BINNING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

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
/** Abstract class for all type of grid binning algorithm.
 * The purpose of binning is to search for which element containing a point quickly
 */
class GridBinning
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GridBinning
    KRATOS_CLASS_POINTER_DEFINITION(GridBinning);

    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef NodeType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GridBinning() {}

    /// Destructor.
    virtual ~GridBinning() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    virtual ElementsContainerType::iterator find(const PointType& P)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    virtual ElementsContainerType::iterator end()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Grid Binning";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
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
    GridBinning& operator=(GridBinning const& rOther);

    ///@}

}; // Class GridBinning

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, GridBinning& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const GridBinning& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GRID_BINNING_H_INCLUDED  defined
