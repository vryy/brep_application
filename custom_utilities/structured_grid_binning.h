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


#if !defined(KRATOS_STRUCTURED_GRID_BINNING_H_INCLUDED )
#define  KRATOS_STRUCTURED_GRID_BINNING_H_INCLUDED



// System includes
#include <set>


// External includes


// Project includes
#include "custom_utilities/grid_binning.h"


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
/** The binning for structured mesh.
 * The structured mesh can be skewed, i.e. having a parallelogram shape
 * The adaptive mesh generated by space-tree sub-division technique is also supported
 */
class StructuredGridBinning
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuredGridBinning
    KRATOS_CLASS_POINTER_DEFINITION(StructuredGridBinning);

    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef Element::GeometryType GeometryType;
    typedef GeometryType::PointType NodeType;
    typedef NodeType::PointType PointType;

    struct SpatialKey
    {
        SpatialKey(int ix, int iy, int iz) : x(ix), y(iy), z(iz) {}
        ~SpatialKey() {}
        bool operator<(const SpatialKey& rOther) const
        {
            if(x == rOther.x)
            {
                if(y == rOther.y)
                {
                    return z < rOther.z;
                }
                else
                    return y < rOther.y;
            }
            else
                return x < rOther.x;
        }
        int kx() const {return x;}
        int ky() const {return y;}
        int kz() const {return z;}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructuredGridBinning(const double& O_1, const double& O_2, const double& O_3, //origin point of the grid
        const double& V1_1, const double& V1_2, const double& V1_3, // alignment vector V1
        const double& V2_1, const double& V2_2, const double& V2_3, // alignment vector V2
        const double& V3_1, const double& V3_2, const double& V3_3) // alignment vector V3
    : mTolerance(1.0e-6), mpElements(NULL)
    {
        mOrigin[0] = O_1; mOrigin[1] = O_2; mOrigin[2] = O_3;
        mV1[0] = V1_1; mV1[1] = V1_2; mV1[2] = V1_3;
        mV2[0] = V2_1; mV2[1] = V2_2; mV2[2] = V2_3;
        mV3[0] = V3_1; mV3[1] = V3_2; mV3[2] = V3_3;
    }

    /// Destructor.
    virtual ~StructuredGridBinning()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Initialize the binning. It's important for all element edges to align with alignment vectors.
    void Initialize(const ElementsContainerType& rElements)
    {
        mpElements = &Elements;

        std::set<long long> V1_coords_set;
        std::set<long long> V2_coords_set;
        std::set<long long> V3_coords_set;

        // compute all the projection of element nodes onto alignment vectors
        double x, y, z, dx, dy, dz;
        for (ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            for (std::size_t i = 0; i < (*it)->GetGeometry().size(); ++i)
            {
                dx = ((*it)->GetGeometry()[i].X0() - mOrigin[0]);
                dy = ((*it)->GetGeometry()[i].Y0() - mOrigin[1]);
                dz = ((*it)->GetGeometry()[i].Z0() - mOrigin[2]);

                x = dx*mV1[0] + dy*mV1[1] + dz*mV1[2];
                y = dx*mV2[0] + dy*mV2[1] + dz*mV2[2];
                z = dx*mV3[0] + dy*mV3[1] + dz*mV3[2];

                V1_coords_set.insert(floor(x/mTolerance));
                V1_coords_set.insert(floor(y/mTolerance));
                V1_coords_set.insert(floor(z/mTolerance));
            }
        }

        mV1_coords.clear();
        mV2_coords.clear();
        mV3_coords.clear();

        std::copy(V1_coords_set.begin(), V1_coords_set.end(), mV1_coords.begin());
        std::copy(V2_coords_set.begin(), V2_coords_set.end(), mV2_coords.begin());
        std::copy(V3_coords_set.begin(), V3_coords_set.end(), mV3_coords.begin());

        // assign the spatial key for each element
        for (ElementsContainerType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            // TODO
        }
    }

    ///@}
    ///@name Access
    ///@{

    ElementsContainerType::iterator find(const PointType& P) override
    {
        // TODO
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
    }

    ElementsContainerType::iterator end() override
    {
        // TODO
        KRATOS_THROW_ERROR(std::logic_error, "Error calling base class function", __FUNCTION__)
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
        return "Grid Binning";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
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

    PointType mOrigin, mV1, mV2, mV3;
    const ElementsContainerType* mpElements;

    double mTolerance;
    std::set<long long> mV1_coords;
    std::set<long long> mV2_coords;
    std::set<long long> mV3_coords;
    // std::map<IndexType, SpatialKey> mElementKeyMap;

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
    StructuredGridBinning& operator=(StructuredGridBinning const& rOther);

    ///@}

}; // Class StructuredGridBinning

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, StructuredGridBinning& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const StructuredGridBinning& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_STRUCTURED_GRID_BINNING_H_INCLUDED  defined
