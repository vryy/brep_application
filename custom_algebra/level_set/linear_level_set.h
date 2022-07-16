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
//  Date:            14 Feb 2017
//


#if !defined(KRATOS_LINEAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_LINEAR_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/brep_mesh_utility.h"


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

///
/**
 * Level set representing signed-distance function to a line
 */
class LinearLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(LinearLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearLevelSet(const double& A, const double& B, const double& C)
    : BaseType()
    {
        double aux = sqrt(pow(A, 2) + pow(B, 2));

        if (aux < 1.0e-10)
            KRATOS_THROW_ERROR(std::logic_error, "The line is degenerated", "")

        mA = A / aux;
        mB = B / aux;
        mC = C / aux;
    }

    /// Copy constructor.
    LinearLevelSet(LinearLevelSet const& rOther)
    : BaseType(rOther), mA(rOther.mA), mB(rOther.mB), mC(rOther.mC)
    {}

    /// Destructor.
    virtual ~LinearLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new LinearLevelSet(*this));
    }


    std::size_t WorkingSpaceDimension() const final
    {
        return 2;
    }


    double GetValue(const PointType& P) const final
    {
        return mA*P(0) + mB*P(1) + mC;
    }


    Vector GetGradient(const PointType& P) const final
    {
        Vector grad(2);
        grad(0) = mA;
        grad(1) = mB;
        return grad;
    }


    /// inherit from BRep
    void GetTangent(const PointType& P, std::vector<PointType>& rTangentialVectors) const final
    {
        rTangentialVectors.resize(1);

        rTangentialVectors[0](0) = -mB;
        rTangentialVectors[0](1) = -mA;
        rTangentialVectors[0](2) = 0.0;
    }


    /// inherit from BRep
    void GetTangentDerivatives(const PointType& P, std::vector<Matrix>& Derivatives) const final
    {
        Derivatives.resize(1);
        Derivatives[0] = ZeroMatrix(3, 3);
    }


    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double t = -(mA*P[0] + mB*P[1] + mC) / (pow(mA, 2) + pow(mB, 2));
        Proj[0] = P[0] + mA*t;
        Proj[1] = P[1] + mB*t;
        Proj[2] = 0.0;
        return 0;
    }


    /***********EXCLUSIVE INTERFACE****************/


    /// Generate the sampling points on the level set surface
    void GeneratePoints(std::vector<PointType>& points,
        const PointType& StartPoint, const PointType& EndPoint,
        const std::size_t& nsampling) const
    {
        // check if the start point and end point are on the line
        double f;
        f = this->GetValue(StartPoint);
        if (f > this->Tolerance())
            KRATOS_THROW_ERROR(std::logic_error, "The start point does not lie on the line", "")
        f = this->GetValue(EndPoint);
        if (f > this->Tolerance())
            KRATOS_THROW_ERROR(std::logic_error, "The end point does not lie on the line", "")

        BRepMeshUtility::GenerateSamplingPoints(points, StartPoint, EndPoint, nsampling);
    }


    /// Create the conditions based on sampling points on the line
    std::pair<ModelPart::NodesContainerType, ModelPart::ConditionsContainerType> CreateLineConditions(ModelPart& r_model_part,
        const std::string& sample_condition_name,
        Properties::Pointer pProperties,
        const PointType& StartPoint,
        const PointType& EndPoint,
        const int& type,
        const std::size_t& nsampling) const
    {
        // firstly create the sampling points on surface
        std::vector<PointType> sampling_points;
        this->GeneratePoints(sampling_points, StartPoint, EndPoint, nsampling);
        // KRATOS_WATCH(sampling_points.size())
        const bool close = false;
        BRepMeshUtility::ConditionMeshInfoType Info = BRepMeshUtility::CreateLineConditions(r_model_part, sampling_points, sample_condition_name, type, close, pProperties);
        return std::make_pair(std::get<0>(Info), std::get<1>(Info));
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
        return "Linear Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "A: " << mA << ", B: " << mB << ", C: " << mC;
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


    double mA, mB, mC;


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
    LinearLevelSet& operator=(LinearLevelSet const& rOther);

    ///@}

}; // Class LinearLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                LinearLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const LinearLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LINEAR_LEVEL_SET_H_INCLUDED  defined
