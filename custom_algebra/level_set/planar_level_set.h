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


#if !defined(KRATOS_PLANAR_LEVEL_SET_H_INCLUDED )
#define  KRATOS_PLANAR_LEVEL_SET_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
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
/** Detail class definition.
*/
class PlanarLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PlanarLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(PlanarLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::PointType PointType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PlanarLevelSet(const double& A, const double& B, const double& C, const double& D)
    : BaseType(), mA(A), mB(B), mC(C), mD(D)
    {}

    /// Contructor with origin and normal vector
    PlanarLevelSet(const PointType& Origin, const PointType& Normal)
    : BaseType()
    {
        double norm = norm_2(Normal);
        mA = Normal[0] / norm;
        mB = Normal[1] / norm;
        mC = Normal[2] / norm;
        mD = -(mA*Origin[0] + mB*Origin[1] + mC*Origin[2]);
    }

    /// Copy constructor.
    PlanarLevelSet(PlanarLevelSet const& rOther)
    : BaseType(rOther), mA(rOther.mA), mB(rOther.mB), mC(rOther.mC), mD(rOther.mD)
    {}

    /// Destructor.
    virtual ~PlanarLevelSet() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new PlanarLevelSet(*this));
    }


    std::size_t WorkingSpaceDimension() const final
    {
        return 3;
    }


    void NormalVector(array_1d<double, 3>& V) const
    {
        V[0] = mA;
        V[1] = mB;
        V[2] = mC;
    }


    double GetValue(const PointType& P) const final
    {
        return mA*P(0) + mB*P(1) + mC*P(2) + mD;
    }


    Vector GetGradient(const PointType& P) const final
    {
        Vector grad(this->WorkingSpaceDimension());
        grad(0) = mA;
        grad(1) = mB;
        grad(2) = mC;
        return grad;
    }


    Matrix GetGradientDerivatives(const PointType& P) const final
    {
        Matrix Jac(this->WorkingSpaceDimension(), 3);
        noalias(Jac) = ZeroMatrix(this->WorkingSpaceDimension(), 3);
        return Jac;
    }


    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        double t = -(mA*P[0] + mB*P[1] + mC*P[2] + mD) / (pow(mA, 2) + pow(mB, 2) + pow(mC, 2));
        Proj[0] = P[0] + mA*t;
        Proj[1] = P[1] + mB*t;
        Proj[2] = P[2] + mC*t;
        return 0;
    }


    /// compute the derivatives of the projection point w.r.t to the original point.
    /// The derivatives are organized as;
    ///     [d Proj[0] / d P[0], d Proj[0] / d P[1], d Proj[0] / d P[2]]
    ///     [d Proj[1] / d P[0], d Proj[1] / d P[1], d Proj[1] / d P[2]]
    ///     [d Proj[2] / d P[0], d Proj[2] / d P[1], d Proj[2] / d P[2]]
    void ProjectionDerivatives(const PointType& P, Matrix& Derivatives) const final
    {
        if (Derivatives.size1() != 3 || Derivatives.size2() != 3)
            Derivatives.resize(3, 3, false);

        double aux = (pow(mA, 2) + pow(mB, 2) + pow(mC, 2));

        Derivatives(0, 0) = 1.0 - pow(mA, 2) / aux;
        Derivatives(0, 1) = -mA*mB / aux;
        Derivatives(0, 2) = -mA*mC / aux;

        Derivatives(1, 0) = -mB*mA / aux;
        Derivatives(1, 1) = 1.0 - pow(mB, 2) / aux;
        Derivatives(1, 2) = -mB*mC / aux;

        Derivatives(2, 0) = -mC*mA / aux;
        Derivatives(2, 1) = -mC*mB / aux;
        Derivatives(2, 2) = 1.0 - pow(mC, 2) / aux;
    }


    /// inherit from BRep
    /// Compute the intersection of the level set with a line connect by 2 points.
    /// Note that, the checking of the intersection of the level set with the line is not performed. Hence one should ensure that before calling this function.
    int Bisect(PointType& P, const PointType& P1, const PointType& P2, const double& Tol) const final
    {
        double f1 = this->GetValue(P1);
        double f2 = this->GetValue(P2);
        if(f1*f2 > 0.0)
        {
            // KRATOS_THROW_ERROR(std::logic_error, "Bisect does not work with two ends at the same side", "")
            return -1;
        }

        noalias(P) = (f2*P1-f1*P2)/(f2-f1);
        return 0;
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
        return "Planar Level Set";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << mA << "*x + " << mB << "*y + " << mC << "*z + " << mD;
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


    double mA, mB, mC, mD;


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
    PlanarLevelSet& operator=(PlanarLevelSet const& rOther);


    ///@}

}; // Class PlanarLevelSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                PlanarLevelSet& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PlanarLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PLANAR_LEVEL_SET_H_INCLUDED  defined
