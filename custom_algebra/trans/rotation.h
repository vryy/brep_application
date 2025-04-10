//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_BREP_APPLICATION_ROTATION_H_INCLUDED )
#define  KRATOS_BREP_APPLICATION_ROTATION_H_INCLUDED

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

// System includes
#include <cmath>

// External includes

// Project includes
#include "custom_algebra/trans/transformation.h"

// On Visual Studio M_PI is available in cmath when _USE_MATH_DEFINES is defined
// but for any reason, it's not, then we define it here
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Kratos
{

template<int TAxis, typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        KRATOS_ERROR << "Error calling unimplemented function";
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<0, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle / 180.0 * M_PI);
        TDataType s = std::sin(angle / 180.0 * M_PI);
        trans_mat(1, 1) = c;
        trans_mat(1, 2) = -s;
        trans_mat(2, 1) = s;
        trans_mat(2, 2) = c;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<1, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle / 180.0 * M_PI);
        TDataType s = std::sin(angle / 180.0 * M_PI);
        trans_mat(0, 0) = c;
        trans_mat(0, 2) = s;
        trans_mat(2, 0) = -s;
        trans_mat(2, 2) = c;
    }
};

template<typename TMatrixType, typename TDataType>
struct ComputeRotationalTransformationMatrix_Helper<2, TMatrixType, TDataType>
{
    static void Execute(TMatrixType& trans_mat, const TDataType& angle)
    {
        TDataType c = std::cos(angle / 180.0 * M_PI);
        TDataType s = std::sin(angle / 180.0 * M_PI);
        trans_mat(0, 0) = c;
        trans_mat(0, 1) = -s;
        trans_mat(1, 0) = s;
        trans_mat(1, 1) = c;
    }
};

/**
 * Represent a Rotation in homogeneous coordinates
 * TAxis represent the axis of rotation
 */
template<int TAxis, typename TDataType>
class Rotation : public Transformation<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Rotation);

    /// Type definitions
    typedef Transformation<TDataType> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

    /// Default constructor
    Rotation(const TDataType& angle) : BaseType()
    {
        ComputeRotationalTransformationMatrix_Helper<TAxis, MatrixType, TDataType>::Execute(BaseType::mTransMat, angle);
    }

    /// Destructor
    ~Rotation() override {}

    /// Information
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "Homogeneous Rotation";
        if (TAxis == 0)
        {
            rOStream << "_X";
        }
        else if (TAxis == 1)
        {
            rOStream << "_Y";
        }
        else if (TAxis == 2)
        {
            rOStream << "_Z";
        }
    }

};

}// namespace Kratos.

#ifdef _MSC_VER
#undef _USE_MATH_DEFINES
#endif

#endif // KRATOS_BREP_APPLICATION_ROTATION_H_INCLUDED
