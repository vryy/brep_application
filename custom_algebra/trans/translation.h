//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 21 Nov 2017 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_BREP_APPLICATION_TRANSLATION_H_INCLUDED )
#define  KRATOS_BREP_APPLICATION_TRANSLATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_algebra/trans/transformation.h"

namespace Kratos
{

/**
Represent a Translation in homogeneous coordinates
 */
template<typename TDataType>
class Translation : public Transformation<TDataType>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Translation);

    /// Type definitions
    typedef Transformation<TDataType> BaseType;
    typedef typename BaseType::MatrixType MatrixType;

    /// Default constructor
    Translation(const TDataType& tX, const TDataType& tY, const TDataType& tZ) : BaseType()
    {
        BaseType::mTransMat(0, 3) = tX;
        BaseType::mTransMat(1, 3) = tY;
        BaseType::mTransMat(2, 3) = tZ;
    }

    /// Destructor
    ~Translation() override {}

    /// Information
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "Homogeneous Translation";
    }

};

}// namespace Kratos.

#endif // KRATOS_BREP_APPLICATION_TRANSLATION_H_INCLUDED
