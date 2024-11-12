//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 25, 2014$
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   28/7/2015: create brep_application.cpp

// System includes

// External includes

// Project includes
#include "brep_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE( FunctionR3Rn::Pointer, LOAD_FUNCTION )
KRATOS_CREATE_VARIABLE( int, CUT_STATUS )
KRATOS_CREATE_VARIABLE( double, CURVE_SEARCH_TOLERANCE )
KRATOS_CREATE_VARIABLE( int, CURVE_MAX_ITERATIONS )
KRATOS_CREATE_VARIABLE( double, CURVE_LOWER_BOUND )
KRATOS_CREATE_VARIABLE( double, CURVE_UPPER_BOUND )
KRATOS_CREATE_VARIABLE( int, CURVE_NUMBER_OF_SAMPLING )
KRATOS_CREATE_VARIABLE( Vector, PLANE_EQUATION_PARAMETERS )
KRATOS_CREATE_VARIABLE( bool, IS_NODAL_LEVEL_SET )
KRATOS_CREATE_VARIABLE( double, LEVEL_SET_VALUE )

} // namespace Kratos
