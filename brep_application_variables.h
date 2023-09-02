//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Oct 25, 2014 $
//   Revision:            $Revision: 1.0 $
//
//
//Change log:
//  +   28/7/2015: create brep_application.h

#if !defined(KRATOS_BREP_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_BREP_APPLICATION_VARIABLES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"
#include "custom_algebra/function/function.h"

namespace Kratos
{

///@name Kratos Globals
///@{

#ifdef SD_APP_FORWARD_COMPATIBILITY

KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, FunctionR3Rn::Pointer, LOAD_FUNCTION )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, int, CUT_STATUS )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, double, CURVE_SEARCH_TOLERANCE )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, int, CURVE_MAX_ITERATIONS )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, double, CURVE_LOWER_BOUND )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, double, CURVE_UPPER_BOUND )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, int, CURVE_NUMBER_OF_SAMPLING )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, Vector, PLANE_EQUATION_PARAMETERS )
KRATOS_DEFINE_APPLICATION_VARIABLE( BREP_APPLICATION, double, LEVEL_SET_VALUE )

#else

KRATOS_DEFINE_VARIABLE( FunctionR3Rn::Pointer, LOAD_FUNCTION )
KRATOS_DEFINE_VARIABLE( int, CUT_STATUS )
KRATOS_DEFINE_VARIABLE( double, CURVE_SEARCH_TOLERANCE )
KRATOS_DEFINE_VARIABLE( int, CURVE_MAX_ITERATIONS )
KRATOS_DEFINE_VARIABLE( double, CURVE_LOWER_BOUND )
KRATOS_DEFINE_VARIABLE( double, CURVE_UPPER_BOUND )
KRATOS_DEFINE_VARIABLE( int, CURVE_NUMBER_OF_SAMPLING )
KRATOS_DEFINE_VARIABLE( Vector, PLANE_EQUATION_PARAMETERS )
KRATOS_DEFINE_VARIABLE( double, LEVEL_SET_VALUE )

#endif

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

} // namespace Kratos

#endif // KRATOS_BREP_APPLICATION_VARIABLES_H_INCLUDED defined
