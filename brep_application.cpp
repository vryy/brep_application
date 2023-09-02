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
#include "brep_application.h"
#include "brep_application_variables.h"

namespace Kratos
{

KratosBRepApplication::KratosBRepApplication()
#ifdef SD_APP_FORWARD_COMPATIBILITY
    : KratosApplication("KratosBRepApplication")
#else
    : KratosApplication()
#endif
{}

void KratosBRepApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosBRepApplication... " << std::endl;

    // register variables to Kratos kernel
    //KRATOS_REGISTER_VARIABLE( LOAD_FUNCTION )
    KRATOS_REGISTER_VARIABLE( CUT_STATUS )
    KRATOS_REGISTER_VARIABLE( CURVE_SEARCH_TOLERANCE )
    KRATOS_REGISTER_VARIABLE( CURVE_MAX_ITERATIONS )
    KRATOS_REGISTER_VARIABLE( CURVE_LOWER_BOUND )
    KRATOS_REGISTER_VARIABLE( CURVE_UPPER_BOUND )
    KRATOS_REGISTER_VARIABLE( CURVE_NUMBER_OF_SAMPLING )
    KRATOS_REGISTER_VARIABLE( PLANE_EQUATION_PARAMETERS )
    KRATOS_REGISTER_VARIABLE( LEVEL_SET_VALUE )
}

} // namespace Kratos
