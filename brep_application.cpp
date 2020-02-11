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
// #include "geometries/triangle_2d_3.h"


namespace Kratos
{

    KRATOS_CREATE_VARIABLE( boost::python::object, LOAD_FUNCTION )
    KRATOS_CREATE_VARIABLE( int, CUT_STATUS )
    KRATOS_CREATE_VARIABLE( double, CURVE_SEARCH_TOLERANCE )
    KRATOS_CREATE_VARIABLE( int, CURVE_MAX_ITERATIONS )
    KRATOS_CREATE_VARIABLE( double, CURVE_LOWER_BOUND )
    KRATOS_CREATE_VARIABLE( double, CURVE_UPPER_BOUND )
    KRATOS_CREATE_VARIABLE( int, CURVE_NUMBER_OF_SAMPLING )

    KratosBRepApplication::KratosBRepApplication()
    : KratosApplication()
    {}

    void KratosBRepApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosBRepApplication... " << std::endl;

        // register variables to Kratos kernel
        KRATOS_REGISTER_VARIABLE( LOAD_FUNCTION )
        KRATOS_REGISTER_VARIABLE( CUT_STATUS )
        KRATOS_REGISTER_VARIABLE( CURVE_SEARCH_TOLERANCE )
        KRATOS_REGISTER_VARIABLE( CURVE_MAX_ITERATIONS )
        KRATOS_REGISTER_VARIABLE( CURVE_LOWER_BOUND )
        KRATOS_REGISTER_VARIABLE( CURVE_UPPER_BOUND )
        KRATOS_REGISTER_VARIABLE( CURVE_NUMBER_OF_SAMPLING )
    }

} // namespace Kratos

