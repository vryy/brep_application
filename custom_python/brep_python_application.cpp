//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Oct 15, 2014 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "brep_application.h"
#include "brep_application_variables.h"
#include "custom_python/add_custom_algebra_to_python.h"
#include "custom_python/add_brep_and_level_set_to_python.h"
#include "custom_python/add_transformation_to_python.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_python/add_occ_to_python.h"
#endif
#include "custom_python/add_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosBRepApplication)
    {

        class_<KratosBRepApplication, KratosBRepApplication::Pointer, bases<KratosApplication>, boost::noncopyable>
        ("KratosBRepApplication");

        BRepApplication_AddFunctionsToPython();
        BRepApplication_AddBRepAndLevelSetToPython();
        BRepApplication_AddTransformationToPython();
        #ifdef BREP_APPLICATION_USE_OPENCASCADE
        BRepApplication_AddOCCToPython();
        #endif
        BRepApplication_AddUtilitiesToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOAD_FUNCTION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CUT_STATUS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CURVE_SEARCH_TOLERANCE )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CURVE_MAX_ITERATIONS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CURVE_LOWER_BOUND )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CURVE_UPPER_BOUND )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CURVE_NUMBER_OF_SAMPLING )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PLANE_EQUATION_PARAMETERS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LEVEL_SET_VALUE )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
