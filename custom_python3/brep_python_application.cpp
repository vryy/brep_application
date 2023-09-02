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

// Project includes
#include "includes/define_python.h"
#include "brep_application.h"
#include "brep_application_variables.h"
#include "custom_python3/add_custom_algebra_to_python.h"
#include "custom_python3/add_brep_and_level_set_to_python.h"
#include "custom_python3/add_transformation_to_python.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_python3/add_occ_to_python.h"
#endif
#include "custom_python3/add_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosBRepApplication, m)
{

    py::class_<KratosBRepApplication, KratosBRepApplication::Pointer, KratosApplication>
    (m, "KratosBRepApplication")
    .def(py::init<>())
    ;

    BRepApplication_AddFunctionsToPython(m);
    BRepApplication_AddBRepAndLevelSetToPython(m);
    BRepApplication_AddTransformationToPython(m);
#ifdef BREP_APPLICATION_USE_OPENCASCADE
    BRepApplication_AddOCCToPython(m);
#endif
    BRepApplication_AddUtilitiesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LOAD_FUNCTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CUT_STATUS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CURVE_SEARCH_TOLERANCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CURVE_MAX_ITERATIONS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CURVE_LOWER_BOUND )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CURVE_UPPER_BOUND )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, CURVE_NUMBER_OF_SAMPLING )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, PLANE_EQUATION_PARAMETERS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( m, LEVEL_SET_VALUE )

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
