// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Aug 2019 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes


// Project includes
#include "custom_python/add_utilities_to_python.h"
#ifdef BREP_APPLICATION_USE_OPENCASCADE
#include "custom_utilities/occ_utility.h"
#endif

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void BRepApplication_AddUtilitiesToPython()
{

    #ifdef BREP_APPLICATION_USE_OPENCASCADE
    class_<OCCUtility, OCCUtility::Pointer, boost::noncopyable>
    ( "OCCUtility", init<>() )
    .def("MakeBottle", &OCCUtility::MakeBottle)
    .def("MakeSphere", &OCCUtility::MakeSphere)
    .def("ReadSTEP", &OCCUtility::ReadSTEP)
    .def("WriteSTEP", &OCCUtility::WriteSTEP)
    .def(self_ns::str(self))
    ;
    #endif

}
}  // namespace Python.
}  // namespace Kratos.

