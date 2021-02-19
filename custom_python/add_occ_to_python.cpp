// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Aug 2019 $
//   Revision:            $Revision: 1.0 $
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "custom_python/add_occ_to_python.h"
#include "custom_utilities/occ_define.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void BRepApplication_AddOCCToPython()
{
    /**************************************************************/
    /************* EXPORT INTERFACE FOR OCC ***********************/
    /**************************************************************/
    /* SINCE WE DO NOT KNOW IF ANY EXTERNAL PACKAGES DO EXPORT ****/
    /* OCC STUFF, WE APPEND THE PREFIX KRATOS_ TO ALL OCC CLASS ***/
    /* TO MAKE SURE THAT THERE ARE NO NAMING CONFLICTS ************/
    /**************************************************************/

    class_<TopoDS_Shape, OCC::shared_ptr<TopoDS_Shape>, boost::noncopyable>
    ( "Kratos_TopoDS_Shape", init<>() )
    .def("__str__", &OCC::PrintObject<TopoDS_Shape>)
    ;

}
}  // namespace Python.
}  // namespace Kratos.

