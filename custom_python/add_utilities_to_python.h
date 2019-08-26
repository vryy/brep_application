// see brep_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23 Aug 2019 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_BREP_APPLICATION_ADD_UTILITIES_TO_PYTHON_H_INCLUDED )
#define  KRATOS_BREP_APPLICATION_ADD_UTILITIES_TO_PYTHON_H_INCLUDED


// System includes
#include <boost/python.hpp>


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  BRepApplication_AddUtilitiesToPython();

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_BREP_APPLICATION_ADD_UTILITIES_TO_PYTHON_H_INCLUDED  defined
