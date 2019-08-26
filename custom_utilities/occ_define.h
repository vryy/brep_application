#if !defined(KRATOS_BREP_APPLICATION_OPENCASCADE_DEFINE_H_INCLUDED )
#define  KRATOS_BREP_APPLICATION_OPENCASCADE_DEFINE_H_INCLUDED

#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>

#include "includes/define.h"

namespace Kratos
{

namespace OCC
{

template<class T>
using shared_ptr = boost::shared_ptr<T>;

template<class T>
using weak_ptr = boost::weak_ptr<T>;

/// Helper function to print OCC object
template<class T>
std::string PrintObject(const T& rObject)
{
    std::stringstream ss;
    BRepTools::Dump(rObject, ss);
    return ss.str();
}

/// Helper function to print OCC object
template<class T>
void DumpObject(std::ostream& rOStream, const T& rObject)
{
    BRepTools::Dump(rObject, rOStream);
}

shared_ptr<TopoDS_Shape> MakeBottle(const Standard_Real myWidth,
    const Standard_Real myHeight, const Standard_Real myThickness);

shared_ptr<TopoDS_Shape> MakeSphere(const Standard_Real cx, const Standard_Real cy,
    const Standard_Real cz, const Standard_Real r);

shared_ptr<TopoDS_Shape> ReadSTEP(const Standard_CString filename);

void WriteSTEP(shared_ptr<TopoDS_Shape> pShape, const Standard_CString filename);

} // namespace OCC

} // namespace Kratos

#endif // KRATOS_BREP_APPLICATION_OPENCASCADE_DEFINE_H_INCLUDED
