
#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>

#include "includes/define.h"
#include "custom_utilities/occ_define.h"

namespace Kratos
{

namespace OCC
{

shared_ptr<TopoDS_Shape> ReadSTEP(const Standard_CString filename)
{
    STEPControl_Reader reader;
    reader.ReadFile(filename);

    // Load STEP file
    Standard_Integer NbRoots = reader.NbRootsForTransfer();

    // gets the number of transferable roots
    std::cout << "Number of roots in STEP file: " << NbRoots << std::endl;

    Standard_Integer NbTrans = reader.TransferRoots();

    // translates all transferable roots, and returns the number of    //successful translations
    std::cout << "STEP roots transferred: " << NbTrans << std::endl;
    std::cout << "Number of resulting shapes is: " << reader.NbShapes() << std::endl;

    shared_ptr<TopoDS_Shape> result = shared_ptr<TopoDS_Shape>(new TopoDS_Shape());
    *result = reader.OneShape();

    return result;
}

} // namespace OCC

} // namespace Kratos
