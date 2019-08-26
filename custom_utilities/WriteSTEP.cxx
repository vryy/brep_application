
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>

#include "includes/define.h"
#include "custom_utilities/occ_define.h"

namespace Kratos
{

namespace OCC
{

void WriteSTEP(shared_ptr<TopoDS_Shape> pShape, const Standard_CString filename)
{
    STEPControl_Writer writer;
    writer.Transfer(*pShape, STEPControl_ManifoldSolidBrep);

    // Translates TopoDS_Shape into manifold_solid_brep entity
    writer.Write(filename);

    // writes the resulting entity in the STEP file
}

} // namespace OCC

} // namespace Kratos
