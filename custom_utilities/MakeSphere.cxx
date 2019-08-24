
#include <gp_Ax2.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>

#include "includes/define.h"
#include "custom_utilities/occ_define.h"

namespace Kratos
{

namespace OCC
{

shared_ptr<TopoDS_Shape> MakeSphere(const Standard_Real cx, const Standard_Real cy,
    const Standard_Real cz, const Standard_Real r)
{
    Standard_Real sphere_radius = r;
    // Standard_Real sphere_angle = atan(0.5);
    Standard_Real sphere_angle = atan(1.0)*2;

    gp_Ax2 sphere_origin = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(cx, cy, cz+1.0));
    shared_ptr<TopoDS_Shape> sphere = shared_ptr<TopoDS_Shape>(new TopoDS_Shape());
    *sphere = BRepPrimAPI_MakeSphere(sphere_origin, sphere_radius, -sphere_angle, sphere_angle).Shape();
    // *sphere = BRepPrimAPI_MakeSphere(sphere_origin, sphere_radius, sphere_angle).Shape();
    // *sphere = BRepPrimAPI_MakeSphere(sphere_origin, sphere_radius, sphere_angle).Shape();

    return sphere;
}

} // namespace OCC

} // namespace Kratos
