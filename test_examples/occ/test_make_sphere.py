from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

util = OCCUtility()
sphere = util.MakeSphere(0.0, 0.0, 0.0, 1.0)
# print(sphere)
util.WriteSTEP(sphere, "sphere.stp")

brep = OCCBRep()
brep.SetShape(sphere)
print(brep)

P1 = Point3D(0.0, 0.1, 0.2)
print(brep.IsInside(P1)) # shall be True

P2 = Point3D(1.0, 1.0, 1.0)
print(brep.IsInside(P2)) # shall be False
