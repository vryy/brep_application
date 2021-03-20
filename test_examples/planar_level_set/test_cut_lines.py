from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

P1 = Point3D(1.0, 1.0, 2.0)
P2 = Point3D(0.0, 0.0, -1.0)

ls = PlanarLevelSet(0.0, 0.0, 1.0, 0.0)

P = ls.Bisect(P1, P2)
print("level set: " + str(ls))
print("intersection point 1: " + str(P))

P3 = Point3D(0.0, 0.0, 1.0)

P = ls.Bisect(P1, P2, P3)
print("intersection point 2: " + str(P))
