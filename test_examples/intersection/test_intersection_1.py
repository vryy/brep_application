import math

from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

xt = MonomialFunctionR1R1X()

yt = ScalarFunctionR1R1(1.0)

zt = ZeroFunctionR1R1()

alignment_curve = ParametricCurve(xt, yt, zt) # curve y=1.0

R = 1.0
ls1 = DistanceToCurveLevelSet(alignment_curve, R)

ls2 = PlanarLevelSet(1.0, 0.0, 0.0, -1.0) # plane x=1

brep_intersection_util = BRepIntersectionUtility()
points = brep_intersection_util.Intersect(ls1, ls2, 10)
print("intersection points:")
for p in points:
    print(p)

