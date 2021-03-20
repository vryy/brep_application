import math

from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

xt = MonomialFunctionR1R1X()

yt = SinFunctionR1R1(xt)

zt = ZeroFunctionR1R1()

alignment_curve = ParametricCurve(xt, yt, zt) # curve y=sin(x)

ls = PlanarLevelSet(1.0, 0.0, 0.0, -1.0) # plane x=1

P = alignment_curve.ComputeIntersection(ls)
print("intersection: " + str(P))
print("correct intersection: 1.0, " + str(math.sin(1.0)) + ", 0.0")
