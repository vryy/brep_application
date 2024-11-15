import math

from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

def norm(vec):
    return math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

xt = MonomialFunctionR1R1X()

yt = SinFunctionR1R1(xt)

zt = ZeroFunctionR1R1()

alignment_curve = ParametricCurve(xt, yt, zt) # curve y=sin(x)

ls = PlanarLevelSet(1.0, 0.0, 0.0, -1.0) # plane x=1

# test intersection with level set
P = alignment_curve.ComputeIntersection(ls)
print("intersection: " + str(P))
print("correct intersection: 1.0, " + str(math.sin(1.0)) + ", 0.0")
assert(abs(P[1] - math.sin(1.0)) < 1e-10)

# test compute length
alignment_curve[CURVE_NUMBER_OF_SAMPLING] = 20
leng = alignment_curve.ComputeLength(0.0, 1.0)
print("length: " + str(leng))
assert(abs(leng - 1.311442498215547) < 1e-8)

# test equally distance division to two
tolerance = alignment_curve[CURVE_SEARCH_TOLERANCE]
tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 2)
print(tvec)
P1 = alignment_curve.GetValue(tvec[0])
Pm = alignment_curve.GetValue(tvec[1])
P2 = alignment_curve.GetValue(tvec[2])
d1 = P1 - Pm
d2 = P2 - Pm
assert(abs(norm(d1) - norm(d2)) < tolerance)

# test equally distance division to three
tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 3)
print(tvec)
P1 = alignment_curve.GetValue(tvec[0])
Pm1 = alignment_curve.GetValue(tvec[1])
Pm2 = alignment_curve.GetValue(tvec[2])
P2 = alignment_curve.GetValue(tvec[3])
d1 = P1 - Pm1
d2 = Pm2 - Pm1
d3 = P2 - Pm2
assert(abs(norm(d1) - norm(d2)) < tolerance)
assert(abs(norm(d2) - norm(d3)) < tolerance)

# # test equally distance division to ten
# tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 4)
# print(tvec)
# tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 5)
# print(tvec)
# tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 6)
# print(tvec)
# # tvec = alignment_curve.ComputeEquallyDistanceDivisionByRecursiveBisection(0.0, 1.0, 7)
# # print(tvec)

# test uniform division to three
tvec = alignment_curve.ComputeUniformDivision(0.0, 1.0, 3)
print(tvec)
# TODO add check values

# test uniform division to ten
tvec = alignment_curve.ComputeUniformDivision(0.0, 1.0, 10)
print(tvec)
# TODO add check values

# test uniform division to hundred
tvec = alignment_curve.ComputeUniformDivision(0.0, 1.0, 100)
print(tvec)
# TODO add check values
