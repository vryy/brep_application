import math
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

### Compute cross product
def cross(c, a, b):
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c

### Compute dot product
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

### Normalize a vector
def normalize(a):
    norma = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    a[0] = a[0] / norma
    a[1] = a[1] / norma
    a[2] = a[2] / norma
    return a

### Create a list of Frenet frame along a curve. The Frenet frame is stored as a transformation matrix.
### zvec is a reference vector to compute B at the first sampling point. It shall not be parallel with the tangent vector of the first sampling point.
def GenerateLocalFrenetFrame(curve, num_sampling_points, zvec = [1.0, 0.0, 0.0]):
    trans_list = []
    xi_list = []
    B = Array3()
    for i in range(0, num_sampling_points):
        xi = float(i) / (num_sampling_points-1)
        pnt = xi
        P = curve.GetValue(pnt)
        T = curve.GetDerivative(0, pnt)
        T = normalize(T)

        if i == 0:
            cross(B, zvec, T)
            B = normalize(B)
        else:
            B = B - dot(B, T)*T
            B = normalize(B)

        # print("P: " + str(P))
        # print("T: " + str(T))
        # print("B: " + str(B))

        trans = Transformation(B, T, P)
        # print("trans: " + str(trans))
        trans_list.append(trans)
        xi_list.append(xi)

    return [xi_list, trans_list]
