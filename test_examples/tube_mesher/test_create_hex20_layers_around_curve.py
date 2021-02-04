import math
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.StructuralApplication import *

def WriteOutput( model_part, time ):
    write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
    write_elements = WriteConditionsFlag.WriteConditions
    #write_elements = WriteConditionsFlag.WriteElementsOnly
    post_mode = GiDPostMode.GiD_PostBinary
    multi_file_flag = MultiFileFlag.MultipleFiles
    gid_io = GidIO( model_part.Name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
    gid_io.InitializeMesh( time )
    mesh = model_part.GetMesh()
    print("mesh at time " + str(time) + " is ready for printing")
    gid_io.WriteMesh( mesh )
    print("mesh written...")
    gid_io.FinalizeMesh()

t_list = [0.0, 1.0]
x_list = [0.0, 10.0]
y_list = [0.0, 0.0]
z_list = [0.0, 0.0]

xt = CubicSplineFunctionR1R1()
xt.SetLeftBoundary(1, (x_list[1]-x_list[0]) / (t_list[1]-t_list[0]))
xt.SetRightBoundary(1, (x_list[len(x_list)-1]-x_list[len(x_list)-2]) / (t_list[len(t_list)-1]-t_list[len(t_list)-2]))
xt.SetPoints(t_list, x_list)

yt = CubicSplineFunctionR1R1()
yt.SetLeftBoundary(1, (y_list[1]-y_list[0]) / (t_list[1]-t_list[0]))
yt.SetRightBoundary(1, (y_list[len(y_list)-1]-y_list[len(y_list)-2]) / (t_list[len(t_list)-1]-t_list[len(t_list)-2]))
yt.SetPoints(t_list, y_list)

zt = CubicSplineFunctionR1R1()
zt.SetLeftBoundary(1, (z_list[1]-z_list[0]) / (t_list[1]-t_list[0]))
zt.SetRightBoundary(1, (z_list[len(z_list)-1]-z_list[len(z_list)-2]) / (t_list[len(t_list)-1]-t_list[len(t_list)-2]))
zt.SetPoints(t_list, z_list)

alignment_curve = ParametricCurve(xt, yt, zt)
# alignment_curve.Export("alignment.txt", 0.0, 1.0, 100, 2);
alignment_curve.SetValue(CURVE_LOWER_BOUND, 0.0)
alignment_curve.SetValue(CURVE_UPPER_BOUND, 1.0)
alignment_curve.SetValue(CURVE_NUMBER_OF_SAMPLING, 20)
print("curve is created")

mesher = BRepMeshUtility()
r_list = [1.0, 2.0, 3.0, 4.0]
nsampling_layers = [1, 2, 3]
nsampling_axial = 5
nsampling_radial = 10
tmin = 0.0
tmax = 1.0
thetype = 2 # H20
rotate_angle = 0.0
start_angle = 0.0
end_angle = 2.0*math.pi
last_node_id = 0
mesher = TubeMesher(alignment_curve, r_list, nsampling_layers, nsampling_axial, nsampling_radial, rotate_angle, start_angle, end_angle, tmin, tmax, thetype, last_node_id)
points = mesher.GetPoints()
# print("points: ", len(points))
elements = mesher.GetElements()
# print("elements:", elements)
conditions = mesher.GetConditions()
# print("conditions:", conditions)
# slice000 = mesher.GetSlice(0, 0, 0) # slice, layer, sub-layer
# print("slice000:", slice000)
slice00 = mesher.GetSlice(0, 0) # slice, layer
slice01 = mesher.GetSlice(0, 1) # slice, layer
slice02 = mesher.GetSlice(0, 2) # slice, layer
slice10 = mesher.GetSlice(1, 0) # slice, layer
slice11 = mesher.GetSlice(1, 1) # slice, layer
slice12 = mesher.GetSlice(1, 2) # slice, layer
del mesher
# print("points:")
# for p in results[0]:
#     print(str(p))
print("----------------------")
# print("elements:", elements)
print("elemental info:")
print("num layers: " + str(len(elements)))
print("num sub-layers:")
for i in range(0, len(elements)):
    print(len(elements[i]))
print("num rings: " + str(len(elements[0][0])))
print("num segments: " + str(len(elements[0][0][0])))
print("num nodes: " + str(len(elements[0][0][0][0])))
print("----------------------")
print("conditional info:")
# print("conditions:", conditions)
print("num boundary layers: " + str(len(conditions)))
print("num rings: " + str(len(conditions[0])))
print("num segments: " + str(len(conditions[0][0])))
print("num nodes: " + str(len(conditions[0][0][0])))
print("----------------------")

mp = ModelPart("mesh_h20")

idx = 1
for p in points:
    mp.CreateNewNode(idx, p[0], p[1], p[2])
    idx = idx + 1

element_name = "KinematicLinear3D20N"
idx = 1
layer = 1
for tmp1 in elements:
    elem_prop = mp.Properties[layer]
    for tmp2 in tmp1:
        for tmp3 in tmp2:
            for tmp4 in tmp3:
                # print(tmp4)
                mp.CreateNewElement(element_name, idx, tmp4, elem_prop)
                idx = idx + 1
    layer = layer + 1

# elem_prop = mp.Properties[2]
# for tmp3 in elements[1][0]:
#     for tmp4 in tmp3:
#         print(tmp4)
#         mp.CreateNewElement(element_name, idx, tmp4, elem_prop)
#         idx = idx + 1
# print("------------------------")
# elem_prop = mp.Properties[3]
# for tmp3 in elements[1][1]:
#     for tmp4 in tmp3:
#         print(tmp4)
#         mp.CreateNewElement(element_name, idx, tmp4, elem_prop)
#         idx = idx + 1

condition_name = "FaceForce3D8N"
layer = 10
for tmp1 in conditions:
    cond_prop = mp.Properties[layer]
    for tmp2 in tmp1:
        for tmp3 in tmp2:
            # print(tmp3)
            mp.CreateNewCondition(condition_name, idx, tmp3, cond_prop)
            idx = idx + 1
    layer = layer + 1

condition_name = "FaceForce3D8N"
cond_prop = mp.Properties[20]
for tmp1 in slice00:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

cond_prop = mp.Properties[21]
for tmp1 in slice01:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

cond_prop = mp.Properties[22]
for tmp1 in slice02:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

cond_prop = mp.Properties[30]
for tmp1 in slice10:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

cond_prop = mp.Properties[31]
for tmp1 in slice11:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

cond_prop = mp.Properties[32]
for tmp1 in slice12:
    for tmp2 in tmp1:
        # print(tmp2)
        mp.CreateNewCondition(condition_name, idx, tmp2, cond_prop)
        idx = idx + 1

print(mp)

# for node in mp.Nodes:
#     print(str(node.X0) + "," + str(node.Y0) + "," + str(node.Z0))

WriteOutput(mp, 0.0)

