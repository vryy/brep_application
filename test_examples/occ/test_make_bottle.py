from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *

util = OCCUtility()
bottle = util.MakeBottle(70.0, 50.0, 30.0)
print(bottle)

util.WriteSTEP(bottle, "bottle.stp")
