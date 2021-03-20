##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 14. Mar 00:15:32 CET 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./cube1.gid')
import cube1_include
from cube1_include import *
# calculate insitu-stress for geology_virgin.gid
model = cube1_include.Model('cube1',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

# =====================
# | USER SCRIPT FOR CALCULATION OF EKATE.GID |
# vvvvvvvvvvvvvvvvvvvvv
ls = PlanarLevelSet(0.0, 0.0, 1.0, -0.5)
sec = ls.Intersect(model.model_part.Elements[1])
print("section: " + str(sec))

tria = sec.Triangulation()
print("points:")
for p in tria[0]:
    print(p)
print("triangulation:", tria[1])

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
