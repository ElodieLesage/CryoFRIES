#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
------------------------------------------------------------------------

CryoFRIES

Written by Elodie Lesage for 
"Simulation of freezing cryomagma reservoirs in viscoelastic ice shells"
elodie.lesage@jpl.nasa.gov

(c) 2022 California Institute of Technology. All rights Reserved.

This software models the eruption of a sperical cryomagma reservoir
embeded in the ice shell of Europa.

------------------------------------------------------------------------
"""

import reservoir_definitions as rd
import graphical_outputs_definitions as gd

import matplotlib.pyplot as plt
import numpy as np



"""---------------------------------------------------------------------
#       INPUTS
#--------------------------------------------------------------------"""

# Choose between pure water or briny cryomagma and ice:
# If you take into account impurities such as salts in the cryomagma and 
# ices, select 1. If you consider pure water, select 0.
salts = 0

# reservoir radius (m)
# radius = 1905.46071796
radius = 500

# reservoir depth (m) (< 5500m to stay in the elastic zone)
depth = 2000

"""------------------------------------------------------------------'"""




#---------------------------------------------------------------------------
# Initialize arrays
#---------------------------------------------------------------------------

# Thermal properties class
TP = rd.thermalParameters() 

# Physical properties class
PP = rd.physicalParameters()

# Body properties class
BP = rd.bodyParameters()

# Modell derived values
RES = rd.reservoirModelDerivedValues()
RES.R = radius
RES.h = depth

# Graphical outputs
PC = gd.plotChoice

#---------------------------------------------------------------------------
# Set the cryomagma and ice composition and properties
#---------------------------------------------------------------------------
    
# Constants depending on the liquid/ice composition :
TP, PP = rd.cryoComp(salts,TP,PP)

#---------------------------------------------------------------------------
# Freezing 
#---------------------------------------------------------------------------

RES = rd.cryoFreeze(BP,PP,TP,RES)
print('t_c = ', RES.t_c/3600/24, ' days')

#---------------------------------------------------------------------------
# Deformation from Dragoni & Maganensi :
#---------------------------------------------------------------------------

RES = rd.findR2(TP, PP, BP, RES)

RES = rd.cryoRheol(BP,PP,TP,RES)
print('tau = ', RES.tau/3600/24, ' days')

#---------------------------------------------------------------------------
# Stress field and deformation of 1 reservoir:
#---------------------------------------------------------------------------


""" graphical output? 0 = no, 1 = yes """
PC.graph1 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save1 = 0



t = RES.t_c
p = RES.deltaP_c
R1 = RES.R1
R2 = RES.R2
gd.plotGraph1(t, p, R1, R2, RES, PC)


#---------------------------------------------------------------------------
# Iterative model for 1 reservoir:
#---------------------------------------------------------------------------
RES = rd.iterate(PP, TP, RES)

""" graphical output? 0 = no, 1 = yes """
PC.graph2 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save2 = 0

gd.plotGraph2(RES, PC)





