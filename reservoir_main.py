#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:40:25 2020

@author: elodielesage

"""

import reservoir_definitions as rd
import graphical_outputs_definitions as gd

import numpy as np
from numpy import linspace
from numpy import zeros
import matplotlib.pyplot as plt
from pylab import savefig
import os



"""---------------------------------------------------------------------
#       INPUTS
#--------------------------------------------------------------------"""

# Choose between pure water or briny cryomagma and ice:
# If you take into account impurities such as salts in the cryomagma and 
# ices, select 1. If you consider pure water, select 0.
salts = 0

# reservoir radius (m)
radius = 500

# reservoir depth (m) (< 5500m to stay in the elastic zone)
depth = 2000

"""------------------------------------------------------------------'"""




#---------------------------------------------------------------------
# Initialize arrays
#---------------------------------------------------------------------

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

#---------------------------------------------------------------------
# Deformation from Dragoni & Maganensi :
#---------------------------------------------------------------------

RES = rd.cryoRheol(BP,PP,TP,RES)
print('tau = ', RES.tau/3600/24, ' days')

"""
#---------------------------------------------------------------------
# Results for 1 reservoir, without graphical output
#---------------------------------------------------------------------

epsilon = RES.deltaP_c/100000
i=0

dP = rd.time2deformation(RES.t_c, PP, RES)
dP_temp = dP

dP = rd.pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)

while dP_temp - dP > epsilon :
    dP_temp = dP
    dP = rd.pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)
    print(RES.deltaP_c + dP_temp - dP)
    i = i+1

print('--------------------------------------------')
print('number of iterations: ', i)

"""


#---------------------------------------------------------------------
# Results for 1 reservoir, with graphical output
#---------------------------------------------------------------------


epsilon = RES.deltaP_c/100000
t_max = RES.t_c*2
i=0

PL = gd.plotLegend()

gd.initGraph()

gd.plotPressure(i, RES.t_c, RES.deltaP_c, t_max, PL)

dP = rd.time2deformation(RES.t_c, PP, RES)
dP_temp = dP

[dP, t_c] = rd.pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)
while dP_temp - dP > epsilon :
    i = i+1
    if i<5:
        gd.plotPressure(i, t_c, RES.deltaP_c - dP_temp, t_max, PL)
    dP_temp = dP
    [dP, t_c] = rd.pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)
    

print('--------------------------------------------')
print('number of iterations: ', i)

plt.legend(bbox_to_anchor=(1,0.5), loc="center left")
#plt.xlim((0,1.5))
plt.xlim((0,t_max))
#plt.ylim((0,1.01))
plt.xlabel(r'Time (s)')
plt.ylabel('Pressure (MPa)')

#gd.savePDF()
plt.show()


"""

#baisse_pression = deform_pressure(deltaP_c + 2.2971e6)
# For a 500m in radius reservoir, with h = 2000, R2 = 1.3R1, T = 200 K, 
# deltaP_c converges to deltaP_c+2.2971 MPa approximately, which corresponds to 
# t_c = 1715 days (4.7 years) instead of 1051 days without deformation.
# -> to compare with the Maxwell time

"""