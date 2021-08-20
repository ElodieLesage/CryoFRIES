#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:40:25 2020
@author: elodielesage
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

# reservoir radii (m)
r_val = 10** np.linspace(1, 3, 31)

# reservoir depth (m) (< 5500m to stay in the elastic zone)
h_val = np.linspace(1000, 10000, 51)

"""------------------------------------------------------------------'"""

#---------------------------------------------------------------------
#       OUTPUTS
#---------------------------------------------------------------------

OUT = gd.outputParameters()

OUT.r_val = r_val
OUT.h_val = h_val
OUT.tcFix = np.zeros((len(h_val), len(r_val)))
OUT.tcDeform = np.zeros((len(h_val), len(r_val)))



#---------------------------------------------------------------------
#       MAIN LOOP
#---------------------------------------------------------------------

i = 0

for depth in h_val:
    
    j = 0
    
    for radius in r_val:
        
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
        #print('t_c = ', RES.t_c/3600/24, ' days')
        
        #---------------------------------------------------------------------------
        # Deformation from Dragoni & Maganensi :
        #---------------------------------------------------------------------------
        
        RES = rd.findR2(TP, PP, BP, RES)
        RES = rd.cryoRheol(BP,PP,TP,RES)
        
        #---------------------------------------------------------------------------
        # Iterative model for 1 reservoir:
        #---------------------------------------------------------------------------
        print('Depth: ', depth, ' / Radius: ', radius)
        RES = rd.iterate(PP, TP, RES)
        
        OUT.tcFix[i,j] = RES.t_val[0]
        
        if RES.isConverging == 1:
            OUT.tcDeform[i,j] = RES.t_val[-1]
        else:
            OUT.tcDeform[i,j] = float("NAN")
        
        j = j+1
        
    i = i+1

gd.plotGraph3(OUT)





















