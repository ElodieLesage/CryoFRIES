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

#  reservoir radii (m) (Radius of the sphere of equivalent volume)
r_val = np.logspace(1, 5, 71)

# reservoir depth (m)
h_val = np.linspace(1000, 10000, 61)

# aspect ratio of the ellipsoid (a=b=Fc, F=1 for a sphere):
F = 500

"""------------------------------------------------------------------'"""

#---------------------------------------------------------------------
#       OUTPUTS
#---------------------------------------------------------------------

OUT = rd.outputParameters()

OUT.r_val = r_val
OUT.h_val = h_val
OUT.tcFixFilter = np.zeros((len(h_val), len(r_val)))
OUT.VeFixFilter = np.zeros((len(h_val), len(r_val)))


# Graphical outputs
PC = gd.plotChoice

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
        RES.F = F
        
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

        # Maxwell time:
        RES.t_max = RES.eta / PP.G  
        
        if RES.t_c > RES.t_max :
            RES.isConverging = 0 
        
        #---------------------------------------------------------------------------
        # Outputs:
        #---------------------------------------------------------------------------        

        OUT.tcFixFilter[i,j] = RES.t_c
        OUT.VeFixFilter[i,j] = RES.V_i*(1-(1-RES.n*(PP.rho_w/PP.rho_i)))
        
        if RES.isConverging == 0:
            OUT.tcFixFilter[i,j] = np.nan
            OUT.VeFixFilter[i,j] = np.nan
        
        j = j+1
        
    i = i+1


#---------------------------------------------------------------------------
# Graph of the freezing time + erupted volume + eruption duration (to do)
#---------------------------------------------------------------------------

""" graphical output? 0 = no, 1 = yes """
PC.graph5 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save5 = 1

gd.plotGraph5(OUT, PC, RES)
















