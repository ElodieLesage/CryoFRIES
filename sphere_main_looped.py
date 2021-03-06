#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
------------------------------------------------------------------------

CryoFRIES

Written by Elodie Lesage for 
"Simulation of freezing cryomagma reservoirs in viscoelastic ice shells"
elodie.lesage@jpl.nasa.gov

(c) 2022 California Institute of Technology. All rights Reserved.

This software models the eruption of sperical cryomagma reservoirs
with a range of radius and depth, embeded in the ice shell of Europa.

------------------------------------------------------------------------
"""

import reservoir_definitions as rd 
import graphical_outputs_definitions as gd
import numpy as np



"""---------------------------------------------------------------------
#       INPUTS
#--------------------------------------------------------------------"""

# Choose between pure water or briny cryomagma and ice:
# If you take into account impurities such as salts in the cryomagma and 
# ices, select 1. If you consider pure water, select 0.
salts = 0

#  reservoir radii (m)
r_val = np.logspace(1, 4, 61)

# reservoir depth (m) (< 5500m to stay in the elastic zone)
h_val = np.linspace(1000, 10000, 61)

"""------------------------------------------------------------------'"""

#---------------------------------------------------------------------
#       OUTPUTS
#---------------------------------------------------------------------

OUT = rd.outputParameters()

OUT.r_val = r_val
OUT.h_val = h_val
OUT.tcFix = np.zeros((len(h_val), len(r_val)))
OUT.tcFixFilter = np.zeros((len(h_val), len(r_val)))
OUT.tcDeform = np.zeros((len(h_val), len(r_val)))
OUT.tcDeformFilter = np.zeros((len(h_val), len(r_val)))
OUT.VeFix = np.zeros((len(h_val), len(r_val)))
OUT.VeFixFilter = np.zeros((len(h_val), len(r_val)))
OUT.VeDeform = np.zeros((len(h_val), len(r_val)))
OUT.VeDeformFilter = np.zeros((len(h_val), len(r_val)))


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
        
        #---------------------------------------------------------------------------
        # Outputs:
        #---------------------------------------------------------------------------        
        
        OUT.tcFix[i,j] = RES.t_val[0]
        OUT.VeFix[i,j] = RES.V_i*(1-(1-RES.n*(PP.rho_w/PP.rho_i)))

        if RES.isConverging == 1:
            OUT.tcDeform[i,j] = RES.t_val[-1]
            OUT.VeDeform[i,j] = RES.V_i*(1-(1-RES.nFinal*(PP.rho_w/PP.rho_i)))            
        else:
            OUT.tcDeform[i,j] = np.nan
            OUT.VeDeform[i,j] = np.nan
            
        if RES.isMaxwell == 1:
            OUT.tcDeformFilter[i,j] = RES.t_val[-1]
            OUT.VeDeformFilter[i,j] = RES.V_i*(1-(1-RES.nFinal*(PP.rho_w/PP.rho_i)))
            OUT.tcFixFilter[i,j] = OUT.tcFix[i,j]
            OUT.VeFixFilter[i,j] = OUT.VeFix[i,j]
        else:
            OUT.tcDeformFilter[i,j] = np.nan
            OUT.VeDeformFilter[i,j] = np.nan
            OUT.tcFixFilter[i,j] = np.nan
            OUT.VeFixFilter[i,j] = np.nan         
                
                

         
        j = j+1
        
    i = i+1


#---------------------------------------------------------------------------
# Graph of the freezing time (fix Vs viscoelastic)
#---------------------------------------------------------------------------

""" graphical output? 0 = no, 1 = yes """
PC.graph3 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save3 = 1

gd.plotGraph3(OUT, PC)


#---------------------------------------------------------------------------
# Graph of the erupted volume (fix Vs viscoelastic)
#---------------------------------------------------------------------------

""" graphical output? 0 = no, 1 = yes """
PC.graph4 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save4 = 1

gd.plotGraph4(OUT, PC)

#---------------------------------------------------------------------------
# Graph of the freezing time (fix Vs viscoelastic) for a res. of given radius
#---------------------------------------------------------------------------

""" graphical output? 0 = no, 1 = yes """
PC.graph5 = 1
""" save as PDF? 0 = no, 1 = yes """
PC.save5 = 0

gd.plotGraph5(OUT, PC)
















