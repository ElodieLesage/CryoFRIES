#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 16:46:38 2021

@author: lesage

"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import savefig
import os

class plotLegend: # or PL
    def __init__(self): 
        #blue set
        self.colorList = ['black', 'darkblue', 'dodgerblue', 'turquoise', 'mediumspringgreen']
        #red set
        #self.colorList = ['black', 'saddlebrown', 'orangered', 'orange', 'gold']
        self.legendList = ['i = 0', 'i = 1', 'i = 2', 'i = 3', 'i = 4']


#---------------------------------------------------------------------
# Basic functions
#---------------------------------------------------------------------
        
def initGraph():        
    plt.rc('font', family='Serif')
    plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
    plt.figure()

def savePDF():
    os.chdir("/Users/elodie/Documents/2020-2021/reservoir_deformation/numerical_model_results/results_trapeze_pressure_with_deformation")
    savefig('pressure_source.pdf', bbox_inches='tight')
    os.chdir("../")


#---------------------------------------------------------------------
# Graph of the pressure as a function of the time for 1 reservoir, 
# 1 iteration
#---------------------------------------------------------------------

def plotPressure(i, t_c, p, t_max, PL):
    t_val = np.linspace(0, t_max, 1000)
    x = [0, t_c, t_c, t_c]
    y = [0, p, 0, 0]
    p_val = np.interp(t_val, x, y)
    
    color = PL.colorList[i]
    legend = PL.legendList[i]
    
    plt.plot(t_val, p_val/1e6, color, label = legend)
 