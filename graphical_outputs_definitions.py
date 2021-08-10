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


class plotChoice: # or PC
    def __init__(self):
        self.graph1 = 0
        self.save1 = 0


#---------------------------------------------------------------------
# Basic functions
#---------------------------------------------------------------------
        
def initGraph():        
    plt.rc('font', family='Serif')
    plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
    plt.figure()

# /!\ Needs to be adapted to any dir/file name
def savePDF(path):
    os.chdir(path)
    savefig('pressure_source.pdf', bbox_inches='tight')
    os.chdir("../")

 
#---------------------------------------------------------------------
# Graph 1 : 5 iterations of P(t) in 1 reservoir
#---------------------------------------------------------------------

def plotGraph1(RES, PC):
    
    if PC.graph1 == 1:
        
        #blue set
        colorList = ['black', 'darkblue', 'dodgerblue', 'turquoise', 'mediumspringgreen']
        #red set
        #colorList = ['black', 'saddlebrown', 'orangered', 'orange', 'gold']
        legendList = ['i = 0', 'i = 1', 'i = 2', 'i = 3', 'i = 4']
        
        t_max = RES.t_c*2 
        initGraph()
        for i in RES.i_val[0:5]:
            t_val = np.linspace(0, t_max, 1000)
            x = [0, RES.t_val[i], RES.t_val[i], RES.t_val[i]]
            y = [0, RES.p_val[i], 0, 0]
            p_val = np.interp(t_val, x, y)
            plt.plot(t_val, p_val/1e6, colorList[i], label = legendList[i])
            
        plt.legend(bbox_to_anchor=(1,0.5), loc="center left")
        #plt.xlim((0,1.5))
        plt.xlim((0,t_max))
        #plt.ylim((0,1.01))
        plt.xlabel(r'Time (s)')
        plt.ylabel('Pressure (MPa)')
        
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results/iterative_pressure")
            
        plt.show()
    







































