#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 16:46:38 2021
@author: lesage
"""

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import savefig
import os


import reservoir_definitions as rd
RES = rd.reservoirModelDerivedValues()

class plotChoice: # or PC
    def __init__(self):
        self.graph1 = 0
        self.save1 = 0

class outputParameters: # or OUT
    def __init__(self): 
        # radii list
        self.r_val = 0
        # depths list
        self.h_val = 0
        # list of the freezing times with fixed walls
        self.tcFix = 0
        # list of the freezing times with reservoir deformation
        self.tcDeform = 0
        
#---------------------------------------------------------------------
# Basic functions
#---------------------------------------------------------------------
        
def initGraph():        
    plt.rc('font', family='Serif')
    plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
    plt.figure()

def savePDF(path, name, RES):
    os.chdir(path)
    newpath = path + "/R_" + str(RES.R) + "_H_" + str(RES.h)
    isDir = os.path.isdir(newpath)
    if isDir == False:
        os.mkdir(newpath)
    os.chdir(newpath)
    savefig(name, bbox_inches='tight')
    os.chdir("../../")


#---------------------------------------------------------------------
# Graph 1 : stress and displacement fields at varying r and t
# inputs: time t over which the pressure p is applied
#---------------------------------------------------------------------

def plotGraph1(t_c, p_c, R1, R2, RES, PC):
    
    if PC.graph1 == 1:
        
        r, t, t0, t1, t2, t3, p0 = sym.symbols('r t t0 t1 t2 t3 p0')
        
        r1_val = np.linspace(R1, R2, 500)
        r2_val = np.linspace(R2, R2*3, 500) 
        
        # times used in the stress/deformation field equations (= time scale
        # of mechanical constrain)
        t1_temp = 1*t_c
        t2_temp = 1.1*t_c
        t3_temp = 1.2*t_c
        
        u1_init = RES.u1.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        u2_init = RES.u2.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        sigma1_rr_init = RES.sigma1_rr.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        sigma2_rr_init = RES.sigma2_rr.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        sigma1_tt_init = RES.sigma1_tt.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        sigma2_tt_init = RES.sigma2_tt.subs(p0,p_c).subs(t1,t1_temp).subs(t2,t2_temp).subs(t3,t3_temp)
        
        # times at which we sample the stress/deformation fields for
        # graphical outputs
        tFix_0 = 0
        tFix_1 = t_c
        tFix_2 = t_c*2
        tFix_inf = t_c*100

        time_list = [tFix_0, tFix_1, tFix_2, tFix_inf]
     
        # dictionnaries to store the valutes at each time
        u1 = {}
        u2 = {}
        sigma1_rr = {}
        sigma2_rr = {}
        sigma1_tt = {}
        sigma2_tt = {}
        
        u1_val = {}
        u2_val = {}
        sigma1_rr_val = {}
        sigma2_rr_val = {}
        sigma1_tt_val = {}
        sigma2_tt_val = {}
        
        for time in time_list :
            u1[time] = u1_init.subs(t,time)
            u2[time] = u2_init.subs(t,time)
            sigma1_rr[time] = sigma1_rr_init.subs(t,time)
            sigma2_rr[time] = sigma2_rr_init.subs(t,time)
            sigma1_tt[time] = sigma1_tt_init.subs(t,time)
            sigma2_tt[time] = sigma2_tt_init.subs(t,time)
            
            u1_val[time] = np.zeros(len(r1_val))
            u2_val[time] = np.zeros(len(r2_val))
            sigma1_rr_val[time] = np.zeros(len(r1_val))
            sigma2_rr_val[time] = np.zeros(len(r2_val))
            sigma1_tt_val[time] = np.zeros(len(r1_val))
            sigma2_tt_val[time] = np.zeros(len(r2_val))
            
            for i in range(len(r1_val)) :
                u1_val[time][i] = u1[time].subs(r,r1_val[i])
                sigma1_rr_val[time][i] = sigma1_rr[time].subs(r,r1_val[i])
                sigma1_tt_val[time][i] = sigma1_tt[time].subs(r,r1_val[i])
            for i in range(len(r2_val)) :
                u2_val[time][i] = u2[time].subs(r,r2_val[i])
                sigma2_rr_val[time][i] = sigma2_rr[time].subs(r,r2_val[i])
                sigma2_tt_val[time][i] = sigma2_tt[time].subs(r,r2_val[i])
       
        #blue set
        color_list = ['k:', 'darkblue', 'dodgerblue', 'turquoise']
        #red set
        #color_list = ['k:', 'saddlebrown', 'orangered', 'orange']
        legend_list = ['t = 0', r't = $\tau_c$', r't = 2$\tau_c$', r't = 100$\tau_c$']
        
        # Graph of u(r,t)
        initGraph()
        for (time, color, legend) in zip(time_list,color_list, legend_list) :
            plt.plot(r1_val/R1, u1_val[time], color, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, u2_val[time], color)
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'k:')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'k:')
        plt.legend()
        plt.xlabel('$r/R_1$')
        plt.ylabel('$u(r)$ (m)')
        plt.legend(loc='upper right')
        #plt.xlim((1,4))
        #plt.ylim((0,60))
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "u_r_t", RES)
        plt.show()
         
        # Graph of sigma_rr(r,t)
        initGraph()
        for (time, color, legend) in zip(time_list,color_list, legend_list) :
            plt.plot(r1_val/R1, abs(sigma1_rr_val[time])/p_c, color, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, abs(sigma2_rr_val[time])/p_c, color)
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'k:')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'k:')
        plt.legend()
        plt.xlabel('$r/R_1$')
        plt.ylabel(r'$\sigma_{rr}(r)/\Delta P_c$')
        plt.legend(loc='upper right')
        #plt.xlim((1,4))
        #plt.ylim((0,1.01))
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "sigma_rr_r_t", RES)
        plt.show()
       
        # Graph of sigma_tt(r,t) = sigma_pp(r,t)
        initGraph()
        for (time, color, legend) in zip(time_list,color_list, legend_list) :
            plt.plot(r1_val/R1, (sigma1_tt_val[time])/p_c, color, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, (sigma2_tt_val[time])/p_c, color)
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'k:')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'k:')
        plt.legend()
        plt.xlabel('$r/R_1$')
        plt.ylabel(r'$\sigma_{\theta\theta}(r)/\Delta P_c$ = $\sigma_{\phi\phi}(r)/\Delta P_c$')
        plt.legend(loc='upper right')
        #plt.xlim((1,4))
        #plt.ylim((-1.01,0.55))
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "sigma_tt_r_t", RES)     


 
#---------------------------------------------------------------------
# Graph 2 : 5 iterations of P(t) in 1 reservoir
#---------------------------------------------------------------------

def plotGraph2(RES, PC):
    
    if PC.graph2 == 1:
        
        #blue set
        colorList = ['black', 'darkblue', 'dodgerblue', 'turquoise', 'mediumspringgreen']
        #red set
        #colorList = ['black', 'saddlebrown', 'orangered', 'orange', 'gold']
        legendList = ['i = 0', 'i = 1', 'i = 2', 'i = 3', 'i = 4']
        
        t_max = RES.t_val[-1]*1.1 
        initGraph()
        for i in RES.i_val[0:5]:
            t_val = np.linspace(0, t_max, 1000)
            x = [0, RES.t_val[i], RES.t_val[i], RES.t_val[i]]
            y = [0, RES.p_val[i], 0, 0]
            p_val = np.interp(t_val, x, y)
            plt.plot(t_val, p_val/1e6, colorList[i], label = legendList[i])
        
        t_val = np.linspace(0, t_max, 1000)
        x = [0, RES.t_val[-1], RES.t_val[-1], RES.t_val[-1]]
        y = [0, RES.p_val[-1], 0, 0]
        p_val = np.interp(t_val, x, y)
        plt.plot(t_val, p_val/1e6, 'k--', label = 'Last iteration')
            
        plt.legend(bbox_to_anchor=(1,0.5), loc="center left")
        #plt.xlim((0,1.5))
        plt.xlim((0,t_max))
        #plt.ylim((0,1.01))
        plt.xlabel(r'Time (s)')
        plt.ylabel('Pressure (MPa)')
        
        if PC.save2 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "iterative_model", RES)
            
        plt.show()
        
        
#---------------------------------------------------------------------
# Graph 3: Loop on several radii and depths
#---------------------------------------------------------------------

def plotGraph3(OUT):
    
    # axis:
    x = OUT.r_val 
    y = OUT.h_val/1000
    
    tcFixYears = OUT.tcFix /3600 /24 /365.25
    tcDeformYears = OUT.tcDeform /3600 /24 /365.25
    
    plt.rc('font', family='Serif')
    plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
    plt.figure(figsize=(15,7))
    plt.subplot(1,2,1)
    
    ax = plt.subplot(1,2,1)
    ax.set_xscale('log')
    background = plt.pcolor(x, y, tcFixYears, cmap='RdYlBu', vmin = 1e-4, vmax = 1e4, norm=LogNorm())
    contour_dashed = plt.contour(x, y, tcFixYears, norm=LogNorm(), colors='black', linestyles='dashed')
    plt.clabel(contour_dashed, inline=True, fontsize=14, fmt= '%.2f')
    cbar=plt.colorbar(background)
    cbar.set_label(r"$t_{c}$ (years)", labelpad=-40, y=1.05, rotation=0)
    plt.title(u"(a) Freezing time""\n""with fix wall")
    plt.xlabel(u"Reservoir radius (m)")
    plt.ylabel("Reservoir depth (km)")
    ttl = ax.title
    ttl.set_position([.5, 1.02])
    
    plt.subplot(1,2,2)
    ax = plt.subplot(1,2,2)
    ax.set_xscale('log')
    #plt.gca().invert_yaxis()
    background = plt.pcolor(x, y, tcDeformYears, cmap='RdYlBu', vmin = 1e-4, vmax = 1e4, norm=LogNorm())
    contour_dashed = plt.contour(x, y, tcDeformYears, norm=LogNorm(), colors='black', linestyles='dashed')
    plt.clabel(contour_dashed, inline=True, fontsize=14, fmt= '%.2f')
    cbar=plt.colorbar(background)
    cbar.set_label(r"$t_{cv}$ (years)", labelpad=-40, y=1.05, rotation=0)
    plt.title(u"(b) Freezing time""\n""with deformation")
    plt.xlabel(u"Reservoir radius (m)")
    plt.ylabel("Reservoir depth (km)")
    ttl = ax.title
    ttl.set_position([.5, 1.02])

    plt.show()




































