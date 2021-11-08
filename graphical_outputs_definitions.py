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
from math import floor, log10


import reservoir_definitions as rd
RES = rd.reservoirModelDerivedValues()

class plotChoice: # or PC
    def __init__(self):
        self.graph1 = 0
        self.save1 = 0
        self.graph2 = 0
        self.save2 = 0
        self.graph3 = 0
        self.save3 = 0
        self.graph4 = 0
        self.save4 = 0
        
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
    
def savePDF2(path, name):
    os.chdir(path)
    newpath = path + "_looped"
    isDir = os.path.isdir(newpath)
    if isDir == False:
        os.mkdir(newpath)
    os.chdir(newpath)
    savefig(name, bbox_inches='tight')
    os.chdir("../../")
    
# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)


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
        tFix_1 = t_c
        tFix_2 = t_c*1.11
        tFix_3 = t_c*1.15
        tFix_4 = t_c*1.2

        time_list = [tFix_1, tFix_2, tFix_3, tFix_4]
     
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
       
        color_list = ['k', 'k', 'k', 'k']
        linestylelist = ['solid', 'dashed', 'dotted', (0, (1, 5))]
        
        legend1 = r't = $\tau_c$ = ' + sci_notation(RES.t_c,1) + ' s'
        legend2 = r't = 1.11*$\tau_c$ = ' + sci_notation(1.11*RES.t_c,1) + ' s'
        legend3 = r't = 1.15*$\tau_c$ = ' + sci_notation(1.15*RES.t_c,1) + ' s'
        legend4 = r't = 1.2*$\tau_c$ = ' + sci_notation(1.2*RES.t_c,1) + ' s'
        legend_list = [legend1, legend2, legend3, legend4]
        
        # Graph of u(r,t)
        initGraph()
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'lightgrey', label = 't = 0')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'lightgrey')
        for (time, color, linestyle, legend) in zip(time_list, color_list, linestylelist, legend_list) :
            plt.plot(r1_val/R1, u1_val[time], color, linestyle=linestyle, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, u2_val[time], color, linestyle=linestyle)
        plt.xlabel('$r/R_1$')
        plt.ylabel('$u(r)$ (m)')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #plt.xlim((1,4))
        #plt.ylim((0,60))
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "u_r_t", RES)
        plt.show()
         
        # Graph of sigma_rr(r,t)
        initGraph()
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'lightgrey', label = 't = 0')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'lightgrey')
        for (time, color, linestyle, legend) in zip(time_list, color_list, linestylelist, legend_list) :
            plt.plot(r1_val/R1, abs(sigma1_rr_val[time])/p_c, color, linestyle=linestyle, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, abs(sigma2_rr_val[time])/p_c, color, linestyle=linestyle)
        plt.xlabel('$r/R_1$')
        plt.ylabel(r'$\sigma_{rr}(r)/\Delta P_c$')
        #plt.legend(loc='upper right')
        #plt.xlim((1,4))
        #plt.ylim((0,1.01))
        if PC.save1 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "sigma_rr_r_t", RES)
        plt.show()
       
        # Graph of sigma_tt(r,t) = sigma_pp(r,t)
        initGraph()
        plt.plot(r1_val/R1, np.zeros(len(r1_val)), 'lightgrey', label = 't = 0')
        if RES.isElastic == 1:
            plt.plot(r2_val/R1, np.zeros(len(r2_val)), 'lightgrey')
        for (time, color, linestyle, legend) in zip(time_list, color_list, linestylelist, legend_list) :
            plt.plot(r1_val/R1, (sigma1_tt_val[time])/p_c, color, linestyle=linestyle, label = legend)
            if RES.isElastic == 1:
                plt.plot(r2_val/R1, (sigma2_tt_val[time])/p_c, color, linestyle=linestyle)
        plt.xlabel('$r/R_1$')
        plt.ylabel(r'$\sigma_{\theta\theta}(r)/\Delta P_c$ = $\sigma_{\phi\phi}(r)/\Delta P_c$')
        #plt.legend(loc='upper right')
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
        plt.xlabel('log10(Time (s))')
        plt.ylabel('Pressure (MPa)')
        
        if PC.save2 == 1:
            savePDF("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "iterative_model", RES)
            
        plt.show()
        
        
#---------------------------------------------------------------------
# Graph 3: Freezing time, fix Vs viscoelastic
#          Loop on several radii and depths
#---------------------------------------------------------------------

def plotGraph3(OUT, PC):
       
    if PC.graph3 == 1:
        
        # axis:
        x = OUT.r_val 
        y = OUT.h_val/1000
        
        tcFixYears = OUT.tcFix /3600 /24 /365.25
        tcFixFilterYears = OUT.tcFixFilter /3600 /24 /365.25
        tcDeformYears = OUT.tcDeform /3600 /24 /365.25

        j = 0
        for R in x:
            i = 0
            for H in y*1000:
                if H+2*R > 10000:
                    tcFixYears[i,j] = np.nan
                    tcDeformYears[i,j] = np.nan
                i = i+1
            j = j+1
        
        plt.rc('font', family='Serif')
        plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
        plt.figure(figsize=(25,7))
        
        ax = plt.subplot(1,3,1)
        ax.set_xscale('log')
        background = plt.pcolormesh(x, y, tcFixYears, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e-3, vmax = 1e3),shading='gouraud')
        contour_dashed = plt.contour(x, y, tcFixYears, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt= '%.2f')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$t_{c}$ (years)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(a) Freezing time \n with fix wall")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        ax = plt.subplot(1,3,2)
        ax.set_xscale('log')
        background = plt.pcolormesh(x, y, tcFixFilterYears, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e-3, vmax = 1e3), shading='gouraud')
        contour_dashed = plt.contour(x, y, tcFixFilterYears, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt= '%.2f')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$t_{c}$ (years)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(b) Freezing time with \n fix wall and viscoelastic filter")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        ax = plt.subplot(1,3,3)
        ax.set_xscale('log')
        background = plt.pcolormesh(x, y, tcDeformYears, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e-3, vmax = 1e3), shading='gouraud')
        contour_dashed = plt.contour(x, y, tcDeformYears, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt= '%.2f')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$t_{cv}$ (years)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(c) Freezing time \n with viscoelastic deformation")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        if PC.save3 == 1:
            savePDF2("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "freezing_time_allres")
            
        plt.show()


#---------------------------------------------------------------------
# Graph 4: Erupted volume, fix Vs viscoelastic
#          Loop on several radii and depths
#---------------------------------------------------------------------

def plotGraph4(OUT, PC):
           
    if PC.graph4 == 1:
        
        # axis:
        x = OUT.r_val 
        y = OUT.h_val/1000
        
        VeFix = OUT.VeFix
        VeFixFilter = OUT.VeFixFilter
        VeDeform = OUT.VeDeform
   
        j = 0
        for R in x:
            i = 0
            for H in y*1000:
                if H+2*R > 10000:
                    VeFix[i,j] = np.nan
                    VeFixFilter[i,j] = np.nan
                    VeDeform[i,j] = np.nan
                i = i+1
            j = j+1
        
        plt.rc('font', family='Serif')
        plt.rc('font', **{'serif' : 'Times New Roman', 'family' : 'serif', 'size' : 16})
        plt.figure(figsize=(25,7))
        
        ax = plt.subplot(1,3,1)
        ax.set_xscale('log')
        background = plt.pcolormesh(x, y, VeFix, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e3, vmax = 1e10), shading='gouraud')
        contour_dashed = plt.contour(x, y, VeFix, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt='%.1e')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$V_{e}$ (m$^3$)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(a) Erupted volume \n with fix wall")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        ax = plt.subplot(1,3,2)
        ax.set_xscale('log')
        #plt.gca().invert_yaxis()
        background = plt.pcolormesh(x, y, VeFixFilter, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e3, vmax = 1e10), shading='gouraud')
        contour_dashed = plt.contour(x, y, VeFixFilter, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt='%.1e')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$V_{e}$ (m$^3$)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(b) Erupted volume with \n fix wall and viscoelastic filter")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        
        ax = plt.subplot(1,3,3)
        ax.set_xscale('log')
        #plt.gca().invert_yaxis()
        background = plt.pcolormesh(x, y, VeDeform, cmap='RdYlBu_r', norm=LogNorm(vmin = 1e3, vmax = 1e10), shading='gouraud')
        contour_dashed = plt.contour(x, y, VeDeform, norm=LogNorm(), colors='black', linestyles='dashed')
        plt.clabel(contour_dashed, inline=True, fontsize=14, fmt='%.1e')
        cbar=plt.colorbar(background)
        cbar.set_label(r"$V_{ev}$ (m$^3$)", labelpad=-40, y=1.1, rotation=0)
        plt.title(u"(c) Erupted volume \n with viscoelastic deformation")
        plt.xlabel(u"Reservoir radius (m)")
        plt.ylabel("Reservoir depth (km)")
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.gca().invert_yaxis()
        ax = plt.gca()
        ax.set_facecolor('k')
        
        if PC.save4 == 1:
            savePDF2("/Users/lesage/Documents/2020-2021/reservoir_deformation/numerical_model/results", "erupted_volume_allres")
            
        plt.show()



















