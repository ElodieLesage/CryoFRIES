# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 09:21:32 2021

@author: elodielesage

"""


import numpy as np
import scipy as sp
import sympy as sym
import scipy.optimize
from sympy.abc import s,t,r
from sympy.integrals import inverse_laplace_transform
from sympy.integrals import laplace_transform



class thermalParameters: # or TP
    def __init__(self): 
        # melting temperature
        self.T_m = 0
        # ice specific heat capacity
        self.cp_s = 2e3
        # ice thermal conductivity
        self.k_ice = 2.3
        # ice thermal diffusivity
        self.kappa = 0
        # fusion latent heat
        self.Lh = 3e5
        # temperatures at the top and bottom of the ice shell
        self.T_min = 100
        self.T_max = 250

class physicalParameters: # or PP
    def __init__(self): 
        # water compressibility
        self.beta = 5e-10
        # water kinematic viscosity
        self.nu = 0.325
        # h_min = surface, h_max = ocean top
        self.h_min = 0
        self.h_max = 1e4
        # ice tensile strenght at top and bottom of the ice shell
        self.sigma_top = 1.7e6
        self.sigma_bot = 1e6
        # Rigidities :
        self.mu1 = 9e9
        self.mu2 = 9e9
        # Young's moduus
        self.E = 9e9
        # ice density	
        self.rho_i = 920.0
        # water density (at rest)
        self.rho_w = 1000.0 

        
class bodyParameters: # or BP
    def __init__(self): 
        # gravity
        self.g = 1.315
        
class reservoirModelDerivedValues: # or RES
    def __init__(self): 
        # reservoir depth (m) (< 5500m to stay in the elastic zone)
        self.h = 0
        # reservoir radius (m)
        self.R = 0
        # initial reservoir volume
        self.V_i = 0
        # Temperature far from the reservoir
        self.T = 0
        # lithostatique pressure around the reservoir
        self.P0 = 0
        #---------------------------------
        # Lesage et al.'s freezing model :
        #---------------------------------
        # ice tensile strenght in the crust
        self.sigma_c = 0
        # overpressure required to break the reservoir wall 
        self.deltaP_c = 0
        # cryomagma frozen fraction
        self.n = 0
        # thickness of the critical frozen layer
        self.Sc = 0
        # Stefan problem solving
        self.Lambda = 0
        # time required to freeze nc :
        self.t_c = 0
        #---------------------------------
        # Dragoni and Magnanensi's model :
        #---------------------------------
        # Radii of the reservoir and the viscoelastic shell (m)
        self.R1 = 0
        self.R2 = 0
        # Temperature at the edge of the viscoelastic shell
        self.T2 = 0
        # Compressibilities
        self.K1 = 0
        self.K2 = 0
        # Viscosity (only for the viscoelastic zone)
        self.eta = 0
        # viscoelastic timescale
        self.tau = 0
        # outputs 
        self.u1 = 0
        self.u2 = 0
        self.sigma1_rr = 0
        self.sigma2_rr = 0 
        self.sigma1_tt = 0
        self.sigma2_tt = 0
        self.sigma1_pp = 0
        self.sigma2_pp = 0
        #---------------------------------
        # Iterative model
        #---------------------------------
        self.i_val=[]
        self.t_val=[]
        self.p_val=[]
        self.converging = 1 # 0 if the iterative model converges 
        # toward a pressure value, i.e. if the reservoir behaves
        # elastically ; 0 if not, i.e. the reservoir behaves
        # as a viscous material



#---------------------------------------------------------------------
# Basic maths functions
#---------------------------------------------------------------------

# Laplace transform
def L(f):
    return laplace_transform(f, t, s, noconds=True)

# Inverse Laplace transform
def invL(F):
    return inverse_laplace_transform(F, s, t)

# Heaviside function 
def H(x):
    if x != 0:
        return sym.Heaviside(x)
    else:
        return 1
    
# Helpful to debug :
def replaceHeaviside(string):
    string = str(string)
    new_string = string.replace("Heaviside", "H")
    new_string = eval(new_string)
    return new_string

        
#---------------------------------------------------------------------
# Cryomagma composition
#---------------------------------------------------------------------

def cryoComp(salts,TP,PP):
   
    if salts == 0:
        # ice density	
        PP.rho_i = 920.0
        # water density (at rest)
        PP.rho_w = 1000.0 
        # melting temperature
        TP.T_m = 273
    
    elif salts == 1:
        # ice density	
        PP.rho_i = 1130.0	
        # water density (at rest)
        PP.rho_w = 1180.0
        # melting temperature (eutectic)
        TP.T_m = 268

    # Calculate thermal diffusivity (m^2/s)
    TP.kappa = TP.k_ice / (PP.rho_i * TP.cp_s)

    return TP, PP


#---------------------------------------------------------------------
# Reservoir freezing related functions
#---------------------------------------------------------------------

# Stefan's problem solving
def lambdaFunction(l, TP, RES):
    return (l*np.exp(l*l)*(1+sp.special.erf(l))-(TP.cp_s*(273-RES.T)/(TP.Lh*np.sqrt(np.pi))))

# Function that calculates the lambda constant required to solve the Stefan
# problem and returns it:
def lambdaStefan(TP, RES):
    # Stefan problem solving
    Lambda = sp.optimize.fsolve(lambdaFunction, 0.5, args=(TP, RES))
    return Lambda[0]

# Function that calculates the thermal gradient around the reservoir thank to
# the Stefan problem, and returns the thickness of the frozen layer in the
# reservoir S :
def frozenLayer(t, TP, RES):
    # Stefan problem solving
    Lambda = lambdaStefan(TP, RES)
    # thickness of the frozen layer at time t
    S_temp = 2*Lambda*np.sqrt(TP.kappa*t)
    return S_temp

# Inverse funcion of "frozen_layer": this one calculates the time required 
# to freeze a given layer of thickness e
def freezingTime(e, TP, RES):
    # Stefan problem solving
    Lambda = lambdaStefan(TP, RES)
    # time to freeze the thickness e :
    t_temp = (e/(2*Lambda*np.sqrt(TP.kappa)))**2
    return t_temp


#---------------------------------------------------------------------------
# Reservoir freezing
#---------------------------------------------------------------------------

def cryoFreeze(BP,PP,TP,RES):
    
    # lithostatique pressure around the reservoir
    RES.P0 = BP.g*PP.rho_i*RES.h
    
    # temperature around the reservoir
    RES.T = TP.T_min+(((TP.T_max-TP.T_min)/(PP.h_max-PP.h_min))*RES.h)
    
    # ice tensile strenght in the crust
    RES.sigma_c = PP.sigma_top+((PP.sigma_bot-PP.sigma_top)/(PP.h_max-PP.h_min))*RES.h
    
    # overpressure in the reservoir required to open it 
    RES.deltaP_c = 2*(RES.sigma_c + RES.P0)
    print('critical overpressure at 1st order: ', RES.deltaP_c/1e6, ' MPa')
    
    # frozen fraction
    RES.n = (np.exp(PP.beta*RES.deltaP_c)-1)/((PP.rho_w/PP.rho_i)*np.exp(PP.beta*RES.deltaP_c)-1)
    
    # thickness of the critical frozen layer in the reservoir :
    RES.Sc = RES.R*(1-(1-RES.n)**(1/3))
    
    # time required to freeze nc :
    RES.t_c = freezingTime(RES.Sc, TP, RES)
    
    return RES
    
#---------------------------------------------------------------------------
# Deformation model from Dragoni and Magnanensi
#---------------------------------------------------------------------------

def cryoRheol(BP,PP,TP,RES):    
    
    # initial reservoir volume :
    RES.V_i = 4/3*np.pi*RES.R**3
    
    # Compressibilities :
    RES.K1 = PP.E/2*(1+PP.nu)
    RES.K2 = RES.K1
    
    # Radii :
    RES.R1 = RES.R
    RES.R2 = 1.3*RES.R
    
    # Temperature around reservoir :
    #T2 = 100+((150/10000)*h)
    RES.T2 = 200
    
    # Viscosity (only for the viscoelastic zone) :
    RES.eta = 10**14*sym.exp(25.2*(273/RES.T2-1))
    #eta = 10**14*sym.exp(25.2*(273/200-1))
    
    # Constants used for the following calculations
    D = 3*RES.K1*RES.R1**3*(PP.mu1-PP.mu2) - PP.mu1*RES.R2**3*(3*RES.K1+4*PP.mu2)
    D0 = RES.eta*D
    RES.tau = RES.eta * (PP.mu1*(RES.R2/RES.R1)**3*(3*RES.K1+4*PP.mu2)-3*RES.K1*(PP.mu1-PP.mu2))/(3*RES.K1*PP.mu1*PP.mu2)
    
    t, t0, t1, t2, t3, p0 = sym.symbols('t t0 t1 t2 t3 p0')
    
    M = 1 - sym.exp(-t/RES.tau)
    N = H(t-t1)*(1-sym.exp(-(t-t1)/RES.tau))
    O = H(t-t3)*(1-sym.exp(-(t-t3)/RES.tau))
    P = H(t-t2)*(1-sym.exp(-(t-t2)/RES.tau))
    
    f1 = ((t-((t-t1)*H(t-t1)))/t1)+(((t-t3)*H(t-t3)-(t-t2)*H(t-t2))/(t3-t2))
    f2 = RES.tau * ((M-N)/t1+(O-P)/(t3-t2))
    
    # coeff A1 in zones 1 and 2:
    A1_1 = -p0*RES.tau*PP.mu1*RES.R1**3*r*((3*RES.K1+4*PP.mu2)*RES.R2**3/(4*r**3)-PP.mu2)
    A1_2 = -3*p0*RES.tau*RES.K1*PP.mu1*RES.R1**3*RES.R2**3/(4*r**2)
    # coeff A2 in zones 1 and 2:
    A2_1 = -p0*RES.R1**3*r*((3*RES.K1+4*PP.mu2)*(RES.eta-RES.tau*PP.mu1)*RES.R2**3/(4*r**3)+RES.eta*(PP.mu1-PP.mu2)+PP.mu1*PP.mu2*RES.tau)
    A2_2 = -p0*RES.R1**3*RES.R2**3/(4*r**2)*(RES.eta*(3*RES.K1+4*PP.mu2)-3*RES.K1*PP.mu1*RES.tau)
    # coeff B1 in zones 1 and 2:
    B1_1 = 3*p0*RES.tau*RES.K1*PP.mu1*PP.mu2*RES.R1**3
    B1_2 = 3*p0*RES.tau*RES.K1*PP.mu1*PP.mu2*RES.R1**3*RES.R2**3/r**3
    # coeff B2 in zones 1 and 2:
    B2_1 = p0*RES.R1**3*(RES.eta*PP.mu1*(3*RES.K1+4*PP.mu2)*RES.R2**3/r**3 - 3*RES.K1*(RES.eta*(PP.mu1-PP.mu2)+RES.tau*PP.mu1*PP.mu2))
    B2_2 = p0*PP.mu2*RES.R1**3*RES.R2**3/r**3*(RES.eta*(3*RES.K1+4*PP.mu1)-3*RES.K1*PP.mu1*RES.tau)
    # coeff C1 in zones 1 and 2:
    C1_1 = 3*p0*RES.tau*RES.K1*PP.mu1*PP.mu2*RES.R1**3
    C1_2 = -3/2*p0*RES.tau*RES.K1*PP.mu1*PP.mu2*RES.R1**3*RES.R2**3/r**3
    # coeff C2 in zones 1 and 2:
    C2_1 = p0*RES.R1**3*(-RES.eta*PP.mu1*(3*RES.K1+4*PP.mu2)*RES.R2**3/(2*r**3)-3*RES.K1*(RES.eta*(PP.mu1-PP.mu2)+RES.tau*PP.mu1*PP.mu2))
    C2_2 = -1/2*p0*PP.mu2*RES.R1**3*RES.R2**3/r**3*(RES.eta*(3*RES.K1+ 4*PP.mu1)-3*RES.K1*PP.mu1*RES.tau)
    
    RES.u1 = A1_1/D0*f1 + A2_1/D0*f2
    RES.u2 = A1_2/D0*f1 + A2_2/D0*f2
    # print('u(r,t) zone 1 : ', u1)
    # print('u(r,t) zone 2 : ', u2)
    # print('u(r,t) : ok')
    
    RES.sigma1_rr = B1_1/D0*f1 + B2_1/D0*f2
    RES.sigma2_rr = B1_2/0*f1 + B2_2/D0*f2
    # print('sigma_rr(t) zone 1 : ', sigma1_rr)
    # print('sigma_rr(t) zone 2 : ', sigma2_rr)
    # print('simga_rr(r,t) : ok')
    
    RES.sigma1_tt = C1_1/D0*f1 + C2_1/D0*f2
    RES.sigma2_tt = C1_2/D0*f1 + C2_2/D0*f2
    # print('sigma_tt(t) zone 1 : ', sigma1_tt)
    # print('sigma_tt(t) zone 2 : ', sigma2_tt)
    # print('sigma_tt(r,t) : ok')
    
    RES.sigma1_pp = RES.sigma1_tt
    RES.sigma2_pp = RES.sigma2_tt
    # print('sigma_pp(t) zone 1 : ', sigma1_pp)
    # print('sigma_pp(t) zone 2 : ', sigma2_pp)
    # print('sigma_pp(r,t) : ok')

    return RES

# To calculate the reservoir deformation after at a given time:
    
def time2deformation(time, PP, RES):

    # to avoid an error due to identical times:
    t0_temp = 0
    t1_temp = 1*time
    t2_temp = 1.1*time
    t3_temp = 1.2*time

    p0_temp = RES.deltaP_c
    
    t, t0, t1, t2, t3, p0 = sym.symbols('t t0 t1 t2 t3 p0')
    
    u_temp = RES.u1.subs(r,RES.R1).subs(t,t1_temp).subs(t0, t0_temp).subs(t1,t1_temp).subs(t2, t2_temp).subs(t3, t3_temp).subs(p0, p0_temp)
    u_temp = replaceHeaviside(u_temp)
    #print('deformation : ', u_temp*100, ' cm')
    
    R_new = RES.R1 + u_temp
    V_new = 4/3*np.pi*R_new**3
    dP = -1/PP.beta * np.log(float(V_new/RES.V_i))
    #print('pressure drop : ', dP/1e6, ' MPa')
    
    return dP

# To calculate the reservoir deformation after applying a given inner pressure:

def pressure2deformation(deltaP, PP, TP, RES):

    # frozen fraction correspondinf to dP:
    n_temp = (np.exp(PP.beta*deltaP)-1)/((PP.rho_w/PP.rho_i)*np.exp(PP.beta*deltaP)-1)
    # thickness of the critical frozen layer in the reservoir:
    Sc_temp = RES.R*(1-(1-n_temp)**(1/3))
    # time required to freeze nc:
    tc_temp =freezingTime(Sc_temp, TP, RES)
    #print('new freezing time : t_c_new = ', tc_temp/3600/24, ' days')

    t0_temp = 0
    t1_temp = 1*tc_temp
    t2_temp = 1.1*tc_temp
    t3_temp = 1.2*tc_temp
    p0_temp = deltaP
    
    t, t0, t1, t2, t3, p0 = sym.symbols('t t0 t1 t2 t3 p0')
    
    u_temp = RES.u1.subs(r,RES.R1).subs(t,t1_temp).subs(t0, t0_temp).subs(t1,t1_temp).subs(t2, t2_temp).subs(t3, t3_temp).subs(p0, p0_temp)
    u_temp = replaceHeaviside(u_temp)
    #print('deformation : ', u_temp*100, ' cm')
    
    # reservoir radius after deformation and associated pressure drop dP
    R_new = RES.R + u_temp
    V_new = 4/3*np.pi*R_new**3
    dP = -1/PP.beta * np.log(float(V_new/RES.V_i))
    #print('pressure drop : ', dP/1e6, ' MPa')
    
    return (dP, tc_temp)


def iterate(PP, TP, RES):
    # convergence criterion:
    epsilon = RES.deltaP_c/100000
    
    i=0
    RES.i_val.append(0)
    RES.t_val.append(RES.t_c)
    RES.p_val.append(RES.deltaP_c)

    dP = time2deformation(RES.t_c, PP, RES)
    
    
    while i<2:
        i = i+1
        dP_temp = dP
        [dP, t_c] = pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)
        RES.i_val.append(i)
        RES.t_val.append(t_c)
        RES.p_val.append(RES.deltaP_c - dP_temp)
    
    if RES.t_val[2]-RES.t_val[1] > RES.t_val[1]-RES.t_val[0] :
        
        RES.converging = 0
        print('---------- Diverged -----------')
        
    else :
        while dP_temp - dP > epsilon :
            i = i+1
            dP_temp = dP
            [dP, t_c] = pressure2deformation(RES.deltaP_c - dP, PP, TP, RES)
            RES.i_val.append(i)
            RES.t_val.append(t_c)
            RES.p_val.append(RES.deltaP_c - dP_temp)
    
        print('++++++++++ Converged ++++++++++')
        print('number of iterations: ', i-1)

    

        
    return RES

























