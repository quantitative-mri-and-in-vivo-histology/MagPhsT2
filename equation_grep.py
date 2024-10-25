#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 09:37:05 2023

@author: wangd
"""
import numpy as np
import matplotlib.pyplot as plt 

import pyepg

#from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import least_squares 

def greps_signal_pyepg(T1, T2, alpha, delta_phi, TR, ABS_ERR_TOL = 1e-8 ):
    EPG = pyepg.PyEPG()
    EPG.SetParameters(1.0, T1, T2)
    EPG.Equilibrium()
    N = EPG.StepsToSS(alpha, delta_phi, TR, ABS_ERR_TOL)
    #print('N=%d, phase=%7.5f, M=(%7.5f, %7.5f)' % (N,EPG.GetPhase(),EPG.GetReA(),EPG.GetImA() ))
    signal = (EPG.GetReA()+1.0j*EPG.GetImA())*np.exp(-1.0j*np.deg2rad(EPG.GetPhase()))
    return ( np.abs(signal), np.angle(signal) -np.pi/2 )


def greps_diff_signal_err_least_squares(x,fa,tr,dphi,signal):
    t1 = x[0]
    t2 = x[1]
    am = x[2] #amplitude modulation
    er = np.zeros((2*len(dphi),)) # error function of double length (real/imag values)
    n  = len(dphi)   #assume, sorted positive-only dphi values:  [dphi(0) ..., dph_(n-1), dph_n]  
    for i in range(n):
        a,p               = greps_signal_pyepg(t1,t2,fa,dphi[i],tr,1e-6)
        epgsig            = am * a * np.exp(-1.0j * p)
        er[2*i]           = np.real(epgsig-signal[i])
        er[2*i+1]         = np.imag(epgsig-signal[i])

    return er

