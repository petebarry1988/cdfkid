# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 15:35:19 2016
@author: PeteBarry

Module to calculate quasiparticle lifetime.

"""


from scipy.constants import pi, eV, k; kB = k/eV
import numpy as np


def deltaT(T,Tc):
    """ Calculation of temperature dependent gap energy. """
    del0 = 1.76 * Tc * kB
    return del0 * np.sqrt( np.cos((np.pi/2)*(T/Tc)**2) ) 

def calc_tauqp_kaplan(T, Tc, tau0 = 430.e-9):
    """ Calculate the quasiparticle lifetime according to Kaplan (1976). Value return is in seconds. """
    delT = deltaT(T, Tc)
    return tau0 / np.sqrt(pi) * (kB * Tc / 2 / delT)**2.5 * np.sqrt(Tc / T) * np.exp(delT / kB / T)
    
def calc_tauqp(T, Tc, tau0 = 430.e-9, taumax = 1.6e-3):
    """ Calculate the quasiparticle lifetime including the empircal addition  """
    tauqpT = calc_tauqp_kaplan(T, Tc, tau0)    
    return 1. / (1/tauqpT + 1/taumax)
    
def nqp(T, Tc, N0 = 1.72e10):
    """ Calculate expected quasiparticle density. This is calculated from the approximate expression, assuming
    a thermal distribution in the limit kB*T<< delta.  
    Units are in eV - 
        - N0 [eV^-1 um^-3]
        - delta [eV] """
    del0 = deltaT(T, Tc)
    return 2 * N0 * np.sqrt(2 * pi * kB * T* del0) * np.exp(- del0/ (kB*T) ) # eV^-1 * um-3 * sqrt(eV* eV) = um-3

