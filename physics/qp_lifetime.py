# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 15:35:19 2016
@author: PeteBarry

Module to calculate quasiparticle lifetime.

"""


from scipy.constants import eV, k; kB = k/eV
import numpy as np


def deltaT(T,Tc):
    """ Calculation of temperature dependent gap energy. """
    del0 = 1.76 * Tc * kB
    return del0 * np.sqrt( np.cos((np.pi/2)*(T/Tc)**2) ) 

def calc_tauqp(T, Tc, tau0 = 430.e-9):
    """ Calculate the quasiparticle lifetime according to Kaplan (1976). """
    delT = deltaT(T, Tc)
    return 1e6 * tau0 / np.sqrt(np.pi) * (kB * Tc / 2 / delT)**2.5 * np.sqrt(Tc / T) * np.exp(delT / kB / T)
    
    