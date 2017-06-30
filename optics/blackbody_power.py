# Program to calculate the power radiated onto each stage of Aloysius from the 300 K window.
# Filters on current stages (11/02/12)
import numpy as np
from scipy.constants import h, c, pi, k as kb
from filters import generate_dummy_filter, Filter, FilterStack

def Planck(x, Tbb, units = 'Frequency'):
    if units == 'Wavenumber':
        return (2e8*h*c**2*x**3) / (np.exp(100*h*c*x/(kb*Tbb)) - 1)  # in units of W/m^2/cm^-1/sr
    elif units == 'Wavelength':
       return (2*h*c**2/(x**5))*(1/(np.exp((h*c)/(kb*Tbb*x)) - 1))  # in units of W/m^3/sr
    elif units == 'Frequency':
       	return (2*h*x**3/(c**2))*(1/(np.exp((h*x)/(kb*Tbb)) - 1))  # in units of W/m^3/sr
    else:
        raise AttributeError ('Incorrect units entered for Planck function. Check spelling and try again.')

# Todo
# Make array of these for each resonator in array
 
class OpticalPower(object):
    """Class to store and calculate the incident power from a blackbody source
    through a filter stack. Will calculate solid angle and take into account position of array
    with respect to the aperture. Optionally will account for antenna beam pattern (not implemented yet).
    
    Requires: 
        - FilterStack instance
    Optional:
        - Tbb (default = [3, 100, 1000, 1.]): list/array containing [Tmin, Tmax, Tstep, Epsilon]
        - Rsource (1): Radius of the source/limiting aperture
        - Distance (1): Distance to the source/limiting aperture
        - lpixel/rpixel (1): Depending on pixel geometry, lpxiel for square, rpixel is radius of circular pixel
        - (wavelength, m) (None, 1.): If given, the throughput is calculated as m * lambda**2, instead of geometrically 
        
    """
    def __init__(self, filterstack = None, Tbb = [3., 100, 1000, 1.], **kwargs):
        self._init_filterstack(filterstack)
        
        Rsource = kwargs['Rsource']      if 'Rsource'  in kwargs else 1
        Distance = kwargs['Distance']    if 'Distance' in kwargs else 1
        if 'lpixel' in kwargs.keys():   Adet = kwargs['lpixel']**2
        elif 'rpixel' in kwargs.keys(): Adet = pi*kwargs['rpixel']**2
        else:  Adet = 1.
        
        (wavelen, m) = kwargs['wavelen'] if 'wavelen' in kwargs and type(kwargs['wavelen']) in [list, tuple, np.array] else (None, 1)

        self.Aomega = self.calc_solid_angle(Rsource, Distance, Adet, (wavelen, m))
        self.Ptot = self.Aomega*self.calc_blackbody_power(*Tbb) # in units of Watts
        
        if 'temperatures' in kwargs.keys(): self.Pbb = self.interp_power(kwargs['temperatures'])

    def _init_filterstack(self, filterstack, **kwargs):
        if type(filterstack) == FilterStack:
            print 'Filter Stack Found'
            self.FilterStack = filterstack
        elif type(filterstack) in [str, unicode, np.unicode_]:
            print 'Filename found. Attempting to read into FilterStack object'
            if 'xltabs' in kwargs.keys(): xltabs = kwargs['xltabs']
            else: xltabs = 'all'
            self.FilterStack = FilterStack(filename = filterstack, xltabs=xltabs)    
        else:
            print 'No filter stack found.'
            self.FilterStack = None   # make null filter stack (3 THz LPE for example)        
            
    def calc_solid_angle(self, Rsource, Distance, Adet, mlambda2 = None):
        """Need geometry of cryostat, and pixel
        Need to account for pixel position (off axis solid angle - numerical)
        """
        if all(mlambda2) and len(mlambda2) == 2 : 
            wavelen, m = mlambda2            
            return m * wavelen**2
            
        else:
            #lpixel = pixelgeometry
            #Rsource, D = cryogeometry
            #Adet = lpixel**2
            theta = np.arctan(Rsource/Distance) # half angle of the cone
            print 'Angle = %.2f rad'%theta
            omega = 2*pi*(1 - np.cos(theta))
            return Adet*omega
        
    def calc_blackbody_power(self, Tmin, Tmax, step = 1000., epsilon = 1.):
        """Calculates the power from a blackbody source for range of temperatures 
        with emiisivity = epsilon, across the frequency span defined by the filter stack.
        
        The units returned to 
            - Ptot [W m-2 sr-1]
            - _Pbb [W m-2 sr-1 Hz-1]
        """
        nugrid = self.FilterStack.nugrid[:,np.newaxis] # new axis to perform the integration faster using numpy axis arguement in trapz
        transtot = self.FilterStack.transtot[:,np.newaxis]       
        Tbb = np.linspace(Tmin, Tmax, step)[np.newaxis, :]
        
        Pbb = epsilon*Planck(nugrid, Tbb, units='Wavenumber')*transtot
        self.Tbb = Tbb.flatten(); self._Pbb = Pbb
        return np.trapz(Pbb, nugrid, axis=0)
        
    def interp_power(self, Tbb):
        return np.interp(Tbb, self.Tbb, self.Ptot)


# --- SuperSpec thesis run
f1 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/B330_300cm-1_LPE.txt'
f2 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/Filter_1420_W945_200um_BP.txt'
f3 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/W1543_MARS2_Safari_MCS_Band2_LPE.txt'

# --------- 2mm dummy (05/11/15 Design Calculation) ---------
if False:
    frequency = 150e9
    lamda = c/frequency
    fs = FilterStack([Filter(data = generate_dummy_filter(50., 3e3, frequency/1.e9, tmax=0.5), name = '2mm Dummy BP')])    
    lpixel = lamda
    Rsource = 52.e-3/2 # source radius
    D = 75.e-3
# ---------------------------------------------------------------

# --------- SK-EOR4-2MHz-A (26/08/2015 Measurement) ---------
if False:
    fname = u'/Users/PeteBarry/Documents/Projects/SpaceKIDs/FilterData/SpaceKIDs/350GHz_BB_Configuration/SK_350GHz_BB_stack.xls'
    fs = FilterStack(filename=fname, xltabs='all')    
    lpixel = 1.3e-3
    D = 49.84e-3
    Rsource = 7.5e-3 # source radius
# ---------------------------------------------------------------

# --------- FTS Sanity Check (19/11/2015) ---------
if False:
    fs = FilterStack([Filter(data = generate_dummy_filter(50., 5e3, 3.e3, tmax=1, kind='LPF'), name = '3 THz Edge')])
    lpixel = 50e-3
    D = 101.6e-3
    Rsource = 1.3e-3/2 # source radius
# ---------------------------------------------------------------
    
        

#------------------------------------------------ 
#----- calculation of offaxis solid angle -------
#------------------------------------------------ 
if False:
    ## IN PROGRESS
    # need to calculate if point x is with disk radius, or outside
    x0 = linspace(0, 10.e-3, 10)# location of detector
    x0 = 0.# location of detector
    #rdisk = 23.8e-3/2
    r0 = linspace(0, 1, 11)
    rm = 1
    
    L = 1 # normal distance from x to disk
    R1 = sqrt(L**2 + (r0 - rm)**2) # shortest distance from x to edge of disk
    Rmax = sqrt(L**2 + (r0 + rm)**2) # longest distance from x to edge of disk
    
    m = 1 - R1**2/Rmax**2
    mp = 1 - m
    a = 4*r0*rm/(r0+rm)**2
    
    if any(a) == 0: where
    phi = arcsin(((1 - m/a)/(1-m))**0.5) if a != 0. else 0.
    # if true

    omega = 2*pi - 2*L/Rmax * ellipk(m) - pi*Heumann(phi, m)

def Heumann(phi, m):
    mp = 1 - m
    Ek  = ellipe(m);# print Ek
    Ekp = ellipeinc(phi, mp);# print Ekp
    Kk  = ellipk(m); #print Kk
    Fkp = ellipkinc(phi, mp);# print Fkp

    return 2/pi*( Ek*Fkp + Kk*Ekp - Kk*Fkp )
                   




#--------------------------------------------------------------------------------------------------#





