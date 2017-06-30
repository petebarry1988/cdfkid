# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 14:00:22 2016

@author: PeteBarry
"""
# --- Module imports ---
import numpy as np
#from scipy.constants import h, c, k, eV

def generate_dummy_filter(fmin, fmax, fc, tmax = 0.9, npoints=1000, kind = 'BPF', bandwidth = 0.1):
    """ Function to generate dummy filter profile. 
    
    Required parameters:
        - fmin, fmax: minimum and maximum frequency, in GHz,
        - fc: cut off frequency in GHz,
        - tmax: in-band transmission.
    
    Optional parameters:
        - npoints = 1000: number of points in the filter profile
        - kind = 'BPF': Type of filter response - 'HPF', 'LPF' or 'BPF',
        - bandwidth = 0.1: Bandwidth of passband. Used only for use with kind = 'BPF'.
        
    Returns:
        np.array([Wavenumber, Transmission]) with shape [npoints, 2].
        
    -----------------------------------------------------------------------------------
        """

    nu = np.linspace(fmin, fmax, npoints)/30.
    nuc = fc/30.
    if kind == 'HPF':   trans = tmax*(np.sign(nu-nuc) + 1)/2
    elif kind == 'LPF': trans = tmax*(1 - np.sign(nu-nuc))/2
    elif kind == 'BPF': 
        bw = nuc*bandwidth/2
        trans = np.zeros_like(nu)
        trans[(nu > nuc-bw) & (nu < nuc+bw)] = tmax
    else: print 'Type of filter not recognised.'; return    
    
    return np.array([nu, trans])
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

class Filter(object):
    """Class to store a single a filter object from either data or a file. 
        
    Required parameters:
        - fmin, fmax: minimum and maximum frequency, in GHz,
        - fc: cut off frequency in GHz,
        - tmax: in-band transmission.
    
    Optional parameters:
        - npoints = 1000: number of points in the filter profile
        - kind = 'BPF': Type of filter response - 'HPF', 'LPF' or 'BPF',
        - bandwidth = 0.1: Bandwidth of passband. Used only for use with kind = 'BPF'.
        
    Returns:
        np.array([Frequency, Transmission]) with shape [npoints, 2].
        """
        
    def __init__(self, filename = None, xltab = None, data = None, name = None):
        from os.path import splitext
        self.name = name
        self.date = None
        
        if data is not None:
            if len(data) != 2: print 'Data should be 2d array of wavenumber and transmission'
            else: self.__init_from_data__(*data)
            
        elif filename: 
            if splitext(filename)[1] == '.xls':
                print 'Reading xls file...'
                self.__read_from_xls__(filename, xltab) 
            else:
                print 'Formats other than .xls are not currently supported'
                return
                            
    def __init_from_data__(self, nu, trans):
        from collections import namedtuple
        data = namedtuple('FilterData','nu, fr, trans')  
        self.data = data(nu, nu*30, trans)
    
    def __read_from_xls__(self, filename,  xltab):
        import xlrd
        from numpy import array
        workbook = xlrd.open_workbook(filename)
        worksheet = workbook.sheet_by_name(xltab)
        nu, trans = array(worksheet.col_values(0,1)), array(worksheet.col_values(1,1))
        self.name = worksheet.cell_value(0,1).rstrip()
        try:
            date = xlrd.xldate_as_tuple(worksheet.cell_value(0, 2), workbook.datemode)
            date = '%02d/%02d/%d'%(date[2],date[1],date[0])
        except IndexError:
            date = 'Date not given'
        self.__init_from_data__(nu, trans)
    
    def plot_filter(self, fignum=None, xaxis='nu'):
        from matplotlib.pyplot import fignum_exists, subplots, gca, gcf
        if fignum_exists(fignum): fig = gcf(); ax = gca()
        else: fig, ax = subplots(1, 1, num=fignum)
        
        if xaxis=='nu': ax.plot(self.data.nu, self.data.trans, label = self.name)
        elif xaxis=='fr': ax.plot(self.data.fr, self.data.trans, label = self.name)
        
        return fig, ax

####################################################################################
class FilterStack(object):
    """Class to store and process a number of filters, resulting in total transmission"""

    def __init__(self, filters=None, filename=None, xltabs=None):
        from numpy import array
        self.stack=None
        self.items, self._items = [],[]
        self.nugrid = array([])
        
        if filters:
            self.append(filters); self._refresh()
        if filename and xltabs:
            self.__gen_stack_from_xls__(filename, xltabs); self._refresh()

    def __gen_stack_from_xls__(self, filename, xltabs):
        import xlrd
        if xltabs=='all':
            xltabs = xlrd.open_workbook(fname).sheet_names()
        
        if type(xltabs) not in [list]: xltabs=[xltabs]
        
        self.append([Filter(filename=filename, xltab=xltab) for xltab in xltabs if not xltab.startswith('_')])
            
    def append(self, filters):
        if type(filters) in [list, np.array]:
            print 'Processing a list of Filters...'
            for f in filters:
                self.items.append(f.name)
                self._items.append(f)
        
        elif type(filters) in [Filter]:
            print 'Processing a single Filter...'
            if not filters.name: filters.name = 'Unknown'
            self.items.append(filters.name)
            self._items.append(filters)
        
        print 'Current stack:\n', '\n'.join(self.items)
        self._refresh()

    def remove(self,selection):
        # remove filter at the specified index
        if type(selection)==int:
            if selection<len(self.items):
                self.items.pop(selection)
                self._items.pop(selection)
                self._refresh()
            else:
                print "FILTER STACK ERROR - specified index doesn't exist"
                return
    
    def interpolate(self):
        numin = min([f.data.nu.min() for f in self._items])
        numax = max([f.data.nu.max() for f in self._items])
        self.nugrid = np.linspace(numin, numax, 5000)     
        
    def plot_stack(self, xaxis='nu'):
        # todo - plot total on top
        for f in self._items:
            fig, ax = f.plot_filter(fignum=1, xaxis=xaxis)
        if xaxis == 'nu': m = 1.; xax = u'Wavenumber (cm\u207B\u00B9)'
        elif xaxis == 'fr': m = 30.; xax = 'Optical Frequency (GHz)'
        ax.plot(self.nugrid*m, self.transtot, label='Total')
        ax.set_xlabel(xax, fontsize=14); 
        ax.set_ylabel(u'Normalised Transmission',fontsize=14)

        ax.legend()

    def _refresh(self, npoints = 10000):
        numin = min([f.data.nu.min() for f in self._items])
        numax = max([f.data.nu.max() for f in self._items])
        self.nugrid = np.linspace(numin, numax, npoints)
        transgrid = np.array( [np.interp(self.nugrid, f.data.nu, f.data.trans) for f in self._items] )
        self.transtot = np.prod(transgrid, axis=0)



# --- 2mm dummy (05/11/15 Design Calculation) ---------
if False:
    fs = FilterStack([Filter(data = generate_dummy_filter(50., 3e3, frequency/1.e9, tmax=0.5), name = '2mm Dummy BP')])
    
    lpixel = lamda
    Rsource = 52.e-3/2 # source radius
    D = 75.e-3

# --- SK-EOR4-2MHz-A (26/08/2015 Measurement) ---------
if False:
    fs = FilterStack(filename=fname, xltabs='all')    
    lpixel = 1.3e-3
    D = 49.84e-3
    Rsource = 7.5e-3 # source radius

# --- FTS Sanity Check (19/11/2015) ---------
if False:
    fs = FilterStack([Filter(data = generate_dummy_filter(50., 5e3, 3.e3, tmax=1, kind='LPF'), name = '3 THz Edge')])
    lpixel = 50e-3
    D = 101.6e-3
    Rsource = 1.3e-3/2 # source radius
# ---------------------------------------------------------------
if False:
    Adet = lpixel**2
    #Adet = pi*(lpixel/2.)**2
    theta = arctan(Rsource/D) # half angle of the cone
    omega = 2*pi*(1-cos(theta))
    
    Aomega = Adet*omega
    #Aomega = lamda**2
    
    TBB = linspace(3, 100, 1000)
    Pt = array([trapz(Planck(fs.nugrid, t, units='Wavenumber')*Aomega*fs.transtot, fs.nugrid) for t in TBB])
    
    fig1, ax1 = subplots(1,1)
    
    ax1.semilogy(TBB, Pt*1e12, lw=2, label='LEKID')
    ax1.set_xlabel('Blackbody Temperature (K)', fontsize=14)
    ax1.set_ylabel(u'Incident Power (pW)', fontsize=14) # sr\u207B\u00B9
    ax1.set_ylim(1e-2, 50)
    import matplotlib.ticker
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    
    import matplotlib as mpl
    show()
    mpl.pyplot.get_current_fig_manager().window.raise_()
    
    
# ------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------- Power Calculation -------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------
 
#class OpticalPower(object):
#    """Class to store and calculate the incident power from a blackbody source
#    through a filter stack. Will calculate solid angle and take into account position of array
#    with respect to the aperture. Optionally will account for antenna beam pattern (not implemented yet).
#    
#    Requires: 
#        - FilterStack instance
#    Optional:
#        - Tbb (default = [3, 100, 1000, 1.]): list/array containing [Tmin, Tmax, Tstep, Epsilon]
#        - Rsource (1): Radius of the source/limiting aperture
#        - Distance (1): Distance to the source/limiting aperture
#        - lpixel/rpixel (1): Depending on pixel geometry, lpxiel for square, rpixel is radius of circular pixel
#        - (wavelength, m) (None, 1.): If given, the throughput is calculated as m * lambda**2, instead of geometrically 
#        
#    """
#    def __init__(self, filterstack = None, Tbb = [3., 100, 1000, 1.], **kwargs):
#        self._init_filterstack(filterstack)
#        
#        Rsource = kwargs['Rsource']      if 'Rsource'  in kwargs else 1
#        Distance = kwargs['Distance']    if 'Distance' in kwargs else 1
#        if   'lpixel' in kwargs.keys():  Adet = kwargs['lpixel']**2
#        elif 'rpixel' in kwargs.keys():  Adet = np.pi*kwargs['rpixel']**2
#            
#        (wavelen, m) = kwargs['wavelen'] if 'wavelen' in kwargs and type(kwargs['wavelen']) in [list, tuple, array] else (None, 1)
#
#        self.Aomega = self.calc_solid_angle(Rsource,Distance, Adet, (wavelen, m))
#        self.Ptot = self.Aomega*self.calc_blackbody_power(*Tbb) # in units of Watts
#               
#    def _init_filterstack(self, filterstack):
#        if type(filterstack) == FilterStack:
#            self.FilterStack = filterstack
#        elif type(filterstack) in [str, unicode, unicode_]:
#            if 'xltabs' in kwargs.keys(): xltabs = kwargs['xltabs']
#            else: xltabs = 'all'
#            self.FilterStack = FilterStack(filename = filterstack, xltabs=xltabs)    
#        else:
#            print 'No filter stack found.'
#            self.FilterStack = None   # make null filter stack (3 THz LPE for example)        
#            
#    def calc_solid_angle(self, Rsource, Distance, Adet, mlambda2 = None):
#        """Need geometry of cryostat, and pixel
#        Need to account for pixel position (off axis solid angle - numerical)
#        """
#        if all(mlambda2) and len(mlambda2) == 2 : 
#            wavelen, m = mlambda2            
#            return m * wavelen**2
#            
#        else:
#            #lpixel = pixelgeometry
#            #Rsource, D = cryogeometry
#            #Adet = lpixel**2
#            theta = arctan(Rsource/Distance) # half angle of the cone
#            print 'Angle = %.2f rad'%theta
#            omega = 2*pi*(1-cos(theta))
#            return Adet*omega
#        
#    def calc_blackbody_power(self, Tmin, Tmax, step = 1000., epsilon = 1.):
#        """Calculates the power from a blackbody source for range of temperatures 
#        with emiisivity = epsilon, across the frequency span defined by the filter stack.
#        
#        The units returned to 
#            - Ptot [W m-2 sr-1]
#            - _Pbb [W m-2 sr-1 Hz-1]
#        """
#        nugrid = self.FilterStack.nugrid[:,np.newaxis] # new axis to perform the integration faster using numpy axis arguement in trapz
#        transtot = self.FilterStack.transtot[:,np.newaxis]       
#        Tbb = linspace(Tmin, Tmax, step)[newaxis, :]
#        
#        Pbb = epsilon*Planck(nugrid, Tbb, units='Wavenumber')*transtot
#        self.Tbb = Tbb.flatten(); self._Pbb = Pbb
#        return trapz(Pbb, nugrid, axis=0)
#        
#    def interp_power(self, Tbb):
#        return interp(Tbb, self.Tbb, self.Ptot)
#
#
#




# For SuperSpec thesis run
f1 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/B330_300cm-1_LPE.txt'
f2 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/Filter_1420_W945_200um_BP.txt'
f3 = '/Volumes/Macintosh HD/Users/Pete/Documents/PhD/Data/Filters/W1543_MARS2_Safari_MCS_Band2_LPE.txt'

# SpaceKIDs 350 GHz fitler stack
fname = u'/Users/PeteBarry/Documents/Projects/SpaceKIDs/FilterData/SpaceKIDs/350GHz_BB_Configuration/SK_350GHz_BB_stack.xls'
f1tab = 'T1090r13'
f2tab = 'S3194R16'
f3tab = 'T1018R4'
f4tab = 'S3144R22'

flist = [f1,f2,f3]






