import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
from matplotlib import ticker

#%% MCM j-values

l = 5.804e-6
m = 1.092
n = 0.377

x = np.arange(np.pi,3*np.pi,0.001)
plt.plot(l * np.cos(x)**m * np.exp(-n/np.cos(x)))

#%%
    
# biacet
l = 3.326e-4 	
m = 0.148
n = 0.215

x = np.arange(-np.pi,np.pi,0.001)
cosz = np.cos(x)

plt.plot(x,l * cosz**m * np.exp(-n/cosz))
plt.xlim(-np.pi,np.pi)

#%% MCM cross-sections and quantum-yield

cjx_weff = np.array([187.,      191.,      193.,      196.,      202.,      208.,
          211.,      214.,      261.,      267.,      277.,      295.,
          303.,      310.,      316.,      333.,      383.,      599.])

# BIACET

folder = '/home/willemvc/Desktop/MCM rates/BIACET/'
file = 'biacet_horowitz01_cs_plum83_qy_298.txt'
save_f = '/home/willemvc/Desktop/projects/EMEP/CloudJ/figures/cross_secs/'

qyield = 0.158
txt = np.loadtxt(folder + file,skiprows=16)
fac = 10**-20 # cm2 / molec

wlens = txt[:,0]
xsects = txt[:,1] * qyield * fac

ind_min = np.argmin(np.abs(cjx_weff - wlens[0]))
ind_max = np.argmin(np.abs(cjx_weff - wlens[-1]))

vals = np.interp(cjx_weff,wlens,xsects)
vals[:ind_min] = 0.
if wlens[-1] < cjx_weff[ind_max]:
    vals[ind_max:] = 0.
  
vals = np.round(vals,23)

plt.title('BIACET (1 bar, 298K) q-yield = 0.158')
plt.xlabel('Wavelength (nm)')
plt.plot(wlens,xsects / fac,label='Cross-section data from MCM database')
plt.plot(cjx_weff,vals/fac,'.',label='CJX-bin effective Wavelength')
plt.ylabel(r'cross-sect $\times$ q-yield (1e-20 cm$^2$ molec$^{-1}$)')
plt.legend()

plt.savefig(save_f + 'BIACET.png', bbox_inches='tight',dpi=300)

#%%
