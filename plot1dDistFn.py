import matplotlib.pyplot as plt
import matplotlib as mlib
import numpy as np
import sdf
from matplotlib.colors import LogNorm
from scipy.integrate import simps
#import epochUtil as epu

plt.ion()

my_jet = plt.cm.get_cmap('jet')
my_jet.set_under('w')

#cspList=np.loadtxt("/home/siminos/articles/srsDyn/libraries/matplotlib/ColorSpiral1.dat")
#cSpiral=mlib.colors.LinearSegmentedColormap.from_list('csp_color',cspList, N=256)


lambdaL = 1e-6
xScale = lambdaL/(2*np.pi) #1e-6
xUnits = '$[\lambda_L/2\,\pi]$'
pScale = 2.73092e-22 # me*c
pUnits = '$[m_e c]$'
tScale = 1e-15
tUnits = 'fs'
EkScale = 1.e6/(6.24e18)
EkUnits = '(MeV)'
fUnits = '(a.u.)'
qe = 1.60218e-19
aScale = pScale/qe
c = 2.998e8
wL = c*2*np.pi/lambdaL
Ec = pScale*wL/qe

xStp=1
yStp=1

xmin=-5
xmax=15

pxmin = -65
pxmax = 105

n0 = 4.5
a0 = 15.

tStart=0



#dir ='/data/finite/siminos/epoch/SIT_a015.0_n007.0_test_env'
dir ='./'
#dir = '/scratch/siminos/SIT/SIT_a015.0_n007.0_LP_test_env'

i= 200
fname='%04d.sdf'%i
#a=SDF(dir+'/'+fname)
#b=a.read()
b = sdf.read(dir+fname, dict=True)
xGrid = b['Grid/Grid_mid'].data[0]/xScale
#xGrid = b['Grid/Particles/electrons'].data[0]/xScale

#xGrid=b['Grid/Grid/X']/xScale

t = b['Header']['time']/tScale

#plt.figure()

#plt.plot(xGrid,b['Derived/Number_Density/electrons'])

#plt.show()
#plt.plot(xGrid,b['Derived/Number_Density/ions'])

plt.figure()

plt.plot(xGrid,b['Electric Field/Ex'].data)

plt.show()
plt.figure()

plt.plot(xGrid,np.sqrt(b['Electric Field/Ey'].data**2+b['Electric Field/Ez'].data**2)/Ec)


plt.show()
# Distribution functions

c=b['dist_fn/x_px/electrons'].data

d=np.reshape(c,[c.shape[0],c.shape[1]])

mx=np.max(d)

mlib.rcParams.update({'font.size': 28})


fig=plt.figure(figsize=(16, 16), dpi=150)

extxpx=[b['Grid/x_px/electrons'].data[0][0]/xScale,b['Grid/x_px/electrons'].data[0][-1]/xScale,b['Grid/x_px/electrons'].data[1][0]/pScale,b['Grid/x_px/electrons'].data[1][-1]/pScale]

plt.imshow(np.transpose(d),cmap=my_jet, vmin=1e15, vmax=mx, extent= extxpx, origin = 'lower', aspect=0.1, norm=LogNorm())

plt.plot([0, 0], [pxmin, pxmax],  color = '0.75')

#plt.title('t = ' +'%.3f'%(t-tStart)  +" "+ tUnits)

plt.ylim(pxmin, pxmax)
plt.xlim(xmin,xmax)

clb=plt.colorbar(shrink=0.6)

plt.xlabel(r'x '+xUnits)
plt.ylabel(r'$p_x$'+pUnits)
clb.set_label(r'$f$'+fUnits)

plt.savefig(dir+'/PS.eps',bbox_inches='tight')
