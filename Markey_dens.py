'''
Plotting Markey Fig 2 values. At the moment this is a brute force test adventure. Not good in general.
We need t0 + 100fs. For the dptest run from Hebbe that's 0968.sdf (time = 1.7e-12s)
And for the split pulse we need t0 + 50fs (note, t0 here is for the second pulse). This is at 3.6e-12, which is 1253.sdf in the dpulse dir.
so we want to call
	python Markey_dens.py -d1 ~/epoch/hebbe/doublepulse/dptest -f1 0968 -d2 ~/epoch/hebbe/doublepulse/dpulse -f2 1253
These files are now stored on the external drive under hebbe/splitpulse/Markey_tests/
'''

import matplotlib as mpl
import argparse
from pathlib import *
import matplotlib.pyplot as plt
import numpy as np
import sdf
import sys
import re

parser = argparse.ArgumentParser(description="Print a single laser/plasma interaction frame showing plasma density.")
parser.add_argument('-v','--verbose', help="Noisy output", action='store_true', dest='loud', required=False, default=False)
parser.add_argument('-d1', help="Directory pointing to target files of top run.", action="store", dest="tdir", required=True)
parser.add_argument('-d2', help="Directory pointing to target files of bottom run.", action="store", dest="tdir2", required=True)
parser.add_argument('-f1', help="Frame of top run.", action="store", dest="frame", type=int, required=True)
parser.add_argument('-f2', help="Frame of bottom run.", action="store", dest="frame2", type=int, required=True)
args = parser.parse_args()

work_dir = Path.expanduser(Path(args.tdir)) #expand user is needed here in case we get thrown a '~'
work2_dir = Path.expanduser(Path(args.tdir2))

input_file = work_dir.joinpath('input.deck') #Gets current input deck
frame_file = work_dir.joinpath('{:04}'.format(args.frame)+'.sdf') #gets file we want to print data from
frame2_file = work2_dir.joinpath('{:04}'.format(args.frame2)+'.sdf')
output_file = Path.cwd().joinpath('markey.png')

labmbdaL = 0.0
with input_file.open(encoding='utf-8') as inputFile:
	deck = inputFile.read().replace('\n','')

	res = re.search(r'lambda\s+=\s+([0-9.]+)\s+\*\s+micron',deck) #NOTE: This isn't that robust really. But I think it'll be fine so long as we stick to the micron convention.
	lambdaL = np.float(res.group(1))*1e-6

#Set up constants. NOTE: Constants directly from Evangelos' file
xScale = lambdaL #1e-6
xUnits = '$\lambda_L$'
pScale = 2.73092e-22 # me*c
pUnits = '$(m_e c)$'
qe = 1.60218e-19
me = 9.10938291e-31
eps0 = 8.854187817e-12
aScale = pScale/qe
c = 2.998e8
wL = c*2*np.pi/lambdaL
tScale =  2.*np.pi/wL
tUnits = '$\\tau_L$'
Ec = pScale*wL/qe
nC =  eps0*me*wL*wL/(qe*qe)
nScale = nC

#Read all frame data
fdata = sdf.read(str(frame_file), dict=True)

#Calculate required plot variables
xGrid = fdata['Grid/Grid_mid'].data[0]/xScale
tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0])/c
t = (fdata['Header']['time']-tStart)#/tScale
dx = (fdata['Grid/Grid_mid'].data[0][1]-fdata['Grid/Grid_mid'].data[0][0])

f2data = sdf.read(str(frame2_file), dict=True)
x2Grid = f2data['Grid/Grid_mid'].data[0]/xScale
t2Start = np.abs(f2data['Grid/Grid_mid'].data[0][0])/c
t2 = (f2data['Header']['time']-t2Start)#/tScale
dx2 = (f2data['Grid/Grid_mid'].data[0][1]-f2data['Grid/Grid_mid'].data[0][0])

f, axarr = plt.subplots(2,3)
axarr[0][0].plot(fdata['Grid/Particles/back_protons'].data[0]/xScale, fdata['Particles/Px/back_protons'].data, label=r'protons')
axarr[0][0].plot(fdata['Grid/Particles/ions'].data[0]/xScale, fdata['Particles/Px/ions'].data, label=r'ions')
axarr[0][0].set_ylabel(r'Electric Field ($V/m$)')
axarr[0][0].set_xlabel(r'$x$ ($\mu m$)')
axarr[0][0].set_ylim(-1e-20,4e-20)
axarr[0][0].set_xlim(5,20)

axarr[0][1].semilogy(xGrid, fdata['Derived/Number_Density/back_protons'].data, label=r'protons')
axarr[0][1].semilogy(xGrid, fdata['Derived/Number_Density/ions'].data, label=r'ions')
axarr[0][1].set_title(r'$t \, =\, t_0+100$ fs')
axarr[0][1].set_ylabel(r'particle density ($m^3$)')
axarr[0][1].set_xlabel(r'$x$ ($\mu m$)')
axarr[0][1].set_ylim(1e24,1e30)
axarr[0][1].set_xlim(5,20)
axarr[0][1].legend(loc='upper right')

axarr[0][2].plot(xGrid, fdata['Electric Field/Ex'].data)
axarr[0][2].set_ylabel(r'Momentum ($Kg m$)')
axarr[0][2].set_xlabel(r'$x$ ($\mu m$)')
axarr[0][2].set_ylim(0,3e12)
axarr[0][2].set_xlim(5,20)

axarr[1][0].plot(f2data['Grid/Particles/back_protons'].data[0]/xScale, f2data['Particles/Px/back_protons'].data, label=r'bs protons')
axarr[1][0].plot(f2data['Grid/Particles/ions'].data[0]/xScale, f2data['Particles/Px/ions'].data, label=r'ions')
axarr[1][0].plot(f2data['Grid/Particles/front_protons'].data[0]/xScale, f2data['Particles/Px/front_protons'].data, label=r'fs protons')
axarr[1][0].set_ylabel(r'Electric Field ($V/m$)')
axarr[1][0].set_xlabel(r'$x$ ($\mu m$)')
axarr[1][0].set_ylim(-2e-20,10e-20)
axarr[1][0].set_xlim(5,35)

axarr[1][1].semilogy(x2Grid, f2data['Derived/Number_Density/back_protons'].data, label=r'bs protons')
axarr[1][1].semilogy(x2Grid, f2data['Derived/Number_Density/ions'].data, label=r'ions')
axarr[1][1].semilogy(x2Grid, f2data['Derived/Number_Density/front_protons'].data, label=r'fs protons')
axarr[1][1].set_title(r'$t \, =\, t_0+50$ fs')
axarr[1][1].set_ylabel(r'particle density ($m^3$)')
axarr[1][1].set_xlabel(r'$x$ ($\mu m$)')
axarr[1][1].set_ylim(1e24,1e30)
axarr[1][1].set_xlim(5,35)
axarr[1][1].legend(loc='upper right')

axarr[1][2].plot(x2Grid, f2data['Electric Field/Ex'].data)
axarr[1][2].set_ylabel(r'Momentum ($Kg m$)')
axarr[1][2].set_xlabel(r'$x$ ($\mu m$)')
axarr[1][2].set_ylim(0,3e12)
axarr[1][2].set_xlim(5,35)


# plt.figure()
# plt.ylabel(r'particle density ($m^3$)')
# plt.xlabel(r'$x$ ($\mu m$)')
# plt.title(r'$t \, =\, t_0$')
#
#
# plt.semilogy(xGrid, fdata['Derived/Number_Density/ions'].data, label=r'ions')
# plt.semilogy(xGrid, fdata['Derived/Number_Density/back_protons'].data, label=r'protons')
#
# plt.ylim(1e23,1e30)
# plt.xlim(5,20) #NOTE: These parameters are arbitrary.
#
# plt.legend(loc='upper right')
plt.show()

#plt.savefig(str(output_file))
#plt.close()
