'''
Plots Electric Field and its FFT for a given frame.

Requires directory and frame. Should then run rust's parallel to get stride information.

Example run case:
	parallel python tt_Ex_fft.py -d ../dptest3 -f ::: {1..15..2}
Here, we're using a relative path, which is safely handled, and a bash range of 1 to 15 with a stride of 2.

ISSUES:
If you get seg faults, it's due to the sdf reader. I've found that if you compile the python bindings with
debug flags on, you don't get the faults anymore. Crazy right? Yeah. I don't know either...
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import scipy.fftpack
import numpy as np
import argparse
import sdf
import sys
import re
from matplotlib2tikz import save as tikz_save
from pathlib import *
from matplotlib.colors import LogNorm

xlimits = [-5, 15] #To speed up fft calculation

parser = argparse.ArgumentParser(
    description="Print a single laser/plasma interaction frame showing plasma density.")
parser.add_argument('-hd', help="1080p output", action='store_true',
                    dest='hd', required=False, default=False)
parser.add_argument('-t', '--tikz', help="Save output as tikz", action='store_true',
                    dest='tikz', required=False, default=False)
parser.add_argument('-v', '--verbose', help="Noisy output",
                    action='store_true', dest='loud', required=False, default=False)
parser.add_argument('-d', help="Directory pointing to target files.",
                    action="store", dest="tdir", required=False, default='.')
parser.add_argument('-f', help="Frame to print to image.",
                    action="store", dest="frame", type=int, required=True)
args = parser.parse_args()

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

# expand user is needed here in case we get thrown a '~'
work_dir = Path.expanduser(Path(args.tdir))
frame_file = work_dir.joinpath('{:04}'.format(args.frame) + '.sdf')
input_file = work_dir.joinpath('input.deck')  # Gets current input deck

save_ext = '.png'
if args.tikz:
    save_ext = '.tex'
output_file = Path.cwd().joinpath('Ex_fft_{:04}'.format(args.frame) + save_ext)

labmbdaL = 0.0
with input_file.open(encoding='utf-8') as inputFile:
    deck = inputFile.read().replace('\n', '')

    # NOTE: This isn't that robust really. But I think it'll be fine so long
    # as we stick to the micron convention.
    res = re.search(r'lambda\s*=\s*([0-9.]+)\s*\*\s*micron', deck)
    lambdaL = np.float(res.group(1)) * 1e-6

# Set up constants. NOTE: Constants directly from Evangelos' file
xScale = lambdaL  # 1e-6
xUnits = '$\lambda_L$'
pScale = 2.73092e-22  # me*c
pUnits = '$(m_e c)$'
qe = 1.60218e-19
me = 9.10938291e-31
eps0 = 8.854187817e-12
aScale = pScale / qe
c = 2.998e8
wL = c * 2 * np.pi / lambdaL
tScale = 2. * np.pi / wL
tUnits = '$\\tau_L$'
Ec = pScale * wL / qe
nC = eps0 * me * wL * wL / (qe * qe)
nScale = nC

fdata = sdf.read(str(frame_file), dict=True)

xGrid = fdata['Grid/Grid_mid'].data[0] / 1e-6 #to μm
tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0]) / c
t = (fdata['Header']['time'] - tStart) / tScale

dx = (fdata['Grid/Grid_mid'].data[0][1]-fdata['Grid/Grid_mid'].data[0][0])
Ex = np.array(fdata['Electric Field/Ex'].data)

lxidx = find_nearest(xGrid, xlimits[0])
hxidx = find_nearest(xGrid, xlimits[1]+1)
xGrid = xGrid[lxidx:hxidx]

Ex = Ex[lxidx:hxidx]

if args.hd:
    monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
    plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
               dpi=monitor_dpi)  # Gives us 1080p output
    mpl.rcParams.update({'font.size': 22})
else:
    plt.figure()

#ax1 = plt.gca()
#ax2 = ax1.twinx()

plt.xlabel(r'$x$ ($\mathrm{\mu}$m)')
plt.title('t/'+tUnits+' = ' +'%.3g'%t)

FFT = abs(scipy.fft(Ex))
freqs = scipy.fftpack.fftfreq(Ex.size, dx/50)

ax1 = plt.subplot(211)
ax1.plot(xGrid, Ex, '#e41a1c')
ax1.set_ylabel(r'$E_x$ (V/m)')
#Bounds of box
ax1.plot((3, 3), (-5e12, 5e12), 'k--')
ax1.plot((6, 6), (-5e12, 5e12), 'k--')
x = np.linspace(0,3,10)
y = (np.exp(2*x)/400)*5e12
expl, = ax1.plot(x, y, '0.75') #Return as a tuple not a list
expl.set_linestyle('--')


ax2 = plt.subplot(212)
#ax2.plot(freqs,scipy.log10(FFT),'x', color='#377eb8')
ax2.semilogy(freqs[int(freqs.size/2):-1],FFT[int(freqs.size/2):-1], color='#377eb8')
#plt.show()

#ax1.plot(xGrid, Ex, '#e41a1c')

#ax2.plot(xGrid, densPF, '#984ea3')
#ax2.plot(xGrid, densE, '#984ea3')
#ax2.plot(xGrid, densPB, '#377eb8')
#ax2.plot(xGrid, densI, '#4daf4a')
#ff7f00

#plt.tight_layout()
ax1.set_xlim(xlimits)
ax1.set_ylim(-5e12, 5e12)
ax2.set_ylim(1e9, 1e16)
ax2.set_xlim(-5e8, 5e8)

if args.tikz:
    if args.hd:
        tikz_save(str(output_file), dpi=monitor_dpi)
    else:
        tikz_save(str(output_file))
else:
    if args.hd:
        plt.savefig(str(output_file), dpi=monitor_dpi)
    else:
        plt.savefig(str(output_file))
plt.close()
