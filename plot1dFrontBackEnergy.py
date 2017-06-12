'''
Same as plot1DEnergy, however we have an assumption that we're plotting separate front and backside
proton layers here.

ASSUMPTIONS:
requires two proton species, named protons_front and protons_back included in the en_cut dist_fn

Example run case:
parallel -j 8 python plot1DEnergy.py -d ../dptest3 -f ::: {1..15..2}

Here, we're using a relative path, which is safely handled, and a bash range of 1 to 15 with a stride of 2.
Forcing the hyperthreads to not be invoked seems to give us optimal coverage on this one.

ISSUES:
If you get seg faults, it's due to the sdf reader. I've found that if you compile the python bindings with
debug flags on, you don't get the faults anymore. Crazy right? Yeah. I don't know either...

'''

import os
import matplotlib as mpl
#Don't invoke X if we're running remotely.
if os.environ.get('DISPLAY','') == '':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sdf
import sys
import re
from matplotlib2tikz import save as tikz_save
from statsmodels.nonparametric.smoothers_lowess import lowess
from pathlib import *

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


parser = argparse.ArgumentParser(
    description="Print a single laser/plasma interaction frame showing particle energy.")
parser.add_argument('-v', '--verbose', help="Noisy output",
                    action='store_true', dest='loud', required=False, default=False)
parser.add_argument('-hd', help="1080p output", action='store_true',
                    dest='hd', required=False, default=False)
parser.add_argument('-t', '--tikz', help="Save output as tikz", action='store_true',
                    dest='tikz', required=False, default=False)
parser.add_argument('-d', help="Directory pointing to target files.",
                    action="store", dest="tdir", required=False, default='.')
parser.add_argument('-f', help="Frame to print to image.",
                    action="store", dest="frame", type=int, required=True)
args = parser.parse_args()

# expand user is needed here in case we get thrown a '~'
work_dir = Path.expanduser(Path(args.tdir))

input_file = work_dir.joinpath('input.deck')  # Gets current input deck
# TODO: This :04 system only works with file outputs of the order I've
# been looking at it recently. But this isn't always the case.
# gets file we want to print data from
frame_file = work_dir.joinpath('{:04}'.format(args.frame) + '.sdf')

save_ext = '.png'
if args.tikz:
    save_ext = '.tex'

# Get constants ready
labmbdaL = 0.0
n0 = 0.0
a0 = 0.0
c = 2.998e8
qe = 1.60218e-19
me = 9.10938291e-31
eps0 = 8.854187817e-12

# NOTE: We don't pull ramp -> tauR here, because I have NFI what that represents, nor do I use it.
# Get values we need from the deck
with input_file.open(encoding='utf-8') as inputFile:
    deck = inputFile.read().replace('\n', '')

    # NOTE: This isn't that robust really. But I think it'll be fine so long
    # as we stick to the micron convention.
    res = re.search(r'lambda\s*=\s*([0-9.]+)\s*\*\s*micron', deck)
    lambdaL = np.float(res.group(1)) * 1e-6

    omega = 2.0 * np.pi * c / lambdaL
    nc = omega**2 * me * eps0 / qe**2
    # NOTE: This only works if we have the density in terms of n_c
    res = re.search(r'den_init\s*=\s*([0-9.]+)\s*\*\s*den_crit', deck)
    n0 = np.float(res.group(1)) * nc
    #n0 = 115.2 * nc

    res = re.search(r'a0\s+=\s+([0-9.]+)', deck)
    a0 = np.float(res.group(1))

# Set up scales and units
xScale = lambdaL  # 1e-6
xUnits = '$\lambda_L$'
pScale = 2.73092e-22  # me*c
pUnits = '$(m_e c)$'
aScale = pScale / qe
wL = c * 2 * np.pi / lambdaL
tScale = 2. * np.pi / wL
tUnits = '$\\tau_L$'
EkScale = 1.e6 / (6.24e18)  # Convert to Mev: M (1e6), ev (C = 6.24e18*qe)
EkUnits = '(MeV)'
Ec = pScale * wL / qe
nC = eps0 * me * wL * wL / (qe * qe)
nScale = nC
qScale = 1e-12
qUnits = '(pC)'

xWindow = 9.  # plot region length
xForward = 6.  # plot region to show preceding xAf
linear = 0  # Linear scale plot on/off

# Read all frame data
fdata = sdf.read(str(frame_file), dict=True)

# Calculate required plot variables
xGrid = fdata['Grid/Grid_mid'].data[0] / xScale
tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0]) / c
t0 = 8e-14 + tStart #80 fs is 2*40fs which should be laser peak
t = (fdata['Header']['time']-t0) / tScale
tcomp = fdata['Header']['time']

output_file = Path.cwd().joinpath(
    'nrg_protons_' + '{:04}'.format(args.frame) + save_ext)

if args.hd:
   monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
   plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
              dpi=monitor_dpi)  # Gives us 1080p output
   mpl.rcParams.update({'font.size': 22})
else:
   plt.figure()

back = np.array(np.squeeze(fdata['dist_fn/en_cut/protons_back'].data))
front = np.array(np.squeeze(fdata['dist_fn/en_cut/protons_front'].data))

En_back = np.array(fdata['Grid/en_cut/protons_back'].data[0]/EkScale)
En_front = np.array(fdata['Grid/en_cut/protons_front'].data[0]/EkScale)

#idxs = find_closest(En_front, En_back) #get an index of where the front array fits into the back array.
E_common = np.linspace(min(min(En_back),min(En_front)), max(max(En_back),max(En_front)), En_front.shape[0])
back_common = np.interp(E_common, En_back, back, left=0, right=0)
front_common = np.interp(E_common, En_front, front, left=0, right=0)

dE = E_common[1]-E_common[0]

total = (front_common + back_common)*qe/qScale/dE

filteredf = lowess(front_common*qe/qScale/dE, E_common, is_sorted=True, frac=0.0075, it=0)
filteredb = lowess(back_common*qe/qScale/dE, E_common, is_sorted=True, frac=0.0075, it=0)
filteredt = lowess(total, E_common, is_sorted=True, frac=0.0075, it=0)

if linear:
    idxE = (np.abs(En-1.)).argmin() # truncate anything smaller than 1 MeV
    plt.plot(En_back[idxE:], total[idxE:]*qe/qScale/dE, label='protons')
else:
    plt.semilogy(E_common,back_common*qe/qScale/dE, label='back', c=(0.894118, 0.101961, 0.109804))
    plt.semilogy(E_common,front_common*qe/qScale/dE, label='front', c=(0.215686, 0.494118, 0.721569))
    plt.semilogy(E_common,total, c=(0.67451, 0.298039, 0.415686), label='total')
    #plt.semilogy(filteredb[:,0], filteredb[:,1], c=(0.894118, 0.101961, 0.109804), label='back')
    #plt.semilogy(filteredf[:,0], filteredf[:,1], c=(0.215686, 0.494118, 0.721569), label='front')
    #plt.semilogy(filteredt[:,0], filteredt[:,1], c=(0.67451, 0.298039, 0.415686), label='total')

plt.xlabel('E '+EkUnits)
plt.ylabel('f(E)  '+ '(arbitrary units)')

#plt.tight_layout()
plt.ylim(5e8,1e15)
plt.xlim(0,8)

plt.title('t0+'+'%.3g'%t + tUnits +' ('+'%.3g'%tcomp +' s),  $a_0=$'+'%.2g'%a0)
plt.legend(loc='best')

# plt.show()
if args.tikz:
    tikz_save(str(output_file))
else:
    if args.hd:
       plt.savefig(str(output_file), dpi=monitor_dpi)
    else:
       plt.savefig(str(output_file))
plt.close()
