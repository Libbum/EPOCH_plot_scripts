'''
Attempt at making Evangelos' 1D compareSpectra function work as I'd like.
This isn't a multifle comparison, and is sloley attempting to get the energy out.
Only requires directory and frame. Should then run rust's parallel to get stride information.

Example run case:
	parallel -j 8 python plot1DEnergy.py -d ../dptest3 -f ::: {1..15..2}
Here, we're using a relative path, which is safely handled, and a bash range of 1 to 15 with a stride of 2.
Forcing the hyperthreads to not be invoked seems to give us optimal coverage on this one.

ISSUES:
If you get seg faults, it's due to the sdf reader. I've found that if you compile the python bindings with
debug flags on, you don't get the faults anymore. Crazy right? Yeah. I don't know either...

'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sdf
import sys
import re
from matplotlib2tikz import save as tikz_save
from pathlib import *

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

    res = re.search(r'a0\s+=\s+([0-9.]+)', deck)
    a0 = np.float(res.group(1))

# Set up scales and units
xScale = 1e-6
xUnits = '$(\mu m)$'
pScale = 2.73092e-22  # me*c
pUnits = '$(m_e c)$'
aScale = pScale / qe
wL = c * 2 * np.pi / lambdaL
tScale = 2. * np.pi / wL
tUnits = '$\\tau_L$'
EkScale = 1.e6 / (6.24e18)  # Convert to M (1e6), ev (C = 6.24e18*qe)
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
t = (fdata['Header']['time'] - tStart) / tScale

# Distribution functions
# Get a string containing all key labels in fdata
fdata_keys = " ".join(fdata.keys())

for species in re.findall(r'dist_fn/en_cut/([a-z,A-Z,0-9_]+)', fdata_keys):
    print('Processing ' + species)
    d = np.squeeze(fdata['dist_fn/en_cut/' + species].data)

    speciesGrid = 'Grid/en_cut/' + species
    output_file = Path.cwd().joinpath(
        'nrg_' + species + '_' + '{:04}'.format(args.frame) + save_ext)

    if args.hd:
       monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
       plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
                  dpi=monitor_dpi)  # Gives us 1080p output
       mpl.rcParams.update({'font.size': 22})
    else:
       plt.figure()

    En = fdata[speciesGrid].data[0]/EkScale

    dE = En[1]-En[0]

    if linear:
        idxE = (np.abs(En-1.)).argmin() # truncate anything smaller than 1 MeV
        plt.plot(En[idxE:], d[idxE:]*qe/qScale/dE, label=r'$n_0=%.2g$'%n0+', $\\tau_R=$??')
        plt.xlim(xmin=En[idxE])
    else:
        plt.semilogy(En,d*qe/qScale/dE, label=r'$n_0=%.2g$'%n0+', $\\tau_R=$??')

    plt.xlabel('E '+EkUnits)
    plt.ylabel('f(E)  '+ '(arbitrary units)')

    plt.tight_layout()

    plt.title('t/'+tUnits+' = ' +'%.3g'%t +',  $a_0=$'+'%.2g'%a0+',  '+species)
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
