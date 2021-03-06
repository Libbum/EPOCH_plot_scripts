'''
Plots Potential (Φ) for a given frame.

Requires directory and frame. Should then run rust's parallel to get stride information.

Example run case:
	parallel python plot1dPhi.py -d ../dptest3 -f ::: {1..15..2}
Here, we're using a relative path, which is safely handled, and a bash range of 1 to 15 with a stride of 2.

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
from matplotlib.colors import LogNorm

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
output_file = Path.cwd().joinpath('phi_{:04}'.format(args.frame) + save_ext)

labmbdaL = 0.0
with input_file.open(encoding='utf-8') as inputFile:
    deck = inputFile.read().replace('\n', '')

    # NOTE: This isn't that robust really. But I think it'll be fine so long
    # as we stick to the micron convention.
    res = re.search(r'lambda\s*=\s*([0-9.]+)\s*\*\s*micron', deck)
    lambdaL = np.float(res.group(1)) * 1e-6

# Set up constants. NOTE: Constants directly from Evangelos' file
xScale = 1e-6
xUnits = '$(\mu m)$'
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

xGrid = fdata['Grid/Grid_mid'].data[0] / xScale #xGrid is common.
tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0]) / c
t = (fdata['Header']['time'] - tStart) / tScale

dx = (fdata['Grid/Grid_mid'].data[0][1]-fdata['Grid/Grid_mid'].data[0][0])
phi = np.zeros(fdata['Electric Field/Ex'].data.shape[0]) # Create matrix

for i in range(1,phi.shape[0]):
    phi[i]=phi[i-1]-dx*fdata['Electric Field/Ex'].data[i-1]/Ec/xScale

phi = np.array(phi)

if args.hd:
    monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
    plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
               dpi=monitor_dpi)  # Gives us 1080p output
    mpl.rcParams.update({'font.size': 22})
else:
    plt.figure()

plt.ylabel(r'$x $' + xUnits)
plt.xlabel(r'$t/$' + tUnits)
plt.title('t/'+tUnits+' = ' +'%.3g'%t)
plt.plot(xGrid, phi)

#plt.tight_layout()
plt.xlim(-100, 100)
plt.ylim(-1.5, 2.5)

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
