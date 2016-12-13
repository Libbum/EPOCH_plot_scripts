'''
Plots Potential (Φ) phase map.

Requires directory and will use all frames available to it.

Example run case:
	python plot1dPhi.py -d ../dptest3
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

# User required massagers
xlimits = [-5, 225] #Time
ylimits = [-100, 250] #Space TODO: Should have ±200 μm, but we have ±250?

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
args = parser.parse_args()

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

# expand user is needed here in case we get thrown a '~'
work_dir = Path.expanduser(Path(args.tdir))

input_file = work_dir.joinpath('input.deck')  # Gets current input deck
frames = np.sort([x for x in work_dir.glob('*.sdf')])

save_ext = '.png'
if args.tikz:
    save_ext = '.tex'
output_file = work_dir.joinpath('phi' + save_ext)

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

#Initialise required variables
phi = []
xGrid = []
tGrid = []

# Read all frame data
for frame in frames:
    fdata = sdf.read(str(frame), dict=True)

    if len(xGrid) == 0:
        xGrid = fdata['Grid/Grid_mid'].data[0] / xScale #xGrid is common.

    tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0]) / c
    tGrid.append((fdata['Header']['time'] - tStart) / tScale)

    dx = (fdata['Grid/Grid_mid'].data[0][1]-fdata['Grid/Grid_mid'].data[0][0])
    phix = np.zeros(fdata['Electric Field/Ex'].data.shape[0]) # Create matrix

    for i in range(1,phix.shape[0]):
        phix[i]=phix[i-1]-dx*fdata['Electric Field/Ex'].data[i-1]/Ec/xScale

    phi.append(phix)

phi = np.array(phi)
tGrid = np.array(tGrid)

lxidx = find_nearest(tGrid, xlimits[0])
hxidx = find_nearest(tGrid, xlimits[1])
tGrid = tGrid[lxidx:hxidx]

lyidx = find_nearest(xGrid, ylimits[0])
hyidx = find_nearest(xGrid, ylimits[1])
xGrid = xGrid[lyidx:hyidx]

phi = np.transpose(phi[lxidx:hxidx,lyidx:hyidx])

if args.hd:
    monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
    plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
               dpi=monitor_dpi)  # Gives us 1080p output
    mpl.rcParams.update({'font.size': 22})
else:
    plt.figure()

plt.ylabel(r'$x/$' + xUnits)
plt.xlabel(r'$t/$' + tUnits)

plt.pcolormesh(tGrid, xGrid, phi, cmap='viridis')
#plt.pcolormesh(xGrid, pGrid, xpx, vmin=1e14, vmax= 5e19, cmap='plasma', norm=LogNorm())
cbar = plt.colorbar()
cbar.set_label(r'$\Phi$')

#plt.tight_layout()
plt.xlim(xlimits)
plt.ylim(ylimits)

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
