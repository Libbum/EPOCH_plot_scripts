'''
Attempt at making Evangelos' 1D density movie function work as I'd like.
Only requires directory and frame. Should then run rust's parallel to get stride information.

Example run case:
	parallel -j 8 python plot1dDens.py -d ../dptest3 -f ::: {1..15..2}
Here, we're using a relative path, which is safely handled, and a bash range of 1 to 15 with a stride of 2.
Forcing the hyperthreads to work seems to give us optimal coverage on this one.

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
output_file = work_dir.joinpath('dens_' + '{:04}'.format(args.frame) + save_ext)

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

# Read all frame data
fdata = sdf.read(str(frame_file), dict=True)

# Calculate required plot variables
xGrid = fdata['Grid/Grid_mid'].data[0] / xScale
tStart = np.abs(fdata['Grid/Grid_mid'].data[0][0]) / c
t = (fdata['Header']['time'] - tStart) / tScale
dx = (fdata['Grid/Grid_mid'].data[0][1] - fdata['Grid/Grid_mid'].data[0][0])

if args.hd:
    monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
    plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
               dpi=monitor_dpi)  # Gives us 1080p output
    mpl.rcParams.update({'font.size': 22})
else:
    plt.figure()

# plt.semilogy(xGrid,fdata['Electric Field/Ex'].data/Ec,
# label='$E_x/E_c$', color='#7ECD0C')
plt.xlabel(r'$x/$' + xUnits)
plt.title(r'$t/$' + tUnits + '$\, =\,$' + '%5.1f' % t)

# Find vector potential
ay = np.zeros(fdata['Magnetic Field/Bz'].data.shape[0])
az = np.zeros(fdata['Magnetic Field/By'].data.shape[0])
a = np.zeros(fdata['Magnetic Field/By'].data.shape[0])

for i in range(ay.shape[0] - 1, 0, -1):
    ay[i - 1] = ay[i] - fdata['Magnetic Field/Bz'].data[i - 1] * dx / aScale
    az[i - 1] = az[i] + fdata['Magnetic Field/By'].data[i - 1] * dx / aScale

a = np.sqrt(ay**2 + az**2)
aE = np.sqrt(fdata['Electric Field/Ey'].data**2 +
             fdata['Electric Field/Ez'].data**2) / Ec

plt.semilogy(xGrid, a, label='$|a|$')

plt.xlim(-100, 10)  # NOTE: These parameters are arbitrary.
plt.ylim(1e0, 1e2)

# Get a string containing all key labels in fdata
fdata_keys = " ".join(fdata.keys())

res = re.findall(r'Derived/Number_Density/[a-z,A-Z,0-9_]+', fdata_keys)
res.sort()
for st in res:
    strLabel = st.replace('Derived/Number_Density/', '')
    plt.semilogy(xGrid, fdata[st].data / nScale,
                 label=r'$n_' + strLabel[0] + '/n_c$')

plt.legend(loc='upper right')
# plt.show()
if args.tikz:
    tikz_save(str(output_file))
else:
    if args.hd:
        plt.savefig(str(output_file), dpi=monitor_dpi)
    else:
        plt.savefig(str(output_file))
plt.close()
