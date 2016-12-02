'''
Builds a maximum energy graph for a number of runs. This requires a good deal of processing, so may take some time.
For the moment this is written specifically to deal with the split pulse experiments, I'll generalise it more if needed in the future.

To run: call this on a directory structure like the following:
split
  |---- 0030
  |---- 0060
  | ...
  |---- 1070
where each folder from the `split` folder is indiciative of the \Delta t time between split pulses.

Requires a start directory. Here this will be the path to `split`, will output one figure in the `split` directory.

ASSUMPTIONS:
File path is similar to the one suggested above, where folder names are \Delta t in femtosectons, and all runs are complete.

LIMITATIONS:
Split pulse specific, only cares about protons. Requires dist_fn/en/protons; however will soft fail for frames that don't carry it.

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
from pathlib import *

parser = argparse.ArgumentParser(
    description="Plot the maximum proton energy for each run within chosen folder.")
parser.add_argument('-v', '--verbose', help="Noisy output",
                    action='store_true', dest='loud', required=False, default=False)
parser.add_argument('-hd', help="1080p output", action='store_true',
                    dest='hd', required=False, default=False)
parser.add_argument('-d', help="Base directory pointing to target directories.",
                    action="store", dest="tdir", required=False, default='.')
args = parser.parse_args()

# expand user is needed here in case we get thrown a '~'
work_dir = Path.expanduser(Path(args.tdir))
#numdirs = sum(1 for x in work_dir.glob('*') if x.is_dir())
dirs = (x for x in work_dir.glob('*') if x.is_dir())

c = 2.998e8
qe = 1.60218e-19
me = 9.10938291e-31
eps0 = 8.854187817e-12
tUnits = '(fs)'
EkScale = 1.e6 / (6.24e18)  # Convert to M (1e6), ev (C = 6.24e18*qe)
EkUnits = '(MeV)'

# Plot value initialisations
dt = []
Emax = []

for cdir in dirs:
    maxEn = 0.0

    dt_dir = work_dir.joinpath(cdir)

    frames = (x for x in dt_dir.glob('*.sdf'))
    for frame in frames:
        fdata = sdf.read(str(frame), dict=True)
        if 'dist_fn/en_full/protons' in fdata:
            # The sdf only writes grid values with associated data. This is our
            # max energry for this frame.
            En = fdata['Grid/en_full/protons'].data[0] / EkScale
            if En[-1] > maxEn:
                maxEn = En[-1]

    # NOTE: it works for now, but we assume sanitary folders here.
    dt.append(int(str(cdir.stem)))
    Emax.append(maxEn)

if args.hd:
    monitor_dpi = 96  # NOTE: This is for Jnana, may be different for Brahman.
    plt.figure(figsize=(1920 / monitor_dpi, 1080 / monitor_dpi),
               dpi=monitor_dpi)  # Gives us 1080p output
    mpl.rcParams.update({'font.size': 22})
else:
    plt.figure()

plt.scatter(dt, Emax)
plt.xlabel('Time delay ' + tUnits)
plt.ylabel('Maximum proton energy  ' + EkUnits)

plt.tight_layout()

if args.hd:
    plt.savefig(str(work_dir.joinpath('maxnrg.png')), dpi=monitor_dpi)
else:
    plt.savefig(str(work_dir.joinpath('maxnrg.png')))
plt.close()
