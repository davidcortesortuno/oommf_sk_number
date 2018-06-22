"""
Script to compute timestep and Q from a directory of OMF files

For example, the charge Q for every time step of a simulation
where snapshots are saved in a folder called my_oommf_sim containing:

    MySim_00-Oxs_TimeDriver-Magnetization-0000000-0010000.omf
    MySim_00-Oxs_TimeDriver-Magnetization-0000005-0060005.omf
    MySim_00-Oxs_TimeDriver-Magnetization-0000010-0110010.omf

can be processed using this cript as:

    python compute_snap-Q_from_folder.py my_oommf_sim --timestep 1e-12

in case timesteps were specified in OOMMF as 1 ps. This script will then
generate a text file with two columns: 

    time        sk_number

    0 * 1e-12      ..
    5 * 1e-12      ..
   10 * 1e-12      ..


You can also pass the --plot_snaps or --plot_charge options to
plot snapshots for every OMF file



Created by: D. Cortes-Ortuno on Tue 27 Mar 2018 12:11:22 BST
            d.cortes@soton.ac.uk
            university of Southampton

Copyright D. Cortes-Ortuno
"""

import os
import argparse
import oommf_sk_number as oskn
import re
import numpy as np

# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser('Save time and Q from a directory of OMF files')

parser.add_argument('path')
parser.add_argument('--timestep', default=1., type=float)
parser.add_argument('--plot_charge', action='store_true')
parser.add_argument('--plot_snaps', action='store_true')

args = parser.parse_args()

# -----------------------------------------------------------------------------

if os.path.isfile(args.path):
    TypeError('Specify a directory')

file_list = sorted(os.listdir(args.path),
                   key=lambda f: int(re.search('(?<=Magnetization-)[0-9]+', f).group(0)) )

data = np.zeros((len(file_list), 2))

if args.plot_snaps:
    if not os.path.exists('pngs'):
        os.makedirs('pngs')
if args.plot_charge:
    if not os.path.exists('pngs_charge'):
        os.makedirs('pngs_charge')

sim_name = re.search('.*(?=-Oxs)', file_list[0]).group(0)
for i, _file in enumerate(file_list):

    snapshot = int(re.search('(?<=Magnetization-)[0-9]+', _file).group(0))

    print('Computing Q for snap: {:06d}'.format(int(snapshot)))

    oommf_data = oskn.SkNumberOOMMF(os.path.join(args.path, _file))

    if args.plot_snaps:
        oommf_data.plot_system(savefig='pngs/{}_{:06d}'.format(sim_name, snapshot) + '.png')
    if args.plot_charge:
        oommf_data.plot_charge_density(savefig='pngs_charge/{}_charge_{:06d}'.format(sim_name, snapshot) + '.png')

    Q = oommf_data.compute_sk_number()
    data[i][0] = snapshot * args.timestep
    data[i][1] = Q

print('-' * 80)
print('Saving data')
np.savetxt(sim_name + '_snap-Q.txt', data)
print('Finished!')
