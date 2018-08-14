import numpy as np
import oommf_sk_number as osn
import re
import os
import glob
import subprocess
import shutil
this_dir = os.path.dirname(os.path.abspath(__file__))


def generate_omfs():

    """
    Simulate a skyrmion reducing the mesh discretisation
    using the OOMMF script: isolated_sk_DMI_Cnv.mif

    It is necessary to have OOMMF installed with the latest DMI
    modules. Recommended to install JOOMMF's OOMMF conda package
    """

    OMF_DIR = os.path.join(this_dir, 'omfs/')
    if os.path.exists(OMF_DIR):
        shutil.rmtree(OMF_DIR)
    os.makedirs(OMF_DIR)

    for n in range(20, 101, 20):
        SIM_NAME = 'omfs/isolated_sk_Cnv_n_{:03d}'.format(n)
        SCRIPT = os.path.join(this_dir, 'isolated_sk_DMI_Cnv.mif')

        job = ('oommf boxsi -threads 2 '
               '-parameters "NX {0} '
               'BASENAME {1}" '
               '"{2}"'.format(n, SIM_NAME, SCRIPT)
               )
        print(job)

        subprocess.call(job, shell=True)


# -----------------------------------------------------------------------------


# Function to sort OMF files according to number of mesh sites n from filename
def get_n(f):
    return int(re.search('(?<=n_)\d+(?=-Oxs)', f).group(0))


def test_sk_number_vs_mesh_discretisation():
    """
    """

    generate_omfs()

    # The files in the right order
    _files = glob.glob(os.path.join(this_dir, 'omfs/*.omf'))
    _files = sorted(_files, key=get_n)

    mesh_length = 60  # nm

    for FILE in _files:
        n = get_n(FILE)
        oommf_data = osn.SkNumberOOMMF(FILE)
        print('dx = {:.2f} nm --> Q = {}'.format(mesh_length / n,
                                                 oommf_data.compute_sk_number())
              )
        oommf_data.plot_system(savefig='test_n_{:03d}.png'.format(n))


if __name__ == '__main__':
    test_sk_number_vs_mesh_discretisation()
