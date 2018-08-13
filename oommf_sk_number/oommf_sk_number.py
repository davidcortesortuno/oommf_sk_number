"""
Class to calculate topological charge from a OMF file

Created by: D. Cortes-Ortuno on Tue 27 Mar 2018 12:11:22 BST
            d.cortes@soton.ac.uk
            university of Southampton

Copyright D. Cortes-Ortuno (A license will be created soon)
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import colorsys


# -----------------------------------------------------------------------------


def convert_to_RGB(hls_color):
    return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                        hls_color[1],
                                        hls_color[2]))


def generate_RGBs(field_data):
    """
    field_data      ::  (n, 3) array
    """
    hls = np.ones_like(field_data)
    hls[:, 0] = np.arctan2(field_data[:, 1],
                           field_data[:, 0]
                           )
    hls[:, 0][hls[:, 0] < 0] = hls[:, 0][hls[:, 0] < 0] + 2 * np.pi
    hls[:, 1] = 0.5 * (field_data[:, 2] + 1)
    rgbs = np.apply_along_axis(convert_to_RGB, 1, hls)

    # Some RGB values can get very small magnitudes, like 1e-10:
    # ???
    # rgbs[rgbs < 0] += 2 * np.pi

    return rgbs

# -----------------------------------------------------------------------------


class SkNumberOOMMF(object):
    """
    Compute Sk Number in a slice of a 3D sample from a OMF file

    input_file      :: Path to OMF file
    z_index         :: Index of the slice, starting from the bottom

    TODO:
        - Calculate coordinates
        - Print which z coordinate corresponds to the chosen z_index
        - Work with PBC?
        - More documentation

    """

    def __init__(self, input_file, z_index=0):

        self.input_file = input_file
        self.z_index = z_index

        # Compute nx, ny, dx, dy
        self.get_headers()
        # Extract data from OMF file and create a grid with spins
        self.compute_data()
        # Generate neighbours array
        self.generate_ngbs()

        self.sk_number = np.zeros((self.nx, self.ny))

    def compute_data(self):
        """
        Extract data from a slice in the XY plane
        """
        self.raw_data = np.loadtxt(self.input_file, comments='#')
        self.raw_data_norm = np.sqrt(self.raw_data[:, 0] ** 2 +
                                     self.raw_data[:, 1] ** 2 +
                                     self.raw_data[:, 2] ** 2)
        no_zero = self.raw_data_norm > 0.0
        self.raw_data[no_zero] = self.raw_data[no_zero] / self.raw_data_norm[no_zero][:, np.newaxis]

        # Spin data in a grid, dimensions are: z, y, x, 3
        self.spin_grid = self.raw_data.reshape(-1, self.ny, self.nx, 3)

        # Get the specified slice along the z dimension
        self.spin_grid = self.spin_grid[self.z_index]

        # Data in a N x 3 matrix (might be unnecesary)
        self.data = self.spin_grid.reshape(self.ny * self.nx, 3)
        self.data_norm = np.sqrt(self.data[:, 0] ** 2 +
                                 self.data[:, 1] ** 2 +
                                 self.data[:, 2] ** 2)

    def get_headers(self):
        """
        Get information of the mesh grid from the header/comments
        of the OMF file
        """

        _file = open(self.input_file)

        # Generate a single string with the whole header up to the line where
        # numerical Data starts
        line = _file.readline()
        data = ''
        while not line.startswith('# Begin: Data Text'):
            data += line
            line = _file.readline()

        attrs = {'xstepsize': 'dx',  'ystepsize': 'dy', 'zstepsize': 'dz',
                 'xmin': 'xmin', 'ymin': 'ymin', 'zmin': 'zmin',
                 'xmax': 'xmax', 'ymax': 'ymax', 'zmax': 'zmax',
                 }

        # Regex search the attributes. Stepsizes are specified as dx, dy, dz
        for k in attrs.keys():
            num_val = float(re.search('(?<={}: )[0-9\-\.e]+'.format(k),
                            data).group(0))
            setattr(self, attrs[k], num_val)

        # Compute number of elements in each direction
        self.nx = int(round((self.xmax - self.xmin) / self.dx))
        self.ny = int(round((self.ymax - self.ymin) / self.dy))
        self.nz = int(round((self.zmax - self.zmin) / self.dz))

        _file.close()

    def set_coordinates(self):
        """
        Compute the coordinates of the system using the information
        from the header of the OMF file. Coordinates are stored in
        1D arrays for every component (x,y,z) using the variables

                self.x, self.y, self.z

        Coordinates are scaled in nm units
        The 1D arrays follow the order of the magnetisation/spin data
        """
        xs, ys, zs = (np.arange(float(self.nx)),
                      np.arange(float(self.ny)),
                      np.arange(float(self.nz))
                      )
        xs *= self.dx
        xs += self.xbase
        ys *= self.dy
        ys += self.ybase
        zs *= self.dz
        zs += self.zbase

        xs = np.tile(np.tile(xs, self.ny), self.nz)
        ys = np.tile(np.repeat(ys, self.nx), self.nz)
        zs = np.repeat(zs, self.nx * self.ny)

        # self.coordinates = np.ravel(np.column_stack((xs, ys, zs))) * 1e9
        self.coordinates = np.column_stack((xs, ys, zs)) * 1e9
        self.x, self.y, self.z = (self.coordinates[:, 0],
                                  self.coordinates[:, 1],
                                  self.coordinates[:, 2])

    def load_system_data(self):
        mag_data = np.loadtxt(self.input_file, comments='#')
        mag_data_grid = mag_data.reshape(-1, self.nx, 3)
        return mag_data_grid

    def _index_2D(self, i, j):
        """
        Returns the index for the cell with ordinals i, j
        or -1 if that cell would be out of bounds.

        """
        if i < 0 or j < 0 or j >= self.ny or i >= self.nx:
            return -1

        return j * self.nx + i

    def generate_ngbs(self):
        nx, ny = self.nx, self.ny
        ngbs = np.zeros((nx * ny, 4), dtype=np.int32)

        for j in range(ny):
            for i in range(nx):
                ngbs[i + j * nx] = [self._index_2D(i, j - 1),
                                    self._index_2D(i, j + 1),
                                    self._index_2D(i - 1, j),
                                    self._index_2D(i + 1, j),
                                    ]

        self.neighbours = ngbs

    def get_spin(self, i, j):
        """
        Get spin components from self.spin_grid given the i,j
        mesh grid indexes
        """
        if i < 0 or j < 0 or j >= self.ny or i >= self.nx:
            return np.zeros(3)
        return self.spin_grid[i, j]

    def compute_sk_number(self):
        """

        Compute the skyrmion number S, defined as:
                                _
                     1         /       dm     dm
             S   =  ---  *    /   m .  --  X  --   dx dy
                    4 PI   _ /         dx     dy

        for a two dimensional layer (we use self.spin_grid)

        A finite differences discretisation of the continuum magnetisation
        field, using central differences, and a simple midpoint rule
        for the 2D integration leads to:

        S =   -(  M_i \cdot ( M_{i+1} \times M_{j+1} )
                + M_i \cdot ( M_{i-1} \times M_{j-1} )
                - M_i \cdot ( M_{i-1} \times M_{j+1} )
                - M_i \cdot ( M_{i+1} \times M_{j-1} )
                ) / (16 * PI)
        """

        # 2nd argument are how many zeroes (before,after) we pad in each axis
        # (we keep 3-spin-components, so we don't pad anything at axis=2)
        spin_pad = np.pad(self.spin_grid, ((1, 1), (1, 1), (0, 0)),
                          mode='constant', constant_values=0.0)
        # Same as doing:
        # spin_pad = np.zeros((nx +2, ny + 2, 3))
        # spin_pad[1:-1, 1:-1, :] = spin_grid

        # Here we vectorise the cross products using the padded matrix to
        # obtain neighbours (which are zero spin components) at the boundary of
        # the sample
        self.sk_number = (np.cross(spin_pad[2:, 1:-1, :],   # s(i+1,j)
                                   spin_pad[1:-1, 2:, :],   # s(i,j+1)
                                   axis=2) +
                          np.cross(spin_pad[:-2, 1:-1, :],  # s(i-1,j)
                                   spin_pad[1:-1, :-2, :],  # s(i,j-1)
                                   axis=2) -
                          np.cross(spin_pad[:-2, 1:-1, :],  # s(i-1,j)
                                   spin_pad[1:-1, 2:, :],   # s(i,j+1)
                                   axis=2) -
                          np.cross(spin_pad[2:, 1:-1, :],   # s(i+1,j)
                                   spin_pad[1:-1, :-2, :],  # s(i,j-1)
                                   axis=2)
                          )

        # The dot product of every site with the cross product between
        # their neighbours that was already computed above
        # We save this quantity to the self.sk_number method

        # self.sk_number = -np.sum(self.spin_grid * self.sk_number,
        #                          axis=2) / (16 * np.pi)
        self.sk_number = -np.einsum('ijk,ijk->ij',
                                    self.spin_grid,
                                    self.sk_number) / (16 * np.pi)

        # Total sk number (integral)
        return np.sum(self.sk_number.flatten())

    def plot_system(self, savefig=None, cmap='hls', cbar=True, dpi=150,
                    cbar_offsets=[-0.15, -0.2]):

        f = plt.figure()
        ax = f.add_subplot(111)

        if cmap == 'hls':
            spin_data = generate_RGBs(self.data[:, :3])
            spin_data[self.data_norm < 1e-10] = [1., 1., 1.]
            spin_data = spin_data.reshape(-1, self.nx, 3)

            p = ax.imshow(spin_data, origin='lower', interpolation='None',
                          vmin=0, vmax=2 * np.pi,
                          extent=[self.xmin * 1e9, self.xmax * 1e9,
                                  self.ymin * 1e9, self.ymax * 1e9]
                          )
            if cbar:
                box = ax.get_position()
                axColor = plt.axes([box.x1 + cbar_offsets[0],
                                    box.y1 + cbar_offsets[1],
                                    0.2, 0.2], projection='polar')
                # For the x axis, we need an extra element so rgb has 1 less
                # element along the y direction (might change in the future)
                azimuths = np.arange(0, 361, 1)
                zeniths = np.arange(20, 50, 1)
                dz = 30

                # colours have 1 less element along x
                rgb = np.ones((dz * 360, 3))
                # Set the HLS hue value from 0 to 2 PI from the azimuth values
                # We tile the circle 30 times:
                #   [0 ... 2PI] -> [0...2PI 0 .. 2PI ...]
                rgb[:, 0] = np.tile(azimuths[:-1] * np.pi / 180, dz)
                # For every circle (360 values) we increase the Light value
                # from 0 to 1, i.e. from black to white, dz times:
                #  [0 .. 1] -> [0 0 ... 0 1 1 ... 1]
                # This code makes a wider outer black area: ----
                greys = np.zeros(360)
                greys[:340] = np.linspace(1, 0, 340)
                greys[340:] = 0
                rgb[:, 1] = np.repeat(greys, dz)
                # Instead of: -----
                # rgb[:, 1] = np.repeat(np.linspace(0, 1, 360), dz)
                # -----
                # Now we convert every row in HLS to RGB values
                rgb = np.apply_along_axis(convert_to_RGB, 1, rgb)

                # And plot in the polar axes:
                axColor.pcolormesh(azimuths * np.pi / 180.0, zeniths,
                                   # only necessary as required n of args:
                                   np.zeros((dz, 360)),
                                   # cmap=plt.cm.hsv
                                   color=rgb
                                   )
                axColor.set_yticks([])
                # axColor.set_xticks([0, np.pi * 0.5, np.pi, 1.5 * np.pi])
                axColor.set_thetagrids([0, 90, 180, 270])
                axColor.tick_params(axis='x', pad=0)
                axColor.set_xticklabels([r'$0$', r'$\pi/2$',
                                         r'$\pi$', r'$3\pi/2$'],
                                        # fontsize=18
                                        )
                axColor.text(0.5, 0.5, r'$\vec{m}$',
                             horizontalalignment='center',
                             verticalalignment='center',
                             transform=axColor.transAxes,
                             # fontsize=20
                             )

        else:
            # Spin data in a grid
            spin_z = self.data[:, 2].reshape(-1, self.nx)

            ax.imshow(spin_z, origin='lower', cmap=cmap, interpolation='None',
                      vmin=-1, vmax=1,
                      extent=[self.xmin * 1e9, self.xmax * 1e9,
                              self.ymin * 1e9, self.ymax * 1e9]
                      )

            if cbar:
                plt.colorbar()

        ax.set_ylabel(r'y (nm)')
        ax.set_xlabel(r'x (nm)')

        if savefig:
            plt.savefig(savefig, bbox_inches='tight', dpi=dpi)
        # plt.show()

    def plot_charge_density(self, savefig=None, dpi=150):

        f = plt.figure()

        charge = self.sk_number
        vmax = np.max(np.abs(charge))

        charge.reshape(-1,)[self.data_norm < 1e-10] = np.nan

        plt.imshow(charge, origin='lower', cmap='RdYlBu', interpolation='None',
                   vmin=-vmax, vmax=vmax,
                   extent=[self.xmin * 1e9, self.xmax * 1e9,
                           self.ymin * 1e9, self.ymax * 1e9]
                   )
        plt.ylabel(r'y (nm)')
        plt.xlabel(r'x (nm)')
        plt.colorbar()

        if savefig:
            plt.savefig(savefig, bbox_inches='tight', dpi=dpi)
        # plt.show()


if __name__ == '__main__':
    _file = 'test/isolated_sk_PdFeIr_PBC-Oxs_TimeDriver-Magnetization-0000000-0000000.omf'
    # _file = 'files/IrCo_J5e11_eprime_012.omf'
    # ngbs = generate_ngbs(h)
    # print(ngbs)
    oommf_data = SkNumberOOMMF(_file)
    print(oommf_data.compute_sk_number())
    oommf_data.plot_system()
