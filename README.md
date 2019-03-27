[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1296536.svg)](https://doi.org/10.5281/zenodo.1296536)

# WARNING:

This tool is not developed anymore and the project has been migrated to OOMMFPy:

https://github.com/davidcortesortuno/oommfpy



# OOMMF Skyrmion Number

Python library for the calculation of the skyrmion number
or topological charge from OMF files produced by OOMMF.

## Instructions 

The OMF files have to be in text format and regular grids. You can obtain
this format directly from OOMMF using the OOMMF tools:

    tclsh oommf.tcl avf2ovf -grid reg my_file.omf my_file_in_text.omf

The format is not important, can be `ovf`, `omf`, etc.

The Python library was designed to be used in external scripts that
can process multiple files or folders. Examples are in the `script` folder.
The main Class for loading and processing the files is called
`SkNumberOOMMF`.

All the calculations are performed using vectorised operations with the Numpy
library.

For now, the Q number can only be computed in a slice of the system in
the XY plane. Works well for 2D systems.

## Basic usage

To test the library, we can use an interactive console as IPython to
analyse a particular `omf` file, say, `my_file.omf`:

```python

import oommf_sk_number as oskn

# Load the class. z_index indicate the layer number
oommf_data = oskn.SkNumberOOMMF('my_file.omf', z_index=0)
print('Sk number =', oommf_data.compute_sk_number())

# Plot the charge density
oommf_data.plot_charge_density('muy_file_charge.png')
```

Options can be passed to the plot methods.

## Version

This is a very first version of the library. There are many options that still
need to be implemented/polished.

## TODO


- [ ] Multiple tests
- [ ] More documentation
- [ ] Create proper system wide installation (setup.py)
- [ ] Check labels on plots
- [ ] Allow slices of the system in different directions
- [ ] Calculate coordinates
- [ ] Allow Periodic boundaries
- [ ] Print z coordinate specified by the z-index

# Citation

Cite this repository as:

```
@Misc{Cortes2018,
  author       = {David Cort{\'e}s-Ortu{\~n}o},
  title        = {OOMMF Skyrmion Number},
  howpublished = {Zenodo doi:10.5281/zenodo.1296536. Github: https://github.com/davidcortesortuno/oommf_sk_number},
  year         = {2018},
  doi          = {10.5281/zenodo.1296536},
  url          = {https://doi.org/10.5281/zenodo.1296536},
}

```
