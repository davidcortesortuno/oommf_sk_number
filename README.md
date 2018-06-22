# OOMMF Skyrmion Number

Python library for the calculation of the skyrmion number
or topological charge from OMF files produced by OOMMF.

## Instructions 

The OMF files have to be in text format and irregular grid. You can obtain
this format directly from OOMMF using the OOMMF tools:

    tclsh oommf.tcl avf2ovf -grid irreg my_file.omf my_file_in_text.omf

The format is not important, can `ovf`, `omf`, etc.

The Python library was designed to be used in external scripts that
can process multiple files or folders. Examples are in the `script` folder.
The main Class for loading and processing the files is called
`SkNumberOOMMF`.

All the calculations are performed using vectorised operations with the Numpy
library.

## Basic usage

To test the library, we can use an interactive console as IPython to
analyse a particular `omf` file, say, `my_file.omf`:

```python

import oommf_sk_number as oskn

# Load the class
oommf_data = oskn.SkNumberOOMMF('my_file.omf')
print('Sk number =', oommf_data.compute_sk_number())

# Plot the charge density
oommf_data.plot_charge_density('muy_file_charge.png')
```

Options can be passed to the plot methods.

## Version

This is a very first version of the library. There are many options that still
need to be implemented/polished.
