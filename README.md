# Photoelectron_refraction
This package examines the refraction effect in the three-step photoemission model.
It has Python and C++ codes; Python version is available for any OS in which Python works, but C++ version works on computers with a C++ compiler such as Linux system.

Since C++ version is much faster than Python version, we recommend to use C++ version for calculations with many (>1000) randomized surfaces.

# Usage

## Python version (/Python)
Before using this package, please copy ```Config_example.py``` to ```Config.py``` and edit if necessary.
```Config.py``` can include personal settings for GUI and is out of Git control.

The following packages are necessary, which can be installed via ```pip install``` command.
- PyQt5
- pyqtgraph
- numpy
- h5py

Execute the following command to open GUI and perform calculations.
```
>python Calculation_GUI.py
```

## C++ version (/CPP)
You need to install HDF5 package from https://www.hdfgroup.org/downloads/hdf5 .
Please don't forget to include ```--enable-cxx``` option when executing ```configure```.
After installation, please modify ```Makefile``` adequately and compile the package by ```make``` command.
The package is based on OpenMP parallelization, so if you don't use OpenMP you need to modify the program.

When compilation has finished successfully, you can perform calculations by
```
>Refraction.o input.dat
```
The description of the input file is writtein in the document above.
