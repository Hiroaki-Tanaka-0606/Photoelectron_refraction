# Photoelectron_refraction
# Makefile

# directory specified by --prefix=DIR in ./configure of HDF5
# written in config.log
# hdf5-x.y.z/hdf5 by default
hdf5 = /home/hiroaki/hdf5-1.12.1/hdf5

# Compile option of HDF5 can be checked by
# >$(hdf5)/bin/h5c++ (filename) -show
# >$(hdf5)/bin/h5c++ (filename) -c -show

commonOption = -DOLD_HEADER_FILENAME -DHDF_NO_NAMESPACE -DNO_STATIC_CAST

# $(hdf5)/include containes "H5Cpp.h"
compileOption = $(commonOption) -I $(hdf5)/include

# In $(hdf5)/lib are HDF5 libraries
linkOption1 = -L $(hdf5)/lib

# Specify necessary libraries for HDF5
linkOption2 = -l hdf5_hl_cpp -l hdf5_cpp -l hdf5_hl -l hdf5 -Wl,-rpath,$(hdf5)/lib

# OpenMP configuration
MKLROOT = /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin
option = -L${MKLROOT} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl $(linkOption1) $(linkOption2)

target = Refraction_main.o HDF5_tools.o

Refraction.o: $(target)
	$(CXX) $(target) $(option) -o Refraction.o

.cpp.o:
	$(CXX) -fopenmp $(compileOption) -c $<

