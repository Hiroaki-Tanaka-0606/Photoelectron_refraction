hdf5 = /home/hiroaki/hdf5-1.12.1/hdf5

commonOption = -DOLD_HEADER_FILENAME -DHDF_NO_NAMESPACE -DNO_STATIC_CAST

compileOption = $(commonOption) -I $(hdf5)/include

linkOption1 = -L $(hdf5)/lib

linkOption2 = -l hdf5_hl_cpp -l hdf5_cpp -l hdf5_hl -l hdf5 -Wl,-rpath,$(hdf5)/lib

MKLROOT = /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/lib/intel64_lin
option = -L${MKLROOT} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl $(linkOption1) $(linkOption2)

target = Refraction_main.o HDF5_tools.o

Refraction.o: $(target)
	$(CXX) $(target) $(option) -o Refraction.o

.cpp.o:
	$(CXX) -fopenmp $(compileOption) -c $<

