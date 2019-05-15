FC=gfortran -O3
config=`python3-config --ldflags`

all : generate_thimble hmc_on_thimble

generate_thimble : forpy_mod.F90 generate_thimble.F90
	$(FC) -c forpy_mod.F90
	$(FC) generate_thimble.F90 forpy_mod.o $(config) -o $@

hmc_on_thimble : mt19937.f90 hmc_on_thimble.F90
	$(FC) $^ -o $@ -fno-range-check

clean :
	rm -f *.o *.mod *~ generate_thimble hmc_on_thimble
