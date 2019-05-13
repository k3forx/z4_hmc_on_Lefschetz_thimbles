FC=gfortran -O3
config=`python3-config --ldflags`

all : generate_thimble

generate_thimble : forpy_mod.F90 generate_thimble.F90
#	$(FC) -c forpy_mod.F90
	$(FC) generate_thimble.F90 forpy_mod.o $(config) -o $@

# cooling : thimble_cooling_z4.f90
# 	$(FC) $^ -o $@

# select : thimble_select.f90
# 	$(FC) $^ -o $@

# MC_thimble : mt19937.f90 MC_thimble_z4.f90
# 	$(FC) $^ -o $@ -I/usr/local/include -llapack95 -llapack -lblas -fno-range-check

# histogram : histogram2d.f90
# 	$(FC) $^ -o $@

# calculation : estimate.f90
# 	$(FC) $^ -o $@

clean :
	rm -f *.o *.mod *~ generate_thimble
