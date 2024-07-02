############################ -*- Mode: Makefile -*- ###########################
## Makefile --- 
## 
## Auteur          : Pierre-François LAVALLEE (CNRS/IDRIS) <lavallee@idris.fr>
## Créé le         : Tue Oct  7 10:54:47 2003
## Dern. mod. par  : Pierre-François LAVALLEE (CNRS/IDRIS) <lavallee@idris.fr>
## Dern. mod. le   : Mon Oct 11 16:40:26 2010
###############################################################################

# Link make.inc with the targeted architecture in Arch directory
# For example : ln -s make_sp6 make.inc
#--------------------------------------
include ./make.inc

FFLAGS += $(INCS)

FFLAGS += #-DADIOS -DHDF5 -DMPIIO -DNETCDF4 -DPNETCDF -DSIONLIB

# Sources directories
#--------------------
SRCDIR = .

# Objects
#------------
OBJ = Ptim.o module_hydro_commun.o mpi_module.o module_hydro_utils.o module_adios.o module_hdf5.o module_mpiio.o module_netcdf.o module_pnetcdf.o module_sionlib.o module_hydro_IO.o module_hydro_principal.o main.o

hydro:  $(OBJ)
	$(MPIF90) $(OBJ) $(LDFLAGS) -o hydro 

clean:
	rm -f *.o *.mod 

cleanall:
	rm -f *.o *.mod hydro

%.o:    $(SRCDIR)/%.f90
	$(MPIF90) $(FFLAGS) -c $<

# Dependencies
#-------------
main.o : main.f90 Ptim.o module_hydro_principal.o mpi_module.o module_hydro_IO.o module_hydro_commun.o
module_hydro_commun.o : module_hydro_commun.f90
module_hydro_IO.o : module_hydro_IO.f90 module_adios.f90 module_hdf5.f90 module_mpiio.o module_netcdf.f90 module_pnetcdf.f90 module_sionlib.f90 module_hydro_commun.o module_hydro_utils.o mpi_module.o module_mpiio.o
module_hydro_principal.o : module_hydro_principal.f90 module_hydro_commun.o module_hydro_utils.o mpi_module.o
module_hydro_utils.o : module_hydro_utils.f90 module_hydro_commun.o mpi_module.o
mpi_module.o : mpi_module.f90 module_hydro_commun.o
Ptim.o : Ptim.f90

module_adios.o : module_adios.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o
module_hdf5.o : module_hdf5.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o
module_mpiio.o : module_mpiio.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o
module_netcdf.o : module_netcdf.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o
module_pnetcdf.o : module_pnetcdf.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o
module_sionlib.o : module_sionlib.f90 mpi_module.o module_hydro_commun.o module_hydro_utils.o

