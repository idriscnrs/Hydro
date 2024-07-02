# C compiler
#------------------
CC = pgcc

# Fortran compilers
#------------------
F90 = pgfortran

# MPI C compiler
#------------------
MPICC = mpicc

# MPI Fortran compiler
#------------------
MPIF90 = mpif90

# C compilation options
#-----------------------------
CFLAGS = -O0 -g

# Fortran compilation options
#-----------------------------
FFLAGS = -O0 -g -acc -m64 -ta=tesla:cc70 -Minfo=accel -Mpreprocess -Dkey_riemann_oldalg

# Linker options
#---------------
LDFLAGS = -O0 -g -acc -ta=tesla:cc70 -Mpreprocess

# OpenMP C compilation options
#---------------------------
CFLAGS_OMP =  -fopenmp

# OpenMP Fortran compilation options
#---------------------------
FFLAGS_OMP =  -fopenmp

# MPI C compilation options
#---------------------------
CFLAGS_MPI =  -m64

# MPI Fortran compilation options
#---------------------------
FFLAGS_MPI =  -m64

# OpenMP link options
#--------------------
LDFLAGS_OMP = -fopenmp

# MPI link options
#--------------------
LDFLAGS_MPI = -m64

