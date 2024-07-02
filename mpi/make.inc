# OpenACC options
#------------------
use_rdma		= ok
use_riemann_oldalg	= ok
use_riemann_bigloop	= nok

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
MPIF90 = mpifort

# C compilation options
#-----------------------------
CFLAGS = -O0 -g

# Fortran compilation options
#-----------------------------
FFLAGS = -O0 -g -Mnoidiom -m64 -acc -ta=tesla:cc70 -Minfo=accel -Mpreprocess -Minline -Dkey_gpu

ifeq ($(use_rdma), ok)
 FFLAGS	+= -Dkey_rdma
endif
ifeq ($(use_riemann_oldalg), ok)
 FFLAGS	+= -Dkey_riemann_oldalg
endif
ifeq ($(use_riemann_bigloop), ok)
 FFLAGS	+= -Dkey_riemann_bigloop
endif

# Linker options
#---------------
LDFLAGS = -O0 -g -acc -ta=tesla:cc70 -m64 -Mpreprocess -Minline

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

