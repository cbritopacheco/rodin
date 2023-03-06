# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Variables corresponding to defines in config.hpp (YES, NO, or value)
MFEM_VERSION           = 40501
MFEM_VERSION_STRING    = 4.5.1
MFEM_SOURCE_DIR        = /Users/carlos/Projects/rodin/third-party/mfem
MFEM_INSTALL_DIR       = /usr/local
MFEM_GIT_STRING        = remotes/origin/HEAD-0-g1e17003af3304e1b74983d4d44dc83bf2e4ca138
MFEM_USE_MPI           = NO
MFEM_USE_METIS         = NO
MFEM_USE_METIS_5       = NO
MFEM_DEBUG             = YES
MFEM_USE_EXCEPTIONS    = NO
MFEM_USE_ZLIB          = NO
MFEM_USE_LIBUNWIND     = NO
MFEM_USE_LAPACK        = NO
MFEM_THREAD_SAFE       = YES
MFEM_USE_LEGACY_OPENMP = NO
MFEM_USE_OPENMP        = NO
MFEM_USE_MEMALLOC      = YES
MFEM_TIMER_TYPE        = 4
MFEM_USE_SUNDIALS      = NO
MFEM_USE_MESQUITE      = NO
MFEM_USE_SUITESPARSE   = YES
MFEM_USE_SUPERLU       = NO
MFEM_USE_SUPERLU5      = NO
MFEM_USE_MUMPS         = NO
MFEM_USE_STRUMPACK     = NO
MFEM_USE_GINKGO        = NO
MFEM_USE_AMGX          = NO
MFEM_USE_GNUTLS        = NO
MFEM_USE_NETCDF        = NO
MFEM_USE_PETSC         = NO
MFEM_USE_SLEPC         = NO
MFEM_USE_MPFR          = NO
MFEM_USE_SIDRE         = NO
MFEM_USE_FMS           = NO
MFEM_USE_CONDUIT       = NO
MFEM_USE_PUMI          = NO
MFEM_USE_HIOP          = NO
MFEM_USE_GSLIB         = NO
MFEM_USE_CUDA          = NO
MFEM_USE_HIP           = NO
MFEM_USE_RAJA          = NO
MFEM_USE_OCCA          = NO
MFEM_USE_CEED          = NO
MFEM_USE_CALIPER       = NO
MFEM_USE_UMPIRE        = NO
MFEM_USE_SIMD          = NO
MFEM_USE_ADIOS2        = NO
MFEM_USE_MKL_CPARDISO  = NO
MFEM_USE_MOONOLITH     = NO
MFEM_USE_ADFORWARD     = NO
MFEM_USE_CODIPACK      = NO
MFEM_USE_BENCHMARK     = NO
MFEM_USE_PARELAG       = NO
MFEM_USE_ENZYME        = NO

# Compiler, compile options, and link options
MFEM_CXX       = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
MFEM_HOST_CXX  = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
MFEM_CPPFLAGS  = 
MFEM_CXXFLAGS  = -std=c++11 -g  -g --coverage
MFEM_TPLFLAGS  = 
MFEM_INCFLAGS  = -I$(MFEM_INC_DIR) $(MFEM_TPLFLAGS)
MFEM_PICFLAG   = 
MFEM_FLAGS     = $(MFEM_CPPFLAGS) $(MFEM_CXXFLAGS) $(MFEM_INCFLAGS)
MFEM_EXT_LIBS  =  -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lmetis -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lumfpack -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lklu -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lamd -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lbtf -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lcholmod -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lcolamd -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lcamd -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lccolamd -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lsuitesparseconfig -Wl,-rpath,/opt/local/lib -L/opt/local/lib -lopenblas
MFEM_LIBS      = -L$(MFEM_LIB_DIR) -lmfem $(MFEM_EXT_LIBS)
MFEM_LIB_FILE  = $(MFEM_LIB_DIR)/libmfem.a
MFEM_STATIC    = YES
MFEM_SHARED    = NO
MFEM_BUILD_TAG = Darwin-21.3.0
MFEM_PREFIX    = /usr/local
MFEM_INC_DIR   = /Users/carlos/Projects/rodin/build2/third-party/mfem
MFEM_LIB_DIR   = /Users/carlos/Projects/rodin/build2/third-party/mfem

# Location of test.mk
MFEM_TEST_MK = /Users/carlos/Projects/rodin/third-party/mfem/config/test.mk

# Command used to launch MPI jobs
MFEM_MPIEXEC    = mpirun
MFEM_MPIEXEC_NP = -np
MFEM_MPI_NP     = 4

# The NVCC compiler cannot link with -x=cu
MFEM_LINK_FLAGS := $(filter-out -x=cu -xhip, $(MFEM_FLAGS))

# Optional extra configuration
MFEM_BUILD_DIR ?= /Users/carlos/Projects/rodin/build2/third-party/mfem
