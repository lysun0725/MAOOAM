#ifndef COMPILER
  COMPILER=gfortran
#endif
ifdef INTMODE
  # Define integrator module
  INTEGRATOR = $(INTMODE)_integrator
  TL_AD_INTEGRATOR = $(INTMODE)_tl_ad_integrator
else
  # Define integrator module
  INTEGRATOR = rk2_integrator
  TL_AD_INTEGRATOR = rk2_tl_ad_integrator
endif
ifdef RES
  ST=_$(RES)
  $(warning Code will $(RES) the inner products in arrays! It may induce a huge memory usage!)
endif
# Define relevant compiler options for debug case and normal case for each compiler
# Normal Case:
IFORT_COMPILER_FLAGS = -assume byterecl -O2
IFORT_LD_FLAGS = -mkl
GCC_COMPILER_FLAGS = -O2 -Wall
#GCC_LD_FLAGS = -L$(LAPACKLIB) -llapack -lblas
GCC_LD_FLAGS = -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64
# Debug Case:
IFORT_COMPILER_FLAGS_DBG = -O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl
IFORT_LD_FLAGS_DBG = -llapack -lblas
GCC_COMPILER_FLAGS_DBG = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow 
#GCC_LD_FLAGS_DBG = -llapack -lblas
GCC_LD_FLAGS_DBG = -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -lm -lmkl_core -lmkl_lapack95_lp64 


SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif


ifeq ($(COMPILER),gfortran)
  #Derive gfortran version from gcc version and check if >=4.6.0.
  GVERSION := $(shell gcc -dumpversion)
  GE40600 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40600)
  
  ifeq ($(COMPILER),gfortran)
    ifeq ($(GE40600), 0)
      $(error GCC/gfortran version $(GVERSION) detected. This program uses features of gfortran >= 4.6.0. Please update your toolchain)
    endif
  else
    $(warning This code was optimized for gfortran >= 4.6.0 and has not been tested with other compilers)
  endif
  COMPILER_FLAGS = $(GCC_COMPILER_FLAGS)
  LD_FLAGS = $(GCC_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(GCC_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(GCC_LD_FLAGS_DBG)
  IFPORT = ifport.mod
else ifeq ($(COMPILER),ifort)
  $(warning This code was optimized for ifort >= 14.0.2 and has not been tested with other compilers)
  COMPILER_FLAGS = $(IFORT_COMPILER_FLAGS)
  LD_FLAGS = $(IFORT_LD_FLAGS)
  COMPILER_FLAGS_DBG = $(IFORT_COMPILER_FLAGS_DBG)
  LD_FLAGS_DBG = $(IFORT_LD_FLAGS_DBG)
else
  $(warning This code was optimized for gfortran >= 4.6.0 and ifort >= 14.0.2 and has not been tested with other compilers)
endif

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .out .test

PROGRAMS = maooam test_tl_ad nature_obs_maooam test_dr #freerun_maooam etkf_maooam wc_etkf_maooam fcst_leading_time fcst_leading_time_wc 

MODULES = m_mt.mod m_da_maooam.mod util.mod params.mod inprod_analytic$(ST).mod aotensor_def.mod ic_def.mod $(INTEGRATOR).mod rk2_dr_integrator.mod stat.mod tensor.mod tl_ad_tensor.mod $(TL_AD_INTEGRATOR).mod $(IFPORT)
MODULE_OBJECTS = $(MODULES:.mod=.o)

TEST_PROGRAMS = test_inprod_analytic test_aotensor
TESTS = $(TEST_PROGRAMS:=.test)

all: $(PROGRAMS) 

debug: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
debug: LD_FLAGS = $(LD_FLAGS_DBG)
debug: all

test: COMPILER_FLAGS = $(COMPILER_FLAGS_DBG) 
test: LD_FLAGS = $(LD_FLAGS_DBG)
test: $(TESTS)

test_inprod_analytic.test: test_inprod_analytic
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

test_aotensor.test: test_aotensor 
	./$< |sort  > tests/$<.out
	diff tests/$<.out tests/$<.ref && echo PASS: $@ || echo FAIL: $@

%.mod: %.f90
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

%.o: %.f90 $(MODULES)
	$(COMPILER) $(COMPILER_FLAGS) -c $< $(LD_FLAGS) 

#%.o: %.f 
#	$(COMPILER) $(LD_FLAGS) $(COMPILER_FLAGS) -c $<

%: %.o $(MODULE_OBJECTS)
	$(COMPILER) $(COMPILER_FLAGS) $^ -o $@ $(LD_FLAGS) 

inprod_analytic$(ST).mod: params.mod util.mod
aotensor_def.mod: tensor.mod inprod_analytic$(ST).mod
$(INTEGRATOR).mod: aotensor_def.mod
rk2_dr_integrator.mod: aotensor_def.mod
ic_def.mod: util.mod aotensor_def.mod
stat.mod: params.mod
tl_ad_tensor.mod: aotensor_def.mod
$(TL_AD_INTEGRATOR).mod: tl_ad_tensor.mod
util.mod: $(IFPORT)
m_da_maooam.mod: m_da_maooam.mod
m_mt.mod: m_mt.mod

clean:
	$(RM) *.o *.mod $(PROGRAMS) $(TEST_PROGRAMS) tests/*.out

.PHONY: clean all test %.test
