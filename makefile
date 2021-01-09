# This is a convenience thing, on ARCHER it will default to GNU whereas locally to local (different mpi wrappers are used)
ifdef CRAYOS_VERSION
.DEFAULT_GOAL :=GNU
else
.DEFAULT_GOAL :=local
endif

CORE_DIR=model_core
COMPONENTS_DIR=components
TESTCASE_DIR=testcases
IO_SERVER_DIR=io
BUILD_DIR=build


ifndef NETCDF_ROOT
	export NETCDF_ROOT=$(shell nf-config --prefix)
endif
ifeq ($(NETCDF_ROOT),)
	export NETCDF_ROOT=$(shell nf-config --prefix)
endif
ifndef NETCDF_DIR
	export NETCDF_DIR=$(NETCDF_ROOT)
endif
ifndef HDF5_DIR
	export HDF5_DIR=$(HDF5_ROOT)
endif

ifndef CASIM_DIR
	export CASIM_DIR=$(CASIM_ROOT)
endif
ifndef PETSC_DIR
	export PETSC_DIR=$(PETSC_ROOT)
endif

FTN=ftn

ifeq ($(strip $(CASIM_DIR)),)
	CASIM_BUILD_LOC=
else
	CASIM_BUILD_LOC=$(CASIM_DIR)/$(BUILD_DIR)/*.o
endif

ifeq ($(strip $(PETSC_DIR)),)
	PETSC_LIB_PATH=
	PETSC_LIBS=
else
	PETSC_LIB_PATH=-L $(PETSC_DIR)/lib
	PETSC_LIBS=-lpetsc
endif

COMPILERFFLAGS=-fopenmp -O3 #
COMPILERRECURSIVE=
ACTIVE=-DU_ACTIVE -DV_ACTIVE -DW_ACTIVE -DUSE_MAKE
DEBUG_FLAGS=-g -fopenmp -fcheck=all -fbacktrace -DDEBUG_MODE -ffpe-trap=invalid,overflow,zero

FFLAGS=-I $(CORE_DIR)/$(BUILD_DIR) -I $(COMPONENTS_DIR)/$(BUILD_DIR) -I $(TESTCASE_DIR)/$(BUILD_DIR) -I $(NETCDF_DIR)/include -I /usr/include $(COMPILERFFLAGS)
LFLAGS=$(shell $(NETCDF_ROOT)/bin/nf-config --flibs)
EXEC_NAME=monc

local: FTN=mpif90
local: GNU

debug: COMPILERFFLAGS = $(DEBUG_FLAGS)
debug: OPT=$(ACTIVE)
ifdef CRAYOS_VERSION
debug: GNU
else
debug: local
endif

nonetcdf: FFLAGS =-I $(CORE_DIR)/$(BUILD_DIR) -I $(COMPONENTS_DIR)/$(BUILD_DIR) -I $(TESTCASE_DIR)/$(BUILD_DIR) -I /usr/include $(COMPILERFFLAGS)
nonetcdf: LFLAGS=
nonetcdf: NONETCDF=TRUE
ifdef CRAYOS_VERSION
nonetcdf: GNU
else
nonetcdf: local
endif

nonetcdf-debug: COMPILERFFLAGS = $(DEBUG_FLAGS)
nonetcdf-debug: OPT=$(ACTIVE)
nonetcdf-debug: FFLAGS =-I $(CORE_DIR)/$(BUILD_DIR) -I $(COMPONENTS_DIR)/$(BUILD_DIR) -I $(TESTCASE_DIR)/$(BUILD_DIR) -I /usr/include $(COMPILERFFLAGS)
nonetcdf-debug: LFLAGS= 
nonetcdf-debug: NONETCDF=TRUE
ifdef CRAYOS_VERSION
nonetcdf-debug: GNU
else
nonetcdf-debug: local
endif

GNU: COMPILERFFLAGS += $(ACTIVE) -cpp -J $(BUILD_DIR) -c
GNU: COMPILERRECURSIVE= -frecursive
GNU: buildmonc

Cray: COMPILERFFLAGS += $(ACTIVE) -e m -J $(BUILD_DIR) -c
Cray: COMPILERRECURSIVE= -e R
Cray: buildmonc

Intel: COMPILERFFLAGS += $(ACTIVE) -fpp -free -c -std03 -module $(BUILD_DIR) -Tf
Intel: COMPILERRECURSIVE= -recursive
Intel: buildmonc

IBM: FTN=mpxlf2003_r
IBM: ACTIVE:= $(foreach option,$(ACTIVE),-WF,$(option))
IBM: COMPILERFFLAGS = $(ACTIVE) -cpp -qmoddir=$(BUILD_DIR) -c
IBM: COMPILERRECURSIVE= -qrecur
IBM: buildmonc

export COMPILERFFLAGS
export COMPILERRECURSIVE
export FTN
export NONETCDF

clean: clean-model_core clean-components clean-testcases
	rm -Rf build/* fort.10 fort.11

clean-build: clean-build-model_core clean-build-components clean-build-testcases
	rm -Rf build fort.10 fort.11

buildmonc: check-vars create-build-dirs compile-model_core compile-components compile-testcases compile-bootstrapper
	$(FTN) -fopenmp -o $(EXEC_NAME) $(BUILD_DIR)/*.o $(CORE_DIR)/$(BUILD_DIR)/*.o $(COMPONENTS_DIR)/$(BUILD_DIR)/*.o $(CASIM_BUILD_LOC)  $(TESTCASE_DIR)/$(BUILD_DIR)/*.o $(LFLAGS)

check-vars:
	$(call check_defined, NETCDF_DIR, Need the path to the NetCDF installation directory as an environment variable - export this before running make)
	$(call check_defined, HDF5_DIR, Need the path to the HDF5 installation directory as an environment variable - export this before running make)

create-build-dirs:
	mkdir -p $(BUILD_DIR)

compile-model_core:
	cd $(CORE_DIR) ; $(MAKE)

clean-model_core:
	cd $(CORE_DIR); $(MAKE) clean

clean-build-model_core:
	cd $(CORE_DIR); $(MAKE) clean-build

compile-components:
	cd $(COMPONENTS_DIR) ; $(MAKE)

clean-components:
	cd $(COMPONENTS_DIR); $(MAKE) clean

clean-build-components:
	cd $(COMPONENTS_DIR); $(MAKE) clean-build

compile-testcases:
	cd $(TESTCASE_DIR) ; $(MAKE)

clean-testcases:
	cd $(TESTCASE_DIR); $(MAKE) clean

clean-build-testcases:
	cd $(TESTCASE_DIR); $(MAKE) clean-build

compile-ioserver:
	cd $(IO_SERVER_DIR) ; $(MAKE)

clean-ioserver:
	cd $(IO_SERVER_DIR) ; $(MAKE) clean

clean-build-ioserver:
	cd $(IO_SERVER_DIR) ; $(MAKE) clean-build

clean-outputs:
	rm -rf parcels_???_????.dat grids_???_????.dat

compile-bootstrapper:
	$(FTN) $(FFLAGS) monc_driver.F90 -o $(BUILD_DIR)/monc_driver.o

# Check that given variables are set and all have non-empty values,
# die with an error otherwise.
#
# Params:
#   1. Variable name(s) to test.
#   2. (optional) Error message to print.
check_defined = \
    $(foreach 1,$1,$(__check_defined))
__check_defined = \
    $(if $(value $1),, \
      $(error Undefined $1$(if $(value 2), ($(strip $2)))))
