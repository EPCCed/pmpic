SRC_COMMON= src/datadefn.F90
SRC_SPECIAL=src/utils/conversions.F90
SRCS_UTIL= src/utils/fileunits.F90 src/utils/logging.F90 src/utils/collections.F90 src/utils/maths.F90 src/utils/q_indices.F90 src/utils/naming_conventions.F90 src/utils/string_utils.F90
SRCS_COMPONENTS=src/components/monc_component.F90 src/components/registry.F90 src/components/timestepper.F90
SRCSF = $(SRCS_UTIL) src/grid/grids.F90 src/grid/prognostics.F90 src/grid/stencil.F90 src/grid/interpolation.F90 src/communication/communicationtypes.F90 src/parcels/parcels.F90 src/state.F90 src/communication/halocommunication.F90 src/physics/saturation.F90 src/configuration/optionsdatabase.F90 src/utils/timer.F90  $(SRCS_COMPONENTS) src/configuration/configurationfileparser.F90 src/physics/scienceconstants.F90 src/communication/MPIC_Haloswap.F90 src/parcels/parcel_interpolation.F90 src/parcels/haloswap.F90 src/fft/ffte_mod.F90 src/fft/pencilfft.F90 src/fft/fft_ops.F90 src/monc.F90
SRCSF77 = src/fft/ffte/factor.f src/fft/ffte/fft235.f src/fft/ffte/kernel.f src/fft/ffte/mfft235.f src/fft/ffte/vzfft1d.f

BUILDDIR=build
IO_BUILDDIR=../io/build

FFLAGS=-I $(BUILDDIR) -I $(NETCDF_DIR)/include $(COMPILERRECURSIVE) $(COMPILERFFLAGS)
FFLAGS=-I $(BUILDDIR)  $(COMPILERRECURSIVE) $(COMPILERFFLAGS) -I/usr/include

F77OBJS = $(patsubst %.f,$(BUILDDIR)/%.o,$(SRCSF77))
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))
OBJS_SPECIAL= $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRC_SPECIAL))
OBJS_COMMON= $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRC_COMMON))
FFLAGS_SPECIAL=$(filter-out -fcheck=%,$(FFLAGS))

all: create-build-dirs $(OBJS_COMMON) $(OBJS_SPECIAL) $(OBJS) $(F77OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS_SPECIAL) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(FFLAGS_SPECIAL) $< -o $(BUILDDIR)/$(notdir $@)

$(OBJS_COMMON) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)

$(F77OBJS) : $(BUILDDIR)/%.o : %.f
	$(FTN) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)

clean:
	rm -Rf $(BUILDDIR)/*

clean-build:
	rm -Rf $(BUILDDIR)
