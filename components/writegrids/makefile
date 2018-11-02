ifeq ($(NONETCDF),TRUE)
        SRCSF = src/writegrids_stub.F90
        NETCDF_INC=
else
        SRCSF = src/writegrids.F90
        NETCDF_INC=-I $(NETCDF_DIR)/include
endif

BUILDDIR=build
COREDIR=../../model_core/build
FFLAGS=-I $(BUILDDIR) -I $(COREDIR) $(NETCDF_INC) $(COMPILERFFLAGS)
OBJS = $(patsubst %.F90,$(BUILDDIR)/%.o,$(SRCSF))

all: create-build-dirs $(OBJS)

create-build-dirs:
	mkdir -p $(BUILDDIR)

$(OBJS) : $(BUILDDIR)/%.o : %.F90
	$(FTN) $(OPT) $(FFLAGS) $< -o $(BUILDDIR)/$(notdir $@)
