# This makefile compiles copygb
# NOTE: This is a serial build as the code is written to be serial
# Usage: copygb [-g "grid [kgds]"] [-i "ip [ipopts]"] [-k "kpds"] [-v "uparms"]
#              [-B mapgrib [-b mapindex] [-A "<> mapthreshold"] [-K "mapkpds"]]
#              [-M "mask"/mergegrib [-m mergeindex]] [-X] [-a] [-N namelist]
#       then either:
#              grib1 index1 grib2
#            or:
#              -x grib1 grib2
#
# Perform destaggering (NMM) and HORIZONTAL interpolation to a defined output grid.
# Create output grid NOT fixed by model integration domain

SHELL    = /bin/sh
TARGET   = copygb.exe

#
# build configuration based on architecture and compiler
include ../../configure.upp

#
# directories of shared resources
MODULES  = ../../sorc/ncep_post.fd.new/kinds_mod.o
LLIBDIR      =    -L/scratch/lus/dev/intel/lib -L/scratch/lus/arw/model/src/post/intel/UPPV4.1/sorc/comlibs/xml
NCEPLIBINC      =    -I/scratch/lus/dev/intel/include
UPPLIBS   =    -lwrfio -lg2_v3.1.0_4 -lg2tmpl_v1.5.0 -lnemsio_d -lsigio_v2.0.1_4 -lsfcio_v1.0.0_4 -lgfsio_4 -lsp_v2.0.2_d -lw3nco_v2.0.6_4 -lw3emc_v2.2.0_4 -lbacio_4 -lip_d -lsp_v2.0.2_d -L/scratch/lus/dev/lib  -ljasper -lpng -lcurl
LOCALINC = -I$(INCMOD) $(NCEPLIBINC)

LIBS     = $(LLIBDIR) $(UPPLIBS)



#
# extra compilation/link flags
#   we are using single word (4 byte) for all libraries
EXTRA_FFLAGS   =
#EXTRA_FFLAGS   = -c -r8
EXTRA_LDFLAGS  = $(LIBS) 

SRCS=	copygb.f
OBJS=   $(SRCS:.f=.o)

#
# TARGETS

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(TARGET)
	$(CP) $(TARGET) $(BINDIR)

$(TARGET): $(SRCS) $(MODULES)
	$(SFC) -o $@ $(SRCS) $(MODULES) $(FFLAGS) $(EXTRA_FFLAGS) $(LDFLAGS) $(EXTRA_LDFLAGS)

clean:	
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	$(RM) $(TARGET) $(OBJS)
	$(RM) $(BINDIR)/$(TARGET)

.IGNORE:
.PHONY: clean
