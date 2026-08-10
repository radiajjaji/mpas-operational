.SUFFIXES: .F .o

all: dummy schemes

dummy:
	echo "****** compile chemistry_schemes ******"

OBJS = \
	module_data_binaer_other.o           \
	module_data_binaer_asect.o           \
	module_data_mosaic_therm.o           \
	module_data_gocart_dust.o            \
	module_data_sorgam.o                 \
	module_data_radm2.o                  \
	module_data_rrtmgaeropt.o            \
	module_mosaic_addemiss.o             \
	module_cal_cldfra.o                  \
	module_gocart_drydep.o               \
	module_dep_simple.o                  \
	module_aer_drydep.o                  \
	module_mixactivate.o                 \
	module_mixactivate_wrappers.o        \
	module_vertmx_wrf.o                  \
	module_mosaic_settling.o             \
	module_cu_gd.o                       \
	module_ctrans_aqchem.o               \
	module_ctrans_kfeta.o                \
	module_ctrans_general.o              \
	module_mosaic_wetscav.o              \
	module_optical_averaging.o

schemes: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)

# DEPENDENCIES:

module_mosaic_addemiss.o:               \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_data_gocart_dust.o           \
	../mpas_atmchem_gocart_dust_init.o  \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_dep_simple.o:                 \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	../mpas_atmchem_constants.o         \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_gocart_drydep.o:                 \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	../mpas_atmchem_constants.o         \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_mixactivate_wrappers.o:          \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_mixactivate.o                \
	module_data_sorgam.o                \
	../mpas_atmchem_constants.o         \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_mixactivate.o:                   \
	module_cal_cldfra.o                 \
	../mpas_atmchem_constants.o         \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_mosaic_settling.o:               \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_data_mosaic_therm.o          \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_ctrans_kfeta.o:                  \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_mosaic_wetscav.o             \
	../mpas_atmchem_utilities.o         \
	../mpas_atmchem_constants.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_ctrans_general.o:                  \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_cu_gd.o                      \
	module_ctrans_aqchem.o              \
	module_mosaic_wetscav.o             \
	../mpas_atmchem_vars.o              \
	../mpas_atmchem_utilities.o         \
	../mpas_atmchem_constants.o         \
	../../physics/mpas_atmphys_vars.o

module_mosaic_wetscav.o:                \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../mpas_atmchem_constants.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_optical_averaging.o:             \
	module_data_binaer_asect.o          \
	module_data_binaer_other.o          \
	module_data_mosaic_therm.o          \
	module_data_rrtmgaeropt.o           \
	../mpas_atmchem_configs.o           \
	../mpas_atmchem_utilities.o         \
	../../physics/mpas_atmphys_vars.o   \
	../mpas_atmchem_vars.o

module_data_sorgam.o:                   \
	module_data_radm2.o

module_data_rrtmgaeropt.o:                   \
	../mpas_atmchem_vars.o
clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I.. -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I.. -I../../../external/esmf_time_f90 -I../../physics -I../../physics/physics_wrf
endif
