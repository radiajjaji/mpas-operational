.SUFFIXES: .F .o

dummy:
	echo "******* make chemistry core *******"

OBJS = \
	mpas_atmchem_todynamics.o	\
	mpas_atmchem_driver.o	        \
	mpas_atmchem_constants.o	\
	mpas_atmchem_driver_emission.o	\
	mpas_atmchem_driver_drydep.o	\
	mpas_atmchem_driver_convection.o	\
	mpas_atmchem_driver_wetscav.o	\
	mpas_atmchem_driver_optical.o	\
	mpas_atmchem_gocart_dust_init.o \
	mpas_atmchem_interface.o	\
	mpas_atmchem_configs.o	        \
	mpas_atmchem_utilities.o	\
	mpas_atmchem_packages.o	        \
	mpas_atmchem_control.o	        \
	mpas_atmchem_binaer_init.o      \
	mpas_atmchem_vars.o

all: core_chemistry core_schemes

core_schemes:
	(cd schemes; make all)

core_chemistry: $(OBJS)
	ar -ru libchem.a $(OBJS)

# DEPENDENCIES:
mpas_atmchem_todynamics.o: \
	mpas_atmchem_vars.o

mpas_atmchem_vars.o: \
	mpas_atmchem_utilities.o	\
	mpas_atmchem_configs.o		\
	../physics/mpas_atmphys_vars.o

mpas_atmchem_binaer_init.o: \
	mpas_atmchem_utilities.o	\
	mpas_atmchem_constants.o	\
	./schemes/module_data_binaer_asect.o \
	./schemes/module_data_binaer_other.o \
	./schemes/module_data_mosaic_therm.o \
	./schemes/module_mosaic_wetscav.o    \
	mpas_atmchem_vars.o

mpas_atmchem_driver.o: \
	mpas_atmchem_binaer_init.o      \
	mpas_atmchem_driver_optical.o	\
	mpas_atmchem_driver_emission.o  \
	mpas_atmchem_driver_drydep.o	\
	mpas_atmchem_driver_convection.o	\
	mpas_atmchem_driver_wetscav.o	\
	mpas_atmchem_interface.o	    \
	mpas_atmchem_utilities.o	    \
	mpas_atmchem_configs.o	        \
	../physics/mpas_atmphys_vars.o  \
	mpas_atmchem_vars.o

mpas_atmchem_interface.o: \
	mpas_atmchem_utilities.o	  \
	../physics/mpas_atmphys_vars.o\
	mpas_atmchem_vars.o

mpas_atmchem_gocart_dust_init.o:  \
	mpas_atmchem_utilities.o	  \
	./schemes/module_data_gocart_dust.o

mpas_atmchem_driver_emission.o: \
	mpas_atmchem_vars.o	   \
	../physics/mpas_atmphys_vars.o\
	mpas_atmchem_utilities.o	  \
	mpas_atmchem_configs.o	\
	mpas_atmchem_interface.o	  \
	mpas_atmchem_gocart_dust_init.o	\
	./schemes/module_mosaic_addemiss.o

mpas_atmchem_driver_drydep.o: \
	mpas_atmchem_utilities.o	  \
	mpas_atmchem_vars.o	   \
	../physics/mpas_atmphys_vars.o\
	mpas_atmchem_configs.o	\
	mpas_atmchem_interface.o	  \
	./schemes/module_dep_simple.o \
	./schemes/module_vertmx_wrf.o \
	./schemes/module_aer_drydep.o \
	./schemes/module_data_sorgam.o \
	./schemes/module_cal_cldfra.o \
	./schemes/module_mixactivate.o \
	./schemes/module_mixactivate_wrappers.o \
	./schemes/module_mosaic_settling.o \
	./schemes/module_gocart_drydep.o

mpas_atmchem_driver_convection.o:       \
	mpas_atmchem_vars.o	        \
	../physics/mpas_atmphys_vars.o  \
	mpas_atmchem_utilities.o        \
	mpas_atmchem_configs.o          \
	mpas_atmchem_interface.o	\
	./schemes/module_ctrans_aqchem.o\
	./schemes/module_cu_gd.o	\
	./schemes/module_ctrans_general.o	\
	./schemes/module_ctrans_kfeta.o

mpas_atmchem_driver_wetscav.o:    \
	mpas_atmchem_vars.o	          \
	../physics/mpas_atmphys_vars.o\
	mpas_atmchem_utilities.o	  \
	mpas_atmchem_configs.o        \
	mpas_atmchem_interface.o	  \
	./schemes/module_mosaic_wetscav.o

mpas_atmchem_driver_optical.o:    \
	mpas_atmchem_vars.o	          \
	../physics/mpas_atmphys_vars.o\
	mpas_atmchem_utilities.o	  \
	mpas_atmchem_configs.o        \
	mpas_atmchem_interface.o	  \
	./schemes/module_data_rrtmgaeropt.o \
	./schemes/module_optical_averaging.o

clean:
	$(RM) *.o *.mod *.f90 libchem.a
	( cd schemes; make clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS)  $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators -I./schemes -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics -I../physics/physics_wrf -I../../framework -I../../operators -I./schemes -I../../external/esmf_time_f90 -L../physics -lphys
endif
