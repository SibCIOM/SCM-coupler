#$Id: Makefile 7 2015-03-31 07:44:54Z plat $
#$HeadURL: $

-include ../Driver/Defs
#F90 = mpif90

#XPATH = ./../../Run
##XPATH = ../../../Run_5/RunNew
#NETCDF_DIR = /usr/local/netcdf-4.1.1
#NETCDF = -L$(NETCDF_DIR)/lib -lnetcdf -lm
DEFS= -DCOMPAQ -DXRES=$(DXRES) -DYRES=$(DYRES) -DM_MODE=$(DM_MODE)
##DEFS=
##DEFS= -Dcoupled -DNPROC_X=1 -DNPROC_Y=1
FFLAGS = -O1 -fpe-all=3 -r8 -w -ip -ipo -132 -Vaxlib -frecord-marker=8 -assume byterecl -mcmodel large -shared-intel $(DEFS) -I$(NETCDF_INCDIR) -I$(CPATH) -cpp
#FFLAGS = -g -inline_debug_info -fpe-all=0 -check all -r8 -w -132 -Vaxlib -frecord-marker=8 -assume byterecl -mcmodel large -shared-intel $(DEFS) -I$(NETCDF_INCDIR) -I$(CPATH) -cpp

SRC = \
atm.F coupler.F flux_ao.F ice_indexes.F atm_indexes.F mpi_tools.F ocn_indexes.F prepare.F \
struct.F ncepdata.F tides.F timer.F tides_data.F read_atmdata.F cpl_indexes.F \
rivers.F conversion.F Module_coare.F msgtypes.F lnd_indexes.F

OBJS = \
atm.o coupler.o flux_ao.o ice_indexes.o atm_indexes.o ocn_indexes.o prepare.o \
struct.o ncepdata.o tides.o ../MPI/mpi_tools.o ../Utils/cdf2D.o \
../Utils/cdftseries.o timer.o tides_data.o read_atmdata.o cpl_indexes.o \
rivers.o conversion.o Module_coare.o msgtypes.o lnd_indexes.o

.F.o:
	$(F90) $(FFLAGS) -c $<

all: cpl

test: test.o srfflux1.o srfflux2.o
	$(F90) $(FFLAGS) -o $@ test.o srfflux1.o srfflux2.o

cpl: $(OBJS)
	$(F90) $(FFLAGS) -o $(XPATH)/$@  $(OBJS)  $(NETCDF) $(MPI_LIB)

clean:
	rm -f $(OBJS) *.mod

../MPI/mpi_tools.o: ../MPI/mpi_tools.F
	$(F90) $(FFLAGS) -c -o../MPI/mpi_tools.o ../MPI/mpi_tools.F

../Utils/cdf2D.o: ../Utils/cdf2D.F
	$(F90) $(FFLAGS) -c -o../Utils/cdf2D.o ../Utils/cdf2D.F

../Utils/cdftseries.o: ../Utils/cdftseries.F
	$(F90) $(FFLAGS) -c -o../Utils/cdftseries.o ../Utils/cdftseries.F

Module_coare.o: Module_coare.F
	$(F90) -free $(FFLAGS) -c $<

atm.o: struct.o

coupler.o: struct.o ice_indexes.o atm_indexes.o ocn_indexes.o lnd_indexes.o timer.o tides_data.o cpl_indexes.o conversion.o msgtypes.o

flux_ao.o: struct.o Module_coare.o

prepare.o: struct.o timer.o tides_data.o conversion.o

get_atmdata.o: struct.o

conversion.o: struct.o
