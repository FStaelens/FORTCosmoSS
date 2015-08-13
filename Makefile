BASE = $(shell /bin/pwd)
EXED = $(BASE)/$(EXECTBL)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src

vpath %.o $(OBJD)
vpath %.f90 $(SRCD)
vpath %.f $(SRCD)

EXECTBL = main

SRCS = scales.f90 inifile.f90 input.f90 grid.f90 math_lib.f90 fparser.f90 potential.f90 common.f90 hydro.f90 sources.f90 constraints.f90 BVP_LA.f BVP_M.f90 output.f90 initial_profiles.f90 init.f90 evolution.f90 main.f90 

OBJS = scales.o inifile.o input.o grid.o math_lib.o fparser.o potential.o common.o hydro.o sources.o constraints.o BVP_LA.o BVP_M.o output.o initial_profiles.o init.o evolution.o main.o

LIBS = #-llapack -lblas

F90 = gfortran

# Mine
LOAD_FLAGS =
LOAD_LIBS = 

all: $(EXECTBL)

$(EXECTBL): $(OBJS)
	@echo '' 
	@echo '============ ...building the executable ============'
	@echo ''
	cd $(OBJD); $(F90) $(F90FLAGS) -o ../$(EXECTBL) $(OBJS) $(LIBS)
	@echo ''
	@echo '============ ...Done! ============'

%.o $(OBJD)/%.o: %.f
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(F90) -c $(F90FLAGS) $(SRCD)/$(*).f

%.o $(OBJD)/%.o: %.f90
	@echo '' 
	@echo '============ ...building' $*.o '============'
	@echo ''
	cd $(OBJD); $(F90) -c $(F90FLAGS) $(SRCD)/$(*).f90

flush:
	rm -f $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.lst $(OBJD)/*.d $(OBJD)/*.pc*

clean:
	rm -f $(OBJD)/*.o $(OBJD)/*.kmo $(OBJD)/*.mod $(OBJD)/*.d $(OBJD)/*.pc* $(OBJD)/*.lst $(OBJD)/core

.SUFFIXES: $(SUFFIXES) .f90 $(SUFFIXES) .f

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(F90FLAGS) -c $<






