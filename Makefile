##
## Makefile for FDPML project
##

DEBUG 		= Yes

#
# List of Fortran modules that should be generated from source files
# (with .f90 suffix) in this directory
#
MODULES		= kinds.mod \
			constants.mod \
			mp_module.mod \
			essentials.mod \
			io_module.mod 

#
# Sources that go into the executable (not modules):
#
SOURCES		= genenerate_domain.f90

#
# Generate list of object files for the SOURCES listed
# above:
#
OBJECTS		= kinds.o \
			constants.o \
			mp_module.o \
			essentials.o \
			io_module.o \
			generate_domain.o
#
# What should the executable produced be named?
#
TARGET		= gendomain.out

##
####
####
##

#
# Fortran Compiler (FC) to be used:
#
FC		= mpifort

#
# Fortran FLAGS (FFLAGS) to be used (not preprocessor definitions,
# just flags that affect compilation dialect, optimizations, etc.):
#
# compiler flags for final executable
FFLAGS_FOR_RUN		= -mkl 

# compiler flags for execulatable with call-back and traceback capabilities 
#(add other flags that you might require for debugging)
FFLAGS_FOR_DEBUG	= -mkl -CB -traceback

#
# Fortran Pre-Processor FLAGS (FPPFLAGS) that affect compilation
# (e.g. "-I../common" or "-DPI=3.14159"):
#
FPPFLAGS	= 

#
# LD FLAGS (LDFLAGS) are flags used during the linking stage that
# are specific to the linker:
#
LDFLAGS		= 

#
# Any libraries needed by the linker (e.g. "-lm"):
#
LIBS		= -lmkl_lapack95_ilp64 -lmkl_blas95_ilp64

##
####
####
##

#
# The default target if nothing is supplied with the "make"
# command:
#
default: $(TARGET)

#
# Run "make all" command to generate target file
#
all : $(TARGET)

#
# The "clean" target just removes all files (modules, object code)
# produced during the build process:
#
.PHONY: clean
clean::
	$(RM) -rf $(TARGET) $(MODULES) $(OBJECTS)

#
# Rule to produce the target executable:
#
$(TARGET): $(MODULES) $(OBJECTS)
ifeq ($(DEBUG), Yes)
	$(FC) $(FFLAGS_FOR_DEBUG) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)
else
	$(FC) $(FFLAGS_FOR_RUN) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)
endif
#
# Pattern rule that produces object code from source code:
#
%.o: %.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -o $@ -c $<
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -o $@ -c $<
endif

#
# Pattern rule that produces a Fortran module from source
# code:
#
%.mod: %.f90
ifeq ($(DEBUG), Yes)
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_DEBUG) -c $< 
else
	$(FC) $(FPPFLAGS) $(FFLAGS_FOR_RUN) -c $< 
endif


