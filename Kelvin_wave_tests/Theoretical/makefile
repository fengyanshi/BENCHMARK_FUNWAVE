#!/bin/make

SRC = $(NAME).c
BINARY = $(NAME)
MEX_SUFFIX = mexa64
 
MATLAB = /pkg/suse11/matlab/R2013a
INTEL = /pkg/suse11/intel/2013/composer_xe_2013.3.163
MKLROOT = $(INTEL)/mkl
 
CC = icpc
OPTS = -DMATLAB_MEX_FILE -openmp -DMKL_ILP64 -fPIC -fno-omit-frame-pointer -DMX_COMPAT_32 -O3 -xHost

OBJS = $(SRC:.c=.o)
 
INC = -I$(MATLAB)/extern/include -I$(MATLAB)/simulink/include -I$(MKLROOT)/include
LIBS = -shared -Wl,--version-script,$(MATLAB)/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined 
LIBS += -Wl,-rpath-link,$(MATLAB)/bin/glnxa64 -L$(MATLAB)/bin/glnxa64 -lmx -lmex -lmat 
LIBS += -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm 

.SUFFIXES: .c .o
 
.c.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@
 
Go: $(OBJS)
	$(CC) $(OBJS) $(OPTS) -o $(BINARY).$(MEX_SUFFIX) $(LIBS)
	@echo Binary created!!
 
clean:
	set nonomatch; rm -f $(BINARY).$(MEX_SUFFIX) *.o