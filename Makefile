# By Jixie Zhang
# Makefile for fortran program, you can specify your executable by command line: 
# "make PROG=your_exe"
# Will look for program file itself under current dir, make sure remove
# all unrelavant source file out of current dir
#
MYOS := $(shell uname -s)
ARCH := $(shell arch)
KEY  := $(ARCH)
ifeq ($(ARCH),x86_64)
  ifeq ($(COMPILER_BIT),32)
    KEY := $(ARCH)_m32
  endif
endif
#####################################################
SRCDIR := src
INCDIR := src
OBJDIR := obj.$(KEY)
#####################################################
FF = f77
LD = f77
FFLAGS = -O2 

#fortran options are crutial, some of them should not be included here
#depends on what options have been used in the cernlib 
ifeq ($(ARCH),i686)
 ifeq ($(COMPILER_BIT),32)
    FFLAGS += -m32
  else
    FF = gfortran
    LD = gfortran
#  FFLAGS  += -fno-automatic -finit-local-zero \
 -fbounds-check -funroll-all-loops -fdollar-ok \
 -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wno-globals -Wuninitialized 
 endif
endif

ifeq ($(ARCH),x86_64)
  ifeq ($(COMPILER_BIT),32)
    FFLAGS += -m32
  else
    FF = gfortran
    LD = gfortran
#    FFLAGS += -fno-automatic -fbounds-check \
 -funroll-all-loops -fdollar-ok \
 -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wuninitialized 
#    FFLAGS += -m64
  endif
endif

FFLAGS += -I$(INCDIR)
ifdef DEBUG
 FFLAGS += -g	
endif

################################ Libraries
ifndef SYSLIBS
 SYSLIBS  = -lnsl -lX11
endif

#for SunOS, cernlib 2001 need ibsocket.a
ifeq ($(MYOS),SunOS)
 SYSLIBS  += -lsocket
endif

ifndef CERNLIBS
 #this program used 2001 cernlib /apps/cernlib/sun4_solaris8/2001 in jlabs2,
 ifndef CERN_ROOT
   CERN_ROOT :=/mnt/site/cernlib/i386_rhel6/2005
   ifeq ($(MYOS),SunOS)
     CERN_ROOT := /apps/cernlib/2006/BUILD
   endif
   ifeq ($(ARCH),x86_64)
     ifeq ($(COMPILER_BIT),64)
       CERN_ROOT := /apps/cernlib/x86_64_rhel5/2005
     endif
   endif
 endif
 CERNLIBS = -L$(CERN_ROOT)/lib -lpdflib804 -lmathlib -lphtools \
 -lgeant321 -lpawlib -lgraflib -lgrafX11 -lpacklib -lkernlib 
endif

#if you are a clas member ......
ifdef USE_CLAS_LIB
 ifndef CLASLIBS
  CLASLIBS = -L$(CLAS_LIB) -lbankdefs -lc_bos_io -lfputil  -lrecutl \
  -lbos -lfpack -lcc -lmapmanager -lclasutil
 endif
else
 CLASLIBS =  # use this # to set a space  
endif

LIBS = $(SYSLIBS) $(CLASLIBS) $(CERNLIBS)
#===========================

###################################
ifndef PROG
  PROG = $(shell grep -i program $(SRCDIR)/*.[fF] | awk '$$2=="program"&&$$4==""{print $$3}' | sed s/\ //g)
endif
ifeq ($(PROG),) 
  PROG = $(shell grep -i program $(SRCDIR)/*.[fF] | awk '$$1=="program"&&$$3==""{print $$2}' | sed s/\ //g)  
endif


PROG := $(shell echo ${PROG} | sed s/\ //g)
EXE   = $(PROG).exe
ifdef DEBUG
EXE  :=$(EXE)_debug
endif

SOURCES  = $(wildcard $(SRCDIR)/*.[fF])
#OBJS     = $(shell echo $(SOURCES) | sed s/f\ /o\ /g)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))
OBJS    := $(patsubst $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))

###################################
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FF) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o:$(SRCDIR)/%.F
	$(FF) $(FFLAGS) -c $< -o $@

%.o:%.f
	$(FF) $(FFLAGS) -c $< -o $@

%.o:%.F
	$(FF) $(FFLAGS) -c $< -o $@

all:  $(OBJDIR) $(EXE)

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

$(EXE): $(OBJDIR) $(OBJS) 
	$(FF) -o $(EXE) $(OBJS) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
exe: lib
	$(FF) -o $(EXE) $(OBJS) -L$(OBJDIR) -l$(PROG) $(LIBS) $(FFLAGS)
	@echo "Linking $(EXE) ... done!"

#this part is to support CLAS_PACK make mechanism
lib: $(OBJDIR) $(OBJS)
	@ar -rf $(OBJDIR)/lib${PROG}.a ${OBJS}

clean:
	rm -f $(OBJS) $(OBJDIR)/lib${PROG}.a core

delete:
	rm -f $(EXE)

help:
	@echo PROG = ${PROG}
	@echo SOURCES = ${SOURCES}
	@echo OBJS = ${OBJS}
	@echo EXE = ${EXE}
	@echo LIBS = ${LIBS}
	@echo FFLAGS = ${FFLAGS}
	@echo CLASLIBS = ${CLASLIBS}
