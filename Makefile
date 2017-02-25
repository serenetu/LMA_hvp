# $Id: Makefile.users.in,v 1.6 2011/03/24 16:20:52 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /home/a_c/ctu/cps/build/..
BUILDDIR = /home/a_c/ctu/cps/build
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = /usr/local/openmpi/bin/mpicc
CXX = /usr/local/openmpi/bin/mpicxx
AS  = as
LDFLAGS =  -L/usr/local/scidac/install/qmp-2.4.2-openmpi/lib -lqmp -L/usr/local/scidac/install/qio/lib -lqio -llime -L/usr/local/fftw/double-openmpi1.8/lib -lfftw3  -L/lqcdproj/rbcdwf/chulwoo/CPS/cps_pp/cps-pi0/local/lib -llapack -lblas -llapacke -lgfortran  -L/usr/local/intel/lib/lib64 -lirc

me = $(notdir $(PWD))
BIN = NOARCH.x

#VPATH :=$(SRCDIR)/tests/$(me)



#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -I/usr/local/scidac/install/qmp-mvapich/include -I/usr/local/scidac/install/qio/include -I/usr/local/fftw/double-mvapich/include

CFLAGS= -g -O2 -O2 -msse -msse2 -msse3
CXXFLAGS= -g -O2 -msse -msse2 -msse3
ASFLAGS= 
DFLAGS +=  -DUSE_SSE -DUSE_QMP -DUSE_QIO -DUSE_BLAS

#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#

.PHONY: cps clean


LIBLIST =\
  $(BUILDDIR)/cps.a \
  /usr/lib64/libgslcblas.a \

#
#  targets
#


all:$(BIN)

.SUFFIXES:
.SUFFIXES:  .o .C .S .c

CSRC :=$(wildcard *.c)
CCSRC :=$(wildcard *.C)
SSRC :=$(wildcard *.S)

COBJ=$(CSRC:.c=.o)
CCOBJ=$(CCSRC:.C=.o)
SOBJ=$(SSRC:.S=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))

$(BIN):  $(OBJS) $(LIBLIST)
	@echo OBJS = $(OBJS)
	$(CXX) $(OBJS) $(LIBLIST) $(LDFLAGS) -o $(BIN)

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<
.S.o:
	$(AS) -o $@ $(ASFLAGS) -c $(INCLIST) $<

cps:
	$(MAKE) -C $(BUILDDIR)

clean:
	rm -f *.dat *.o  $(BIN)
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog
