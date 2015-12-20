# Makefile copied from homework 4 and edited as required
HOME = .
VPATH= . $(HOME)/RectArray $(HOME)/fftTools
DIM=3
system := $(shell uname)
CXX=g++
#CXX=clang++
FFTWDIR=/usr/local/opt/fftw
FFTINCFLAGS = -I$(FFTWDIR)/include
CFLAGS = -g -Wall
#CFLAGS = -O3
CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d

CFLAGS += -I. -I$(HOME)/RectArray -I$(HOME)/fftTools -DDIM=$(DIM) $(FFTINCFLAGS) -D DIM=$(DIM)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 -L./lib -lfft$(DIM)D -lfft1D
OBJS:=atom.o molecule.o system.o parameters.o forcefield.o moves.o Nonbonded.o main.o

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

lib/libfft$(DIM)D.a:$(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
	cd $(HOME)/fftTools;rm *.o;make libfft$(DIM)D.a CXX=$(CXX) DIM=$(DIM) FFTWDIR=$(FFTWDIR) FFTINCFLAGS=$(FFTINCFLAGS)

lib/libfft1D.a:$(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
	cd $(HOME)/fftTools;rm *.o;make libfft1D.a CXX=$(CXX) DIM=$(DIM) FFTWDIR=$(FFTWDIR) FFTINCFLAGS=$(FFTINCFLAGS)

MD: $(OBJS) GNUmakefile ./lib/libfft$(DIM)D.a ./lib/libfft1D.a
	$(CXX) $(CFLAGS) $(LIBS) -o MD.exe $(OBJS)

clean:
	rm  *.o *.exe *.d

realclean:
	rm  *.exe *.d *.o $(OBJS);cd ../fftTools;make clean

-include $(OBJS:.o=.d)
