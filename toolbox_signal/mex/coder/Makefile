# $Id: Makefile,v 1.4 1996/05/30 20:43:29 jmd Exp $
# $Log: Makefile,v $
# Revision 1.4  1996/05/30 20:43:29  jmd
# Added Integrated and TarFile target
#

CXX=g++
CC=gcc
CXXFLAGS= -Wall -O6
#CXXFLAGS= -Wall -g -gstabs+ -ggdb

CCOBJ = \
	Arith.o \
	BitIO.o \
	iHisto.o \
	IntCoding.o \
	allocator.o \
	coder.o \
	coeffset.o \
	entropy.o \
	filter.o \
	global.o \
	image.o \
	quantizer.o \
	transform.o \
	wavelet.o 

MAINOBJ= \
	encode.o \
	decode.o 

OBJ= $(CCOBJ)

INC = \
	Arith.h \
	BitIO.h \
	iHisto.h \
	IntCoding.hh \
	allocator.hh \
	coder.hh \
	coeffset.hh \
	entropy.hh \
	global.hh \
	image.hh \
	metric.hh \
	quantizer.hh \
	transform.hh \
	wavelet.hh

SRC = $(CCOBJ:.o=.cc) $(MAINOBJ:.o=.cc)

#LIBS=  dbgc.o gc.a -lm
LIBS= -lm

DEP= $(OBJ:.o=.d) $(MAINOBJ:.o=.d)

TARGETS = encode decode compare raw2pgm pgm2raw

all:  $(INC) $(SRC) $(DEP) $(TARGETS)

encode: encode.o $(OBJ)
	$(CXX) $(CXXFLAGS)  -o encode encode.o $(OBJ) $(LIBS)

decode: decode.o $(OBJ)
	$(CXX) $(CXXFLAGS)  -o decode decode.o $(OBJ) $(LIBS)

compare: compare.o $(OBJ)
	$(CXX) $(CXXFLAGS)  -o compare compare.o $(OBJ) $(LIBS)

raw2pgm: raw2pgm.o $(OBJ)
	$(CXX) $(CXXFLAGS)  -o raw2pgm raw2pgm.o $(OBJ) $(LIBS)

pgm2raw: pgm2raw.o $(OBJ)
	$(CXX) $(CXXFLAGS)  -o pgm2raw pgm2raw.o $(OBJ) $(LIBS)


TarFile: $(SRC) $(INC) Makefile
	tar -cf - $(SRC) $(INC) Makefile images | gzip > fliit.tar.gz

clean:
	-rm $(OBJ) $(TARGETS) $(DEP) *~ #* core
	-co -u -q $(SRC) $(INC)

TAGS: ${INC} ${SRC}
	etags ${INC} ${SRC}


%.d: %.cc
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $< | sed '\''s/$*.o/& $@/g'\'' > $@'

%.d: %.c
	$(SHELL) -ec '$(CXX) -MM $(CPPFLAGS) $< | sed '\''s/$*.o/& $@/g'\'' > $@'


%.o: %.c
	$(CXX) -c $(CXXFLAGS) $< -o $@

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@


include $(DEP)
