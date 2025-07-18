#####################################################
#
#  Makefile for the whole system:
#  
######################################################
CC = gcc -Wall
CXX = g++ -Wall
CFLAGS = -g -O2
CXXFLAGS= -g -O2
VERSION = 1.0
MODELNAME = CRNLIB$(VERSION)

HEADERS =	src/Epi.h

MAINLIB = 	src/lib/libepi.a

#MAINLINK =	-lepi -lcln -lginac -lglpk
MAINLINK =	-lepi -lcln

OBJDIR = 	src/obj

AR      =       ar crs

OBJECTS =	src/obj/common.o

project: bin/SIRgeneral bin/SIRimmunity bin/SIRcomplete bin/SIRring bin/SIRsfring bin/SIRgrid bin/SIRrandom bin/SIRgamma bin/SIRtorus bin/SIRvax bin/resample bin/SIRoutsize bin/SIRleakeffect

clean: 
	rm bin/* src/lib/libepi.a src/obj/*.o

src/obj/%.o: src/%.c $(HEADERS)
	$(CXX) -c -o $@ $<

src/lib/libepi.a: $(HEADERS) $(OBJECTS)
	$(AR) src/lib/libepi.a $(OBJECTS)



##################

bin/SIRgeneral: src/SIRgeneral.c $(MAINLIB)
	$(CXX) -o bin/SIRgeneral src/SIRgeneral.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRimmunity: src/SIRimmunity.c $(MAINLIB)
	$(CXX) -o bin/SIRimmunity src/SIRimmunity.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRcomplete: src/SIRcomplete.c $(MAINLIB)
	$(CXX) -o bin/SIRcomplete src/SIRcomplete.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRring: src/SIRring.c $(MAINLIB)
	$(CXX) -o bin/SIRring src/SIRring.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRsfring: src/SIRsfring.c $(MAINLIB)
	$(CXX) -o bin/SIRsfring src/SIRsfring.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRgrid: src/SIRgrid.c $(MAINLIB)
	$(CXX) -o bin/SIRgrid src/SIRgrid.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRrandom: src/SIRrandom.c $(MAINLIB)
	$(CXX) -o bin/SIRrandom src/SIRrandom.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRgamma: src/SIRgamma.c $(MAINLIB)
	$(CXX) -o bin/SIRgamma src/SIRgamma.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRtorus: src/SIRtorus.c $(MAINLIB)
	$(CXX) -o bin/SIRtorus src/SIRtorus.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRvax: src/SIRvax.c $(MAINLIB)
	$(CXX) -o bin/SIRvax src/SIRvax.c -lm -Lsrc/lib $(MAINLINK)

bin/resample: src/resample.c $(MAINLIB)
	$(CXX) -o bin/resample src/resample.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRoutsize: src/SIRoutsize.c $(MAINLIB)
	$(CXX) -o bin/SIRoutsize src/SIRoutsize.c -lm -Lsrc/lib $(MAINLINK)

bin/SIRleakeffect: src/SIRleakeffect.c $(MAINLIB)
	$(CXX) -o bin/SIRleakeffect src/SIRleakeffect.c -lm -Lsrc/lib $(MAINLINK)
