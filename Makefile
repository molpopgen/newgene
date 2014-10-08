CC=cc
CXX=c++ 
#OPT=
OPT=-O2
#DEBUG=-g
DEBUG=
DEBUG=-DNDEBUG
PROFILE=
#PROFILE=-pg

CFLAGS= $(OPT) $(DEBUG) $(PROFILE) -Wall -W -ansi -pedantic -I.
CXXFLAGS=$(CFLAGS) -std=c++11
GSL=-lgsl -lgslcblas
SEQ=-lsequence
LIBS=$(GSL) $(SEQ) -lz
TARGETS=dcoal test_new_coal
all: dcoal.o dcoalmk.o cnvcoal.o test_new_coal.o coalesce.o util.o rec.o ca2.o egc.o fixation2.o summstats.o sfs.o jointsfs.o times.o tajima90.o givenT.o fixation2_old.o fixneutral.o transition_sim.o
	$(CXX) $(CXXFLAGS)  -o dcoal dcoal.o coalesce.o util.o rec.o ca2.o egc.o fixation2.o $(LIBS) 
	$(CXX) $(CXXFLAGS)  -o dcoalmk dcoalmk.o coalesce.o util.o rec.o ca2.o egc.o fixation2.o $(LIBS)
	$(CXX) $(CXXFLAGS)  -o cnvcoal cnvcoal.o coalesce.o util.o rec.o ca2.o egc.o fixation2.o $(LIBS)
#	$(CXX) $(CXXFLAGS) $(LIBS) -o fixneutral fixneutral.o
	$(CXX) $(CXXFLAGS)  -o summstats summstats.o $(LIBS)
	$(CXX) $(CXXFLAGS)  -o sfs sfs.o $(LIBS)
	$(CXX) $(CXXFLAGS)  -o jointsfs jointsfs.o $(LIBS)
#	$(CXX) $(CXXFLAGS) -lsequence -o givenT givenT.o tajima90.o
#	$(CXX) $(CXXFLAGS) -lsequence -o times times.o -lgsl -lgslcblas
#	$(CXX) $(CXXFLAGS) $(LIBS) -o test_new_coal test_new_coal.o
#	$(CXX) $(CXXFLAGS) $(GSL) -o transition_sim transition_sim.o	

dcoal.o: ca.hpp rec.hpp egc.hpp util.hpp fixation.hpp constants.hpp
dcoalmk.o: ca.hpp rec.hpp egc.hpp util.hpp fixation.hpp constants.hpp
cnvcoal.o: ca.hpp rec.hpp egc.hpp util.hpp fixation.hpp constants.hpp
coalesce.o: coalesce.hpp
util.o: util.hpp
rec.o: util.hpp
ca.o: coalesce.hpp ca.hpp util.hpp
ca2.o: coalesce.hpp ca.hpp util.hpp
egc.o: egc.hpp coalesce.hpp util.hpp
fixation.o: fixation.hpp rec.hpp ca.hpp egc.hpp
fixation2.o: fixation.hpp rec.hpp ca.hpp egc.hpp
fixation2_old.o: fixation.hpp rec.hpp ca.hpp egc.hpp
tajima90.o: tajima90.hpp

dist:
	mkdir newgene
	cp Makefile *.cc *.hpp newgene
	tar czf newgene.tar.gz newgene
	rm -rf newgene
clean:
	rm -f *.o $(TARGETS)
