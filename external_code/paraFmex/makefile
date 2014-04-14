# makefile for programs and generator

CPP = g++
#CPPFLAGS = -O4 
CPPFLAGS = -g 
SOURCES = paraFmex.cpp Parser.cpp Graph.cpp Solver.cpp Utils.cpp stdafx.cpp
HEADERS = Parser.h Graph.h Solver.h Utils.h stdafx.h
EXECUTABLES = paraF.exe paraF
LIBS = -lm

all: paraF.exe

paraF.exe: $(SOURCES) $(HEADERS)
	$(CPP) $(CPPFLAGS) $(LIBS) -o paraF.exe $(SOURCES)

clean: 
	rm -f $(EXECUTABLES)
