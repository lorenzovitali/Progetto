CXX = g++ -g
CPPFLAGS = -std=c++11 -Wall
DIR = $(mkEigenInc)
LIB = $(mkEigenLib)

EXE = main
OBJECTS = main.o Option.o PDE.o FDM.o SOR.o

all:$(EXE)

$(EXE):$(OBJECTS)
	$(CXX) $(CPPFLAGS) -I$(DIR) $^ -o $@

main.o: main.cc
	$(CXX) -c $(CPPFLAGS) -I$(DIR) $<

FDM.o: FDM.cc FDM.hh
	$(CXX) -c $(CPPFLAGS) -I$(DIR) $<

SOR.o: SOR.hh PDE.hh

PDE.o: PDE.hh Option.hh

Option.o: Option.hh

clean:
	$(RM) *.o
