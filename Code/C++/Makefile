CFLAGS = -O3 -std=c++11 -Wall -Wextra -pedantic
OBJECTS = Atoms.o Cambridge.o CambridgeLJ.o Gravity.o Potts.o RNG.o ScalarType.o SimpleExample.o SpikeSlab.o Utils.o

default:
	g++ -I . $(CFLAGS) -c *.cpp Models/*.cpp
	gfortran -O3 -c Models/Cambridge.f95
	g++ -std=c++11 -o main main.o $(OBJECTS) -lpthread -lgfortran
	rm -f *.o

