CFLAGS = -O3 -std=c++11 -Wall -Wextra -pedantic
OBJECTS = Atoms.o Gravity.o ImageEntropy.o Potts.o PSF.o RNG.o ScalarType.o SimpleExample.o SpikeSlab.o Utils.o

default:
	g++ -I . $(CFLAGS) -c *.cpp Models/*.cpp
	g++ -std=c++11 -o main main.o $(OBJECTS) -lpthread -lgfortran
	rm -f *.o

