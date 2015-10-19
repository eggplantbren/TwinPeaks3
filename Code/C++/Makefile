CFLAGS = -O3 -std=c++11 -Wall -Wextra -pedantic
OBJECTS = Atoms.o Potts.o RNG.o ScalarType.o SimpleExample.o SpikeSlab.o Utils.o

default:
	g++ $(CFLAGS) -c *.cpp Models/*.cpp
	g++ -std=c++11 -o main main.o $(OBJECTS)
	rm -f *.o

