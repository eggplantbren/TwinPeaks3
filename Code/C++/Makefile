CFLAGS = -O3 -std=c++11 -Wall -Wextra -pedantic
OBJECTS = Atoms.o Context.o RNG.o ScalarType.o SimpleExample.o Utils.o

default:
	g++ -I . $(CFLAGS) -c *.cpp Models/Atoms.cpp Models/SimpleExample.cpp
	g++ -std=c++11 -o main main.o $(OBJECTS) -lpthread
	rm -f *.o

