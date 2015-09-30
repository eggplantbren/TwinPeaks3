CFLAGS = -std=c++11 -m64 -O3 -flto -march=native -Wall -Wextra -pedantic -DNDEBUG
OBJECTS = RNG.o SimpleExample.o Utils.o

default:
	g++ $(CFLAGS) -c *.cpp Models/*.cpp
	g++ -std=c++11 -o main main.o $(OBJECTS)
	rm -f *.o

