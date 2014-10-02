CFLAGS = -m64 -Ofast -flto -march=native -funroll-loops -Wall -Wextra -ansi -pedantic -DNDEBUG
LIBS =  -ldnest3 -lgsl -lgslcblas -lboost_system -lboost_thread
OBJECTS = Gravity.o Model.o Potts.o SimpleExample.o

default:
	g++ $(CFLAGS) -c *.cpp Models/*.cpp
	g++ -o main main.o $(OBJECTS) $(LIBS)
	g++ $(CFLAGS) -c sorter.cpp
	g++ -o sorter sorter.o $(OBJECTS) $(LIBS)
	rm -f *.o

