CFLAGS = -m64 -Ofast -flto -march=native -funroll-loops -Wall -Wextra -ansi -pedantic -DNDEBUG
LIBS =  -ldnest3 -lgsl -lgslcblas -lboost_system -lboost_thread

default:
	g++ $(CFLAGS) -c *.cpp
	g++ -o main *.o $(LIBS)
	rm -f *.o

