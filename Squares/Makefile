CFLAGS = -m64 -Ofast -flto -march=native -funroll-loops -Wall -Wextra -ansi -pedantic
LIBS =  -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_system -lboost_thread

default:
	g++ $(CFLAGS) -c *.cpp
	g++ -o main *.o $(LIBS)
	rm -f *.o

