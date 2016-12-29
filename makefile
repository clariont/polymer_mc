all: main

CFLAGS= -O3

clean:
	-rm *.o rosenbluth
        
main: genarray.h
#	/opt/intel/bin/icc $(CFLAGS) -o brown brownianmd.cpp -I/home/angelo/ckt2110/libraries/include -L/home/angelo/ckt2110/libraries/lib/ -lm -lgsl -lgslcblas --static
	g++ $(CFLAGS) -o rosenbluth growpolys.cpp -I/home/angelo/ckt2110/libraries/include -L/home/angelo/ckt2110/libraries/lib/ -lm -lgsl -lgslcblas --static
