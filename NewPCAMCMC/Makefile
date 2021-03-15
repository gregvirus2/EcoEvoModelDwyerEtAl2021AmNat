ravdebug:
	gcc -I/sw/gsl/include -O3 -c -g NewPCAMCMC.c 
	gcc -L/sw/gsl/lib NewPCAMCMC.o -lgsl -lgslcblas -lm -o NewPCAMCMC


clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt
