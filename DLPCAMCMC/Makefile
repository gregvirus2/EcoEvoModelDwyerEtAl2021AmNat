ravdebug:
	gcc -I/sw/gsl/include -O3 -c -g DLPCAMCMC.c 
	gcc -L/sw/gsl/lib DLPCAMCMC.o -lgsl -lgslcblas -lm -o DLPCAMCMC


clean:
	rm -f *.o
	rm -f *.exe
	rm -f FixMCInt
