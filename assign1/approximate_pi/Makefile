all: approximate_pi

clean:
	-rm approximate_pi.o approximate_pi

approximate_pi: approximate_pi.o
	gcc -ansi -pedantic-errors -g -lm -o approximate_pi approximate_pi.o

approximate_pi.o: approximate_pi.c
	gcc -ansi -pedantic-errors -c -Wall -g approximate_pi.c 	
