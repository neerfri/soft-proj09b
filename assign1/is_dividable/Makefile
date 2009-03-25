all: is_dividable

clean:
	-rm is_dividable.o is_dividable

is_dividable: is_dividable.o
	gcc -ansi -pedantic-errors -g -lm -o is_dividable is_dividable.o

is_dividable.o: is_dividable.c
	gcc -ansi -pedantic-errors -c -Wall -g is_dividable.c 	
