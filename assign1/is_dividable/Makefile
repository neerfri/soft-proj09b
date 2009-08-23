ifeq ($(shell uname),Darwin)
PEDANTIC = "" 
else
PEDANTIC = "-pedantic-errors"
endif 
all: is_dividable

clean:
	-rm is_dividable.o is_dividable

is_dividable: is_dividable.o
	gcc -ansi $(PEDANTIC) -g -lm -o is_dividable is_dividable.o

is_dividable.o: is_dividable.c
	gcc -ansi $(PEDANTIC) -c -Wall -g is_dividable.c 	
