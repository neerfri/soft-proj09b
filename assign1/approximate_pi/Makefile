ifeq ($(shell uname),Darwin)
PEDANTIC = "" 
else
PEDANTIC = "-pedantic-errors"
endif 
all: approximate_pi

clean:
	-rm approximate_pi.o approximate_pi

approximate_pi: approximate_pi.o
	gcc -ansi $(PEDANTIC) -g -lm -o approximate_pi approximate_pi.o

approximate_pi.o: approximate_pi.c
	gcc -ansi $(PEDANTIC) -c -Wall -g approximate_pi.c 	
