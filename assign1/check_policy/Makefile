ifeq ($(shell uname),Darwin)
PEDANTIC = "" 
else
PEDANTIC = "-pedantic-errors"
endif 
all: check_policy

clean:
	-rm check_policy.o check_policy

check_policy: check_policy.o
	gcc -g -lm -o check_policy check_policy.o

check_policy.o: check_policy.c
	gcc -c -ansi $(PEDANTIC) -Wall -g check_policy.c 	
