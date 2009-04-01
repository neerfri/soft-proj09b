all: check_policy

clean:
	-rm check_policy.o check_policy

check_policy: check_policy.o
	gcc -g -lm -o check_policy check_policy.o

check_policy.o: check_policy.c
	gcc -c -ansi -pedantic-errors -Wall -g check_policy.c 	
