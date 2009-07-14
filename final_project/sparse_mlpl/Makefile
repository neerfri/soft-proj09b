CFLAGS  =  -Wall -g -ansi -pedantic-errors -I. -I./../shared
CC = gcc
PROJ_NAME = sparse_matrix

all: $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME).o $(PROJ_NAME) main.o

$(PROJ_NAME): $(PROJ_NAME).o main.o
	$(CC) $(CFLAGS) $(PROJ_NAME).o main.o -o $(PROJ_NAME)

$(PROJ_NAME).o:
	$(CC) $(CFLAGS) -c $(PROJ_NAME).c -o $(PROJ_NAME).o
	
shared.o:
	cd ./../shared; make

main.o:
	$(CC) $(CFLAGS) -c main.c -o main.o

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
