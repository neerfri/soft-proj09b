CFLAGS  =  -Wall -g -ansi -pedantic-errors -I ./../sparse_mlpl -L ./../sparse_mlpl
CC = gcc
PROJ_NAME = shared

all: $(PROJ_NAME).o

clean:
	rm -f $(PROJ_NAME).o $(PROJ_NAME) main.o

$(PROJ_NAME).o:
	$(CC) $(CFLAGS) -c $(PROJ_NAME).c -o $(PROJ_NAME).o

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
