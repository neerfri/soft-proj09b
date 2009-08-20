CFLAGS  =  -Wall -g -ansi -pedantic-errors -I. -I./../shared -lm -L../shared
CC = gcc
PROJ_NAME = sparse_matrix
LIB_SOURCES := $(wildcard ./../shared/*.c)
LIBS = $(foreach LIB_NAME,$(basename $(LIB_SOURCES)), $(LIB_NAME).o)

all: $(PROJ_NAME).o
	$(CC) $(CFLAGS) main.c $(PROJ_NAME).o $(LIBS) -o $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME).o $(PROJ_NAME)

#$(PROJ_NAME).o: LIBS 
#	$(CC) $(CFLAGS) -c $(PROJ_NAME).c -o $(PROJ_NAME).o
	
LIBS:
	cd ./../shared; make
	
%.o:%.c LIBS
	$(CC) $(CFLAGS) -c $*.c -o $*.o