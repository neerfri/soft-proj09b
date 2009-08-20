CFLAGS  =  -Wall -g -ansi -pedantic-errors -I../sparse_mlpl -L../sparse_mlpl -I../shared -L../shared -lm
CC = gcc
PROJ_NAME = modularity_mat
LIB_SOURCES := $(wildcard ./../shared/*.c)
LIBS = $(foreach LIB_NAME,$(basename $(LIB_SOURCES)), $(LIB_NAME).o)

all: LIBS sparse_matrix.o
	$(CC) $(CFLAGS) $(PROJ_NAME).c ./../sparse_mlpl/sparse_matrix.o $(LIBS) -o $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME)

LIBS:
	cd ./../shared; make
	
sparse_matrix.o:
	cd ./../sparse_mlpl; make sparse_matrix.o