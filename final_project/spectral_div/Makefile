CFLAGS  =  -Wall -g -ansi -pedantic-errors -I../sparse_mlpl -L../sparse_mlpl -I../shared -L../shared -lm
CC = gcc
PROJ_NAME = spectral_div
LIBS := $(wildcard ./../shared/*.o)

all:
	$(CC) $(CFLAGS) $(PROJ_NAME).c ./../sparse_mlpl/sparse_matrix.o $(LIBS) -o $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME)
