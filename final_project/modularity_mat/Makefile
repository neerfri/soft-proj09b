CFLAGS  =  -Wall -g -ansi -pedantic-errors -I../sparse_mlpl -L../sparse_mlpl
CC = gcc
PROJ_NAME = modularity_mat

all: $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME).o $(PROJ_NAME)

$(PROJ_NAME): $(PROJ_NAME).o sparse_matrix.o
	$(CC) $(CFLAGS) $(PROJ_NAME).o ./../sparse_mlpl/sparse_matrix.o -o $(PROJ_NAME)

$(PROJ_NAME).o:
	$(CC) $(CFLAGS) -c $(PROJ_NAME).c -o $(PROJ_NAME).o

sparse_matrix.o:
	$(CC) $(CFLAGS) -c ./../sparse_mlpl/sparse_matrix.c -o  ./../sparse_mlpl/sparse_matrix.o

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
