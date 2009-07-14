CFLAGS  =  -Wall -g -ansi -pedantic-errors -I../sparse_mlpl -L../sparse_mlpl -I../shared -L../shared
CC = gcc
PROJ_NAME = eigen_pair

all: $(PROJ_NAME)

clean:
	rm -f $(PROJ_NAME).o $(PROJ_NAME)

$(PROJ_NAME): $(PROJ_NAME).o shared.o
	$(CC) $(CFLAGS) $(PROJ_NAME).o ./../sparse_mlpl/sparse_matrix.o ./../shared/shared.o -o $(PROJ_NAME)

$(PROJ_NAME).o:
	$(CC) $(CFLAGS) -c $(PROJ_NAME).c -o $(PROJ_NAME).o
	
shared.o:
	cd ./../shared; make;

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
