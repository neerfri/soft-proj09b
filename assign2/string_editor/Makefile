ifeq ($(shell uname),Darwin)
PEDANTIC = "" 
else
PEDANTIC = "-pedantic-errors"
endif 
CFLAGS  =  -Wall -g -ansi $(PEDANTIC) -I./../sparse_mlpl -L./../sparse_mlpl -I. -I./../shared -lm -L../shared
CC = gcc

all: string_editor

clean:
	rm -f my_list.o string_editor.o string_editor

string_editor: my_list.o string_editor.o
	$(CC) $(CFLAGS) my_list.o string_editor.o -o string_editor

my_list.o string_editor.o: my_list.h

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
