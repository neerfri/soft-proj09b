include Makefile.common
SOURCE_FILES := $(wildcard *.c)
SOURCE_NAMES := $(basename $(SOURCE_FILES))
OBJECT_FILES := $(foreach SOURCE_NAME,$(SOURCE_NAMES), $(SOURCE_NAME).o)


all: $(OBJECT_FILES)

clean:
	rm -f $(OBJECT_FILES)

%.o:%.c
	$(CC) $(CFLAGS) -c $*.c -o $*.o
