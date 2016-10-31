CC = gcc
# CFLAGS = -Wall -std=c99 -pedantic -march=k8 -O2
CFLAGS = -Wall -std=c99 -pedantic -O2
DEBUGFLAGS=-g -lefence
LIBS = -lgsl -lgslcblas -lm

DEPS = header.h
OBJ = 3layer.o model.o convo.o extras.o io.o rti-theory.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

3layer: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) 

clean:
	rm -f 3layer *.o core

