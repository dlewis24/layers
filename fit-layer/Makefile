CC = gcc
# CFLAGS = -Wall -std=c99 -pedantic -march=k8 -O2
CFLAGS = -Wall -std=c99 -pedantic -O2
DEBUGFLAGS=-g -lefence
LIBS = -lgsl -lgslcblas -lm

DEPS = header.h
OBJ = fit-layer.o model.o convo.o extras.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

fit-layer: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) 

clean:
	rm -f fit-layer *.o core

