CC=gcc
CFLAGS=-O3 -std=c99 -fopenmp -march=native
LDFLAGS=-fopenmp -O3

EXE=fidi
OBJS=fd_main.o fd.o paste.o stencils.o

all: $(EXE) 

$(EXE):	$(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXE) *.o
