CC=gcc
CFLAGS=-O3
LIBFLAGS=-lm
OMPFLAGS=-fopenmp

OBJDIR = obj

_OBJ = args.o data.o setup.o vtk.o boundary.o md.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories

all: directories md

obj/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

md: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBFLAGS) $(OMPFLAGS)

clean:
	rm -Rf $(OBJDIR)
	rm -f md

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
