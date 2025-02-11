CC=gcc
NVCC=nvcc
CFLAGS=-O3
LIBFLAGS=-lm

OBJDIR = obj

_OBJ = args.o data.o setup.o vtk.o boundary.o md.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories

all: directories md

obj/%.o: %.cu
	$(NVCC) -c -o $@ $< $(CFLAGS)

obj/md.o: md.cu
	$(NVCC) -arch=sm_86 -c -o $@ $< $(CFLAGS) 

md: $(OBJ)
	$(NVCC) -o $@ $^ $(CFLAGS) $(LIBFLAGS) 

clean:
	rm -Rf $(OBJDIR)
	rm -f md

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
