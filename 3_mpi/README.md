# A Simple Molecular Dynamics Application Using Cell Lists

This is a simple MD simulation application, implemented using a cell list approach with particles stored in linked lists.

By default, the problem is set up with a 50 x 50 domain, and a cut-off distance of 2.5.

At the end of execution, two VTK files are produced for visualisation purposes. They can be loaded into VisIt for analysis.

## Building

The application can be built with the provided Makefile. e.g.

```
$ make
```

This will build an `md` binary.

## Running

The application can be run in its default configuration with:

```
$ ./md
```

This will output its status every 100 iterations. At the end of execution, two VTK files will be produced (The mesh in a .vti file, the particles in a .vtp file). 

There are numerous other options available. These can be queried with:

```
$ ./md --help
```

To write out the simulation state every 100 iterations, you could use:

```
$ mkdir out
$ ./md -c -o out/my_sim
```
