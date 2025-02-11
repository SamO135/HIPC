#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "args.h"
#include "data.h"
#include "vtk.h"

int verbose = 0;
int no_output = 0;
int output_freq = 100;
int enable_checkpoints = 0;

static struct option long_options[] = {
	{"cellx",         required_argument, 0, 'x'},
	{"celly",         required_argument, 0, 'y'},
	{"parts-per-dim", required_argument, 0, 'p'},
	{"cellsize",      required_argument, 0, 's'},
	{"cutoff",        required_argument, 0, 'r'},
	{"endtime",       required_argument, 0, 't'},
	{"iters",         required_argument, 0, 'i'},
	{"del-t",         required_argument, 0, 'd'},
	{"freq",          required_argument, 0, 'f'},
	{"seed",          required_argument, 0, 'e'},
	{"noio",          no_argument,       0, 'n'},
	{"output",        required_argument, 0, 'o'},
	{"checkpoint",    no_argument,       0, 'c'},	
    {"verbose",       no_argument,       0, 'v'},
    {"help",          no_argument,       0, 'h'},
	{0, 0, 0, 0}
};
#define GETOPTS "x:y:p:s:r:t:i:d:f:e:no:cvh"

/**
 * @brief Print a help message
 * 
 * @param progname The name of the current application
 */
void print_help(char *progname) {
	fprintf(stderr, "A simple molecular dynamics simulation using the Lennard-Jones potential and cell lists.\n\n");
	fprintf(stderr, "Usage: %s [options]\n", progname);
	fprintf(stderr, "Options and arguments:\n");
	fprintf(stderr, "  -x N, --cellx=N         Cells in X-dimension\n");
	fprintf(stderr, "  -y N, --celly=N         Cells in Y-dimension\n");
	fprintf(stderr, "  -p N, --parts-per-dim=N Set the number of particles per cell, per dimension\n");
	fprintf(stderr, "  -s N, --cellsize=N      Size of each cell in each dimension\n");
	fprintf(stderr, "  -r N, --cutoff=N        Set the cut off size (must be smaller than cell size)\n");
	fprintf(stderr, "  -t N, --endtime=N       Set the end time\n");
	fprintf(stderr, "  -i, --iters=N           Set the number of iterations\n");
    fprintf(stderr, "  -d, --del-t=DELT        Set the simulation timestep size\n");
	fprintf(stderr, "  -f N, --freq=N          Output frequency (i.e. steps between output)\n");
	fprintf(stderr, "  -e N, --seed=N          Set the seed for the random number generator\n");
	fprintf(stderr, "  -n, --noio              Disable file I/O\n");
	fprintf(stderr, "  -o FILE, --output=FILE  Set base filename for particle output (final output will be in BASENAME.vtp)\n");
	fprintf(stderr, "  -c, --checkpoint        Enable checkpointing, checkpoints will be in BASENAME-ITERATION.vtp\n");
	fprintf(stderr, "  -v, --verbose           Set verbose output\n");
	fprintf(stderr, "  -h, --help              Print this message and exit\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Report bugs to <steven.wright@york.ac.uk>\n");
}

/**
 * @brief Parse the argv arguments passed to the application
 * 
 * @param argc The number of arguments present
 * @param argv An array of the arguments presented
 */
void parse_args(int argc, char *argv[]) {
    int option_index = 0;
    char c;

    while ((c = getopt_long(argc, argv, GETOPTS, long_options, &option_index)) != -1) {
        switch (c) {
			case 'x':
				x = atoi(optarg);
				break;
			case 'y':
				y = atoi(optarg);
				break;
			case 'p':
				num_part_per_dim = atoi(optarg);
				break;
			case 's':
				cell_size = atof(optarg);
				break;
			case 'r':
				r_cut_off = atof(optarg);
				break;
			case 't':
                t_end = atof(optarg);
				break;
			case 'i':
				niters = atof(optarg);
				break;
            case 'd':
                dt = atof(optarg);
                break;
			case 'f':
				output_freq = atoi(optarg);
				break;
			case 'e':
				seed = atol(optarg);
				break;
			case 'n':
				no_output = 1;
				break;
			case 'o':
				set_basename(optarg);
				break;
			case 'c':
				enable_checkpoints = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case '?':
            case 'h':
				print_help(argv[0]);
				exit(1);
        }
    }

	if (r_cut_off > cell_size) {
		fprintf(stderr, "Error: The cell size must be greater than or equal to the cut off distance.\n");
		print_help(argv[0]);
		exit(1);
	}
}

/**
 * @brief Print out the current parameters
 * 
 */
void print_opts() {
    printf("=======================================\n");
    printf("Started with the following options\n");
    printf("=======================================\n");
    printf("  cellx            = %14d\n", x);
	printf("  celly            = %14d\n", y);
   	printf("  cellsize         = %14.12f\n", cell_size);
	printf("  cutoff           = %14.12f\n", r_cut_off);
	printf("  del-t            = %14lf\n", dt);
	printf("  freq             = %14d\n", output_freq);
	printf("  seed             = %14ld\n", seed);
	printf("  endtime          = %14.12lf\n", t_end);
	printf("  noio             = %14d\n", no_output);
	printf("  output           = %s\n", get_basename());
	printf("  checkpoint       = %14d\n", enable_checkpoints);	
    printf("=======================================\n");
}