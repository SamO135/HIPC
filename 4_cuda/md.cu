#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda_runtime.h>

#include "args.h"
#include "boundary.h"
#include "data.h"
#include "setup.h"
#include "vtk.h"

#define N (x+2)


/**
 * @brief This routine calculates the acceleration felt by each particle based on evaluating the Lennard-Jones 
 *        potential with its neighbours. It only evaluates particles within a cut-off radius, and uses cells to 
 *        reduce the search space. It also calculates the potential energy of the system. 
 * 
 * @return double The potential energy
 */

__global__ void comp_accel(struct cell_list* dev_cells, double* results, double* misc_doubles, int* misc_ints) {
	double i_offset;
	double j_offset;
	double ia_offset;
	double jb_offset;

	double p_real_x;
	double p_real_y;
	double q_real_x;
	double q_real_y;

	double dx;
	double dy;
	double r_2;
	
	double r_2_inv;
	double r_6_inv;
	double f;

	int x = misc_ints[0];
	int y = misc_ints[1];
	int num_particles = misc_ints[2];
	double r_cut_off = misc_doubles[0];
	double cell_size = misc_doubles[1];
	double Uc = misc_doubles[2];
	double Duc = misc_doubles[3];
	double r_cut_off_2 = misc_doubles[4];

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid == 0) {return;}

	// zero acceleration for every particle
	// for (int i = 1; i < x+1; i++) {
		int row = tid * (x+2);//i * (x+2);
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < dev_cells[row + j].num_particles; k++) {
				struct particle_t * p = &dev_cells[row + j].particles[k];
				p->ax = 0.0;
				p->ay = 0.0;
			}
		}
	// }

	double pot_energy = 0.0;

	// for (int i = 1; i < x+1; i++) {
		i_offset = (row-1) * cell_size;
		for (int j = 1; j < y+1; j++) {
			j_offset = (j-1) * cell_size;
			for (int k = 0; k < dev_cells[row + j].num_particles; k++) {
				struct particle_t * p = &dev_cells[row + j].particles[k];
				// since particles are stored relative to their cell, calculate the
				// actual x and y coordinates.
				p_real_x = i_offset + p->x;
				p_real_y = j_offset + p->y;
				// Compare each particle with all particles in the 9 cells
				for (int a = -1; a <= 1; a++) {
					int a_row = a * (x+2);
					ia_offset = (row+a_row-1) * cell_size;
					for (int b = -1; b <= 1; b++) {
						jb_offset = (j+b-1) * cell_size;
						for (int c = 0; c < dev_cells[row+a_row + j+b].num_particles; c++) {
							struct particle_t * q = &dev_cells[row+a_row + j+b].particles[c];
							// if p and q are the same particle, skip
							if (p == q) {
								continue;
							}

							// since particles are stored relative to their cell, calculate the
							// actual x and y coordinates.
							q_real_x = ia_offset + q->x;
							q_real_y = jb_offset + q->y;
							
							// calculate distance in x and y, then absolute distance
							dx = p_real_x - q_real_x;
							dy = p_real_y - q_real_y;
							r_2 = dx*dx + dy*dy;
							
							// if distance less than cut off, calculate force and 
							// use this to calculate acceleration in each dimension
							// calculate potential energy of each particle at the same time
							if (r_2 < r_cut_off_2) {
								r_2_inv = 1.0 / r_2;
								r_6_inv = r_2_inv * r_2_inv * r_2_inv;
								
								f = (48.0 * r_2_inv * r_6_inv * (r_6_inv - 0.5));
								p->ax += f*dx;
								p->ay += f*dy;

								pot_energy += 4.0 * r_6_inv * (r_6_inv - 1.0) - Uc - Duc * (sqrt(r_2) - r_cut_off);
							}
						}
					}
				}
			}
		}
	// }
	// return the average potential energy (i.e. sum / number)
	pot_energy = pot_energy / num_particles;
	results[tid] = pot_energy;
}

/**
 * @brief This routine updates the velocity of each particle for half a time step and then 
 *        moves the particle for a whole time step
 * 
 */
void move_particles() {
	// move all particles half a time step
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				// update velocity to obtain v(t + Dt/2)
				p->vx += dth * p->ax;
				p->vy += dth * p->ay;

				// update particle coordinates to p(t + Dt) (scaled to the cell_size)
				p->x += (dt * p->vx);
				p->y += (dt * p->vy);
			}
		}
	}
}

/**
 * @brief This routine updates the cell lists. If a particles coordinates are not within a cell
 *        any more, this function calculates the cell it should be in and performs the move.
 *        If a particle moves more than 1 cell in any direction, this indicates poor settings
 *        and therefore an error is generated.
 * 
 */
void update_cells() {
	// move particles that need to move cell lists
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &((cells[i][j]).particles[k]);
				// if a particles x or y value is greater than the cell size or less than 0, it must have moved cell
				// do a quick check to make sure its not moved 2 cells (since this means our time step is too large, or something else is going wrong)
				if ((p->x < 0.0) | (p->x >= cell_size) | (p->y < 0.0) | (p->y >= cell_size)) {
					if ((p->x < (-cell_size)) || (p->x >= (2*cell_size)) || (p->y < (-cell_size)) || (p->y >= (2*cell_size))) {
						fprintf(stderr, "A particle (%d, %d) (%d) (p%d) has moved more than one cell!\n", i, j, cells[i][j].num_particles, k);
						exit(1);
					}


					// work out whether we've moved a cell in the x and the y dimension
					int x_shift = (p->x < 0.0) ? -1 : (p->x >= cell_size) ? +1 : 0;
					int y_shift = (p->y < 0.0) ? -1 : (p->y >= cell_size) ? +1 : 0;
					
					// the new i and j are +/- 1 in each dimension,
					// but if that means we go out of simulation bounds, wrap it to x and 1
					int new_i = i+x_shift;
					if (new_i == 0) { new_i = x; }
					if (new_i == x+1) { new_i = 1; }
					int new_j = j+y_shift;
					if (new_j == 0) { new_j = y; }
					if (new_j == y+1) { new_j = 1; }

					// update x and y coordinates (i.e. remove the additional cell size)
					p->x = p->x + (x_shift * -cell_size);
					p->y = p->y + (y_shift * -cell_size);

					add_particle(&(cells[new_i][new_j]), p);
					remove_particle(&cells[i][j], k);
				}
			}
		}
	}
}

/**
 * @brief This updates the velocity of particles for the whole time step (i.e. adds the acceleration for another
 *        half step, since its already done half a time step in the move_particles routine). Additionally, this
 *        function calculated the kinetic energy of the system.
 * 
 * @return double The kinetic energy
 */
double update_velocity() {
	double kinetic_energy = 0.0;

	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				// update velocity again by half time to obtain v(t + Dt)
				p->vx += dth * p->ax;
				p->vy += dth * p->ay;

				// calculate the kinetic energy by adding up the squares of the velocities in each dim
				kinetic_energy += (p->vx * p->vx) + (p->vy * p->vy);
			}
		}
	}

	// KE = (1/2)mv^2
	kinetic_energy *= (0.5 / num_particles);
	return kinetic_energy;
}

struct timespec timer;
double get_time() {
	clock_gettime(CLOCK_MONOTONIC, &timer); 
	return (double) (timer.tv_sec + timer.tv_nsec / 1000000000.0);
}

/**
 * @brief This is the main routine that sets up the problem space and then drives the solving routines.
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The exit code of the application
 */
int main(int argc, char *argv[]) {
	double start_time = get_time();
	double elapsed_time = 0;
	// Set default parameters
	set_defaults();
	// parse the arguments 
	parse_args(argc, argv);
	// call set up to update defaults
	setup();

	if (verbose) print_opts();
	
	// set up problem
	problem_setup();

	// apply boundary condition (i.e. update pointers on the boundarys to loop periodically)
	apply_boundary();

	// Create pointers to device memory
	struct cell_list* d_cells;
	double* misc_doubles = (double*) malloc(5 * sizeof(double));
	int* misc_ints = (int*) malloc(3 * sizeof(int));
	double* d_misc_doubles;
	int* d_misc_ints;
	double* d_results;
	double* results = (double*) malloc(N * sizeof(double));

	// Allocate host memory
	misc_doubles[0] = r_cut_off;
	misc_doubles[1] = cell_size;
	misc_doubles[2] = Uc;
	misc_doubles[3] = Duc;
	misc_doubles[4] = r_cut_off_2;

	misc_ints[0] = x;
	misc_ints[1] = y;
	misc_ints[2] = num_particles;

	// Allocate device memory
	cudaMalloc((void **) &d_cells, (x+2) * (y+2) * sizeof(struct cell_list));
	cudaMalloc((void **) &d_misc_doubles, 5 * sizeof(double));
	cudaMalloc((void **) &d_misc_ints, 3 * sizeof(int));
	cudaMalloc((void **) &d_results, N * sizeof(double));


	// Transfer from host to device memory
	cudaMemcpy(d_cells, *cells, (x+2) * (y+2) * sizeof(struct cell_list), cudaMemcpyHostToDevice);
	cudaMemcpy(d_misc_doubles, misc_doubles, 5 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_misc_ints, misc_ints, 3 * sizeof(int), cudaMemcpyHostToDevice);
	

	// Execute kernel
	int block_size = x+2;
	int grid_size = N / block_size;
	comp_accel<<<grid_size, block_size>>>(d_cells, d_results, d_misc_doubles, d_misc_ints);
	cudaDeviceSynchronize();
	
	cudaMemcpy(*cells, d_cells, (x+2) * (y+2) * sizeof(struct cell_list), cudaMemcpyDeviceToHost);

	double potential_energy = 0.0;
	double kinetic_energy = 0.0;

	int iters = 0;
	double t;
	for (t = 0.0; t < t_end; t+=dt, iters++) {
		// move particles half a time step
		move_particles();

		// update cell lists (i.e. move any particles between cell lists if required)
		update_cells();

		// update pointers (because the previous operation might break boundary cell lists)
		apply_boundary();

		// Transfer from host to device memory
		cudaMemcpy(d_cells, *cells, (x+2) * (y+2) * sizeof(struct cell_list), cudaMemcpyHostToDevice);

		// execute kernel
		comp_accel<<<grid_size, block_size>>>(d_cells, d_results, d_misc_doubles, d_misc_ints);
		cudaDeviceSynchronize();
		
		// compute acceleration for each particle and calculate potential energy
		// Transfer 2d array data back to host memory
		// Transfer pot_energy data back to host memory
		cudaMemcpy(*cells, d_cells, (x+2) * (y+2) * sizeof(struct cell_list), cudaMemcpyDeviceToHost);
		cudaMemcpy(results, d_results, N * sizeof(double), cudaMemcpyDeviceToHost);
		potential_energy = 0.0;
		for (int i = 0; i < N; i++)
			potential_energy += results[i];

		// update velocity based on the acceleration and calculate the kinetic energy
		kinetic_energy = update_velocity();

	
		if (iters % output_freq == 0) {
			// calculate temperature and total energy
			double total_energy = kinetic_energy + potential_energy;
			double temp = kinetic_energy * 2.0 / 3.0;

			//calculate elapsed time
			elapsed_time = get_time() - start_time;
			printf("Step %8d, Time: %14.8e (dt: %14.8e), Elapsed time: %lf, Total energy: %14.8e (p:%14.8e,k:%14.8e), Temp: %14.8e\n", iters, t+dt, dt, elapsed_time, total_energy, potential_energy, kinetic_energy, temp);
 
			// if output is enabled and checkpointing is enabled, write out
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+dt);
		}
	}

	// calculate the final energy and write out a final status message
	double final_energy = kinetic_energy + potential_energy;
	elapsed_time = get_time() - start_time;
	printf("Step %8d, Time: %14.8e, Total Time: %lf, Final energy: %14.8e\n", iters, t, elapsed_time, final_energy);
    printf("Simulation complete.\n");

	// if output is enabled, write the mesh file and the final state
	if (!no_output) {
		write_mesh();
		write_result(iters, t);
	}

	// Free device memory
	cudaFree(d_results);
	cudaFree(d_cells);

	free_2d_array((void**) cells);
	return 0;
}

