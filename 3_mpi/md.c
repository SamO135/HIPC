#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "args.h"
#include "boundary.h"
#include "data.h"
#include "setup.h"
#include "vtk.h"

int rank, size;
int moved;

/**
 * @brief This routine calculates the acceleration felt by each particle based on evaluating the Lennard-Jones 
 *        potential with its neighbours. It only evaluates particles within a cut-off radius, and uses cells to 
 *        reduce the search space. It also calculates the potential energy of the system. 
 * 
 * @return double The potential energy
 */
double comp_accel() {
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

	//calculate each rank's start and end indices
	int starti = rank * (x/size) + 1;
    int endi = (rank+1) * (x/size) + 1;
	if (rank == size-1){
		endi = x+1;
	}

	// zero acceleration for every particle
	for (int i = starti; i < endi; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				p->ax = 0.0;
				p->ay = 0.0;
			}
		}
	}

	double pot_energy = 0.0;


	for (int i = starti; i < endi; i++) {
		i_offset = (i-1) * cell_size;
		for (int j = 1; j < y+1; j++) {
			j_offset = (j-1) * cell_size;
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				// since particles are stored relative to their cell, calculate the
				// actual x and y coordinates.
				p_real_x = i_offset + p->x;
				p_real_y = j_offset + p->y;
				// Compare each particle with all particles in the 9 cells
				for (int a = -1; a <= 1; a++) {
					ia_offset = (i+a-1) * cell_size;
					for (int b = -1; b <= 1; b++) {
						jb_offset = (j+b-1) * cell_size;
						for (int c = 0; c < cells[i+a][j+b].num_particles; c++) {
							struct particle_t * q = &cells[i+a][j+b].particles[c];
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
	}
	MPI_Allreduce(MPI_IN_PLACE, &pot_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// return the average potential energy (i.e. sum / number)
	return pot_energy / num_particles;
}

/**
 * @brief This routine updates the velocity of each particle for half a time step and then 
 *        moves the particle for a whole time step
 * 
 */
void move_particles() {
	//calculate each rank's start and end indices
	int starti = rank * (x/size) + 1;
    int endi = (rank+1) * (x/size) + 1;
	if (rank == size-1){
		endi = x+1;
	}

	// move all particles half a time step
	for (int i = starti; i < endi; i++) {
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
	//calculate each rank's start and end indices
	int starti = rank * (x/size) + 1;
    int endi = (rank+1) * (x/size) + 1;
	if (rank == size-1){
		endi = x+1;
	}

	int above = (rank - 1) < 0 ? size-1 : rank - 1;
	int below = (rank + 1) >= size ? 0 : rank + 1;

	// move particles that need to move cell lists
	for (int i = starti-1; i < endi+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &((cells[i][j]).particles[k]);

				// if a particles x or y value is greater than the cell size or less than 0, it must have moved cell
				// do a quick check to make sure its not moved 2 cells (since this means our time step is too large, or something else is going wrong)
				int x_shift = 0;
				int y_shift = 0;
				if ((p->x < 0.0) | (p->x >= cell_size) | (p->y < 0.0) | (p->y >= cell_size)) {
					if ((p->x < (-cell_size)) || (p->x >= (2*cell_size)) || (p->y < (-cell_size)) || (p->y >= (2*cell_size))) {
						fprintf(stderr, "A particle has moved more than one cell!\n");
						exit(1);
					}


					// work out whether we've moved a cell in the x and the y dimension
					x_shift = (p->x < 0.0) ? -1 : (p->x >= cell_size) ? +1 : 0;
					y_shift = (p->y < 0.0) ? -1 : (p->y >= cell_size) ? +1 : 0;
					
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
					
					moved = 1;

					//if the particle in the current cell is not in the top row and moving up OR in 
					// the bottom row and moving down (i.e. not moving out of range)
					if (!((i == starti-1 && x_shift == -1) || (i == endi && x_shift == 1))){
						add_particle(&(cells[new_i][new_j]), p);
						remove_particle(&cells[i][j], k);
					}

					k--;
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
	//calculate each rank's start and end indices
	int starti = rank * (x/size) + 1;
    int endi = (rank+1) * (x/size) + 1;
	if (rank == size-1){
		endi = x+1;
	}

	double kinetic_energy = 0.0;

	for (int i = starti; i < endi; i++) {
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
	MPI_Allreduce(MPI_IN_PLACE, &kinetic_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// KE = (1/2)mv^2
	kinetic_energy *= (0.5 / num_particles);
	return kinetic_energy;
}


/**
 * @brief This is the main routine that sets up the problem space and then drives the solving routines.
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The exit code of the application
 */
int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double start_time = MPI_Wtime();

	// calculate each rank's start and end indices
	int starti = rank * (x/size) + 1;
    int endi = (rank+1) * (x/size) + 1;
	if (rank == size-1){
		endi = x+1;
	}

	// create MPI datatype for particle_t
	MPI_Datatype mpi_particle_t;
	//calculate displacements
	MPI_Aint part_displacements[7];
	struct particle_t dummy_particle;
	MPI_Aint part_base_address;
	MPI_Get_address(&dummy_particle, &part_base_address);
	MPI_Get_address(&dummy_particle.x, &part_displacements[0]);
	MPI_Get_address(&dummy_particle.y, &part_displacements[1]);
	MPI_Get_address(&dummy_particle.ax, &part_displacements[2]);
	MPI_Get_address(&dummy_particle.ay, &part_displacements[3]);
	MPI_Get_address(&dummy_particle.vx, &part_displacements[4]);
	MPI_Get_address(&dummy_particle.vy, &part_displacements[5]);
	MPI_Get_address(&dummy_particle.part_id, &part_displacements[6]);
	for (int i = 0; i < 7; i++){
		part_displacements[i] = MPI_Aint_diff(part_displacements[i], part_base_address);
	}

	int part_blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype part_types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Type_create_struct(7, part_blocklengths, part_displacements, part_types, &mpi_particle_t);
	MPI_Type_commit(&mpi_particle_t);


	// create MPI datatype for cell_list
	MPI_Datatype mpi_cell_list;
	// calculate the displacements
	MPI_Aint cell_displacements[2];
	struct cell_list dummy_cell;
	MPI_Aint cell_base_address;
	MPI_Get_address(&dummy_cell, &cell_base_address);
	MPI_Get_address(&dummy_cell.particles, &cell_displacements[0]);
	MPI_Get_address(&dummy_cell.num_particles, &cell_displacements[1]);
	cell_displacements[0] = MPI_Aint_diff(cell_displacements[0], cell_base_address);
	cell_displacements[1] = MPI_Aint_diff(cell_displacements[1], cell_base_address);

	int cell_blocklengths[2] = {10, 1};
	MPI_Datatype cell_types[2] = {mpi_particle_t, MPI_INT};
	MPI_Type_create_struct(2, cell_blocklengths, cell_displacements, cell_types, &mpi_cell_list);
	MPI_Type_commit(&mpi_cell_list);

	// create MPI datatype for a row of cell_lists
	MPI_Datatype ghost_row;
	MPI_Type_vector(y, 1, 1, mpi_cell_list, &ghost_row);
	MPI_Type_commit(&ghost_row);

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
	
	comp_accel();

	double potential_energy = 0.0;
	double kinetic_energy = 0.0;

	int iters = 0;
	double t;
	for (t = 0.0; t < t_end; t+=dt, iters++) {

		int above = (rank - 1) < 0 ? size-1 : rank - 1;
		int below = (rank + 1) >= size ? 0 : rank + 1;
		
		// move particles half a time step
		move_particles();

		MPI_Sendrecv(&(cells[starti][1]), 1, ghost_row, above, 0, &(cells[endi][1]), 1, ghost_row, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(cells[endi-1][1]), 1, ghost_row, below, 0, &(cells[starti-1][1]), 1, ghost_row, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// update cell lists (i.e. move any particles between cell lists if required)
		update_cells();

		MPI_Sendrecv(&(cells[starti][1]), 1, ghost_row, above, 0, &(cells[endi][1]), 1, ghost_row, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Sendrecv(&(cells[endi-1][1]), 1, ghost_row, below, 0, &(cells[starti-1][1]), 1, ghost_row, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// update pointers (because the previous operation might break boundary cell lists)
		apply_boundary();
		
		// compute acceleration for each particle and calculate potential energy
		potential_energy = comp_accel();

		// update velocity based on the acceleration and calculate the kinetic energy
		kinetic_energy = update_velocity();

	
		if (iters % output_freq == 0 && rank == size-1) {
			// calculate temperature and total energy
			double total_energy = kinetic_energy + potential_energy;
			double temp = kinetic_energy * 2.0 / 3.0;

			//calculate elapsed time
			elapsed_time = MPI_Wtime() - start_time;
			printf("Step %8d, Time: %14.8e (dt: %14.8e), Elapsed time: %lf, Total energy: %14.8e (p:%14.8e,k:%14.8e), Temp: %14.8e\n", iters, t+dt, dt, elapsed_time, total_energy, potential_energy, kinetic_energy, temp);
 
			// if output is enabled and checkpointing is enabled, write out
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+dt);
		}
	}

	if (rank == 0){
		// calculate the final energy and write out a final status message
		double final_energy = kinetic_energy + potential_energy;
		elapsed_time = MPI_Wtime() - start_time;
		printf("Step %8d, Time: %14.8e, Total Time: %lf, Final energy: %14.8e\n", iters, t, elapsed_time, final_energy);
		printf("Simulation complete.\n");

		// if output is enabled, write the mesh file and the final state
		if (!no_output) {
			write_mesh();
			write_result(iters, t);
		}
	}

	free_2d_array((void**) cells);
	MPI_Type_free(&ghost_row);
	MPI_Type_free(&mpi_particle_t);
	MPI_Type_free(&mpi_cell_list);
	
	MPI_Finalize();
	
	return 0;
}

