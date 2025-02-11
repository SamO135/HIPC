#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h> 

#include "setup.h"
#include "data.h"
#include "vtk.h"

/**
 * @brief Set up some default configuration options
 * 
 */
void set_defaults() {
	seed = time(NULL);

	set_default_base();
}

/**
 * @brief Set up configuration options after the arguments have been parsed
 * 
 */
void setup() {
	srand(seed);

	r_cut_off_2 = r_cut_off * r_cut_off;
	r_cut_off_2_inv = 1.0 / r_cut_off_2;
	r_cut_off_6_inv = r_cut_off_2_inv * r_cut_off_2_inv * r_cut_off_2_inv;

	Uc = 4.0 * r_cut_off_6_inv * (r_cut_off_6_inv - 1.0);
	Duc = -48 * r_cut_off_6_inv * (r_cut_off_6_inv - 0.5) / r_cut_off;

	dt = t_end / niters;
	dth = dt / 2.0;
}

/**
 * @brief Set up the problem space, initialise the cells to contain particles,
 *        set the particles to exist on a regular lattice, set their velocities
 *        to be consistent with the initial temperature, but in random orientation.
 * 
 */
void problem_setup() {
	
	// Create a grid of cell lists
	cells = alloc_2d_cell_list_array(x+2, y+2);
	num_particles = x * y * num_part_per_dim * num_part_per_dim;

	double v_sum_x = 0.0;
	double v_sum_y = 0.0;

	// set the normalisation magnitude using the ideal gas law (T = mv^2 / 3)
	double v_magnitude = sqrt(3.0 * init_temp);

	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int a = 0; a < num_part_per_dim; a++) {
				for (int b = 0; b < num_part_per_dim; b++) {
					// set the particles x and y values within the current cell (on a lattice based on number of particles per cell, per dimension)
					double part_x = 0.5 * (1.0 / num_part_per_dim) + ((double) a / num_part_per_dim);
					double part_y = 0.5 * (1.0 / num_part_per_dim) + ((double) b / num_part_per_dim);

					// generate random velocities for the particles, but make sure the overall magnitude is 1.0
					// i.e. generate an angle between 0 and 2*PI then use cos and sin
					double phi = (double) rand() * 2.0 * M_PI / RAND_MAX;
					double rand_vx = cos(phi);
					double rand_vy = sin(phi);

					// create the particle and add it to the current cell list.
					struct particle_t * p = (struct particle_t*) malloc(sizeof(struct particle_t));
					if (p == NULL){
						printf("malloc returned NULL\n");
						exit(EXIT_FAILURE);
					}
					p->x = part_x * cell_size;
					p->y = part_y * cell_size;
					p->vx = rand_vx * v_magnitude;
					p->vy = rand_vy * v_magnitude;
					add_particle(&cells[i][j], p);

					v_sum_x += p->vx;
					v_sum_y += p->vy;
				}
			}	
		}
	}

	// Normalise data to make sure that the total momentum is 0.0 at the start
	double v_avg_x = v_sum_x / num_particles;
	double v_avg_y = v_sum_y / num_particles;

	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				p->vx -= v_avg_x;
				p->vy -= v_avg_y;
			}
		}
	}
}