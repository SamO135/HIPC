#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "vtk.h"
#include "data.h"

char checkpoint_basename[1024];
char result_filename[1024];
char mesh_filename[1024];

/**
 * @brief Set the default basename for file output to out/vortex
 * 
 */
void set_default_base() {
    set_basename((char*) "out/md");
}

/**
 * @brief Set the basename for file output
 * 
 * @param base Basename string
 */
void set_basename(char *base) {
    checkpoint_basename[0] = '\0';
    result_filename[0] = '\0';
    sprintf(checkpoint_basename, "%s-%%d.vtp", base);
    sprintf(result_filename, "%s.vtp", base);
	sprintf(mesh_filename, "%s-mesh.vti", base);
}

/**
 * @brief Get the basename for file output
 * 
 * @return char* Basename string
 */
char *get_basename() {
    return checkpoint_basename;
}

/**
 * @brief Write a checkpoint VTK file (with the iteration number in the filename)
 * 
 * @param iteration The current iteration number
 * @return int Return whether the write was successful
 */
int write_checkpoint(int iters, double t) { 
    char filename[1024];
    sprintf(filename, checkpoint_basename, iters);
    return write_vtk(filename, iters, t);
}

/**
 * @brief Write the final output to a VTK file
 * 
 * @return int Return whether the write was successful
 */
int write_result(int iters, double t) {
    return write_vtk(result_filename, iters, t);
}

/**
 * @brief Write out a particle VTK file (i.e. a .vtp file).
 * 
 * @param filename The filename to use for output
 * @param iters The number of iterations
 * @param t The simulation time
 * @return int Return whether the write was successful
 */
int write_vtk(char * filename, int iters, double t) {
	FILE * f = fopen(filename, "w");
    if (f == NULL) {
        perror("Error");
        return -1;
    }
	
	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(f, "<PolyData>\n");
	fprintf(f, "<FieldData>\n");
    fprintf(f, "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">\n");
    fprintf(f, "%.12e\n", t);
    fprintf(f, "</DataArray>\n");
	fprintf(f, "<DataArray type=\"Int32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">\n");
    fprintf(f, "%d\n", iters);
    fprintf(f, "</DataArray>\n");
    fprintf(f, "</FieldData>\n");
	fprintf(f, "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfCells=\"0\">\n", num_particles);
	fprintf(f, "<Points>\n");
	fprintf(f, "<DataArray type=\"Float64\" Name=\"particles\" NumberOfComponents=\"3\" format=\"ascii\">\n");
	for (int i = 1; i < x+1; i++) {
		for (int j = 1; j < y+1; j++) {
			for (int k = 0; k < cells[i][j].num_particles; k++) {
				struct particle_t * p = &cells[i][j].particles[k];
				double p_real_x = ((i-1) * cell_size) + p->x;
				double p_real_y = ((j-1) * cell_size) + p->y;
				fprintf(f, "%.12e %.12e 0 \n", p_real_x, p_real_y);
			}
		}
	}
	
	fprintf(f, "\n</DataArray>\n");
	fprintf(f, "</Points>\n");
	fprintf(f, "</Piece>\n");
	fprintf(f, "</PolyData>\n");
	fprintf(f, "</VTKFile>\n");
	fclose(f);
	return 0;
}

/**
 * @brief Write out the mesh VTK file (i.e. a .vti file). This mesh can be plotted
 *        alongside the particle data to show the data grid as an overlay.
 * 
 * @return int Return whether the write was successful
 */
int write_mesh() {
	FILE * f = fopen(mesh_filename, "w");
    if (f == NULL) {
        perror("Error");
        return -1;
    }

	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(f, "<ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"%lf %lf 0\">\n", x, y, cell_size, cell_size);
	fprintf(f, "<Piece Extent=\"0 %d 0 %d 0 0\">\n", x, y);
	fprintf(f, "<CellData></CellData>\n");
	fprintf(f, "<PointData></PointData>\n");
	fprintf(f, "<Points></Points>\n");
	fprintf(f, "</Piece>\n");
	fprintf(f, "</ImageData>\n");
	fprintf(f, "</VTKFile>\n");

	fclose(f);
	return 0;
}
