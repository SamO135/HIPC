
#include "boundary.h"
#include "data.h"

/**
 * @brief Apply the boundary conditions. This effectively points the ghost cell areas
 *        to the same cell list as the opposite edge (i.e. wraps the domain).
 *        This has to be done after every cell list update, just to ensure that a destructive
 *        operations hasn't broken things.
 * 
 */
void apply_boundary() {
	// Apply boundary conditions
	for (int j = 1; j < y+1; j++) {
		cells[0][j] = cells[x][j];
		cells[x+1][j] = cells[1][j];
	}

	for (int i = 0; i < x+2; i++) {
		cells[i][0] = cells[i][y];
		cells[i][y+1] = cells[i][1];
	}
	// Corner cells
	// cells[0][0] = cells[x][y];
	// cells[x+1][0] = cells[1][y];
	// cells[0][y+1] = cells[x][1];
	// cells[x+1][y+1] = cells[1][1];
	
	
	// for (int x = 1; x < cols-1; x++){
	// 	cells[rows-1][x] = cells[1][x];
	// 	cells[0][x] = cells[rows-2][x];
	// }

	// for (int y = 1; y < rows-1; y++){
	// 	cells[y][cols-1] = cells[y][1];
	// 	cells[y][0] = cells[y][cols-2];
	// }
	// cells[0][0] = cells[rows-2][cols-2];
	// cells[rows-1][0] = cells[1][cols-2];
	// cells[0][cols-1] = cells[rows-2][1];
	// cells[rows-1][cols-1] = cells[1][1];
}

// [[14, 12, 13, 14, 12], 
//  [ 9,  7,  8,  9,  7], 
//  [14, 12, 13, 14, 12], 
//  [ 9,  7,  8,  9,  7]]