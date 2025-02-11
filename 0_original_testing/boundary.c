
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
		cells[0][j].head = cells[x][j].head;
		cells[x+1][j].head = cells[1][j].head;
	}

	for (int i = 0; i < x+2; i++) {
		cells[i][0].head = cells[i][y].head;
		cells[i][y+1].head = cells[i][1].head;
	}
}