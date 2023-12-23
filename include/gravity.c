#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "point.c"
#include "atom.c"

#ifndef SIM_GRAVITY 
#define SIM_GRAVITY

// values for the size of arrays in the apply_gravity_approx function 
#define OUTER_GRID_WIDTH 20 
#define INNER_GRID_WIDTH 40
#define OUTER_GRID_LEN 400  // DO NOT FORGET TO UPDATE THESE VALUES
#define INNER_GRID_LEN 1600

// minimum number of "grid widths" away we can go to guarantee a 90% accuracy 
// this is used in the apply_gravity_approx function and determined by ../gravity_approx_opt.py
#define GRID_WIDTH_THRESHOLD 3

// modifications for gravitational functions
#define GRAVITATIONAL_CONSTANT 0.8
#define GRAVITATIONAL_DISTANCE_GUARD 20000
#define GRAVITATIONAL_FUNCTION apply_gravity_approx

typedef struct {
        Point center_of_mass;
        int total_mass;
} MassApprox;

// applies gravity towards the mouse for each point 
void apply_gravity_to_mouse(Atom *atoms, int n_atoms, int mouse_x, int mouse_y, int mouse_mass);

// apply gravity of each atom towards the center of mass of every other one
void apply_gravity_to_center(Atom *atoms, int n_atoms);  // O(n)

// my own algorithm for approximating the gravity to each atom 
void apply_gravity_approx(Atom *atoms, int n_atoms);

// apply gravity on each point by considering each other point
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms);  // O(n^2)

// helper funciton to determine gravity from one object to another
Point get_gravity_of_a(int mass_a, Point *com_a, Point *com_b);

/*
 * FUNCTION IMPLEMENTATIONS BELOW
 * END OF PSEUDO HEADER FILE
 */

int _determine_grid_index(Point *position, Point *min_values, Point *max_values, int grid_size) {
        int y_index = (int) ((position->y - min_values->y) / (max_values->y - min_values->y) * grid_size);
        int x_index = (int) ((position->x - min_values->x) / (max_values->x - min_values->x) * grid_size);
        y_index = y_index == grid_size ? grid_size - 1 : y_index;
        x_index = x_index == grid_size ? grid_size - 1 : x_index;
        return y_index * grid_size + x_index;
}

void apply_gravity_approx(Atom *atoms, int n_atoms) {
        
        // determine the min and max value of all the atoms, which shouldn't be too big b/c of removals
        Point min_values, max_values;
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                double pos_x = atom->position.x;
                double pos_y = atom->position.y;
                
                // values for min and max x and y respectively 
                if (min_values.x > pos_x) {
                        min_values.x = pos_x;
                } else if (max_values.x < pos_x) {
                        max_values.x = pos_x; 
                }
                if (min_values.y > pos_y) {
                        min_values.y = pos_y;
                } else if (max_values.y < pos_y){ 
                        max_values.y = pos_y; 
                }
        }

        // create new arrays to hold information
        MassApprox *outer_grid = malloc(sizeof(MassApprox) * OUTER_GRID_LEN);
        MassApprox *inner_grid = malloc(sizeof(MassApprox) * INNER_GRID_LEN);
        AtomNode **atomnode_grid = malloc(sizeof(AtomNode *) * INNER_GRID_LEN);
        memset(outer_grid, 0, sizeof(MassApprox) * OUTER_GRID_LEN);
        memset(inner_grid, 0, sizeof(MassApprox) * INNER_GRID_LEN);

        // initialize the data of each linked list to be null 
        for (int i = 0; i < INNER_GRID_LEN; ++i) {
                atomnode_grid[i] = malloc(sizeof(AtomNode));
                atomnode_grid[i]->data = NULL;
        }
        
        // create the grids 
        for (int i = 0; i < n_atoms; ++i) {

                // determine the right pointers for everything
                Atom *atom = atoms + i;
                Point *pos = &atom->position;
                int mass = atom->mass;
                int outer_index = _determine_grid_index(pos, &min_values, &max_values, OUTER_GRID_WIDTH);
                int inner_index = _determine_grid_index(pos, &min_values, &max_values, INNER_GRID_WIDTH);

                MassApprox *outer_approx = outer_grid + outer_index;
                MassApprox *inner_approx = inner_grid + inner_index;
                AtomNode **head = atomnode_grid + inner_index;
                
                // adjust the values for the approximations
                outer_approx->total_mass += mass;
                outer_approx->center_of_mass.x += pos->x * mass;
                outer_approx->center_of_mass.y += pos->y * mass;

                inner_approx->total_mass += mass;
                inner_approx->center_of_mass.x += pos->x * mass;
                inner_approx->center_of_mass.y += pos->y * mass;

                // add another atom to check for super close points
                prepend_to_atomnode_list(head, atom);
        }

        // divide each center of mass by the total mass to average it
        for (int i = 0; i < INNER_GRID_LEN; ++i) {
                int total_mass = inner_grid[i].total_mass;
                if (total_mass == 0) {
                        continue;
                }
                inner_grid[i].center_of_mass.x /= total_mass;
                inner_grid[i].center_of_mass.y /= total_mass;
        }
        for (int i = 0; i < OUTER_GRID_LEN; ++i) {
                int total_mass = outer_grid[i].total_mass;
                if (total_mass == 0) {
                        continue;
                }
                outer_grid[i].center_of_mass.x /= total_mass;
                outer_grid[i].center_of_mass.y /= total_mass;
        } 

        // go through each atom and apply the correct gravitational values
        for (int i = 0; i < n_atoms; ++i) { 
                Atom *atom = atoms + i;
                int atom_outer_index = _determine_grid_index(&atom->position, &min_values, &max_values, OUTER_GRID_WIDTH);
                int atom_inner_index = _determine_grid_index(&atom->position, &min_values, &max_values, INNER_GRID_WIDTH);
                int atom_outer_y = atom_outer_index / OUTER_GRID_WIDTH; 
                int atom_outer_x = atom_outer_index % OUTER_GRID_WIDTH;
                int atom_inner_y = atom_inner_index / INNER_GRID_WIDTH; 
                int atom_inner_x = atom_inner_index % INNER_GRID_WIDTH;
                
                // go through the appropriate outer grids, checking if they are far enough
                for (int outer_iy = 0; outer_iy < OUTER_GRID_WIDTH; ++outer_iy) {
                        for (int outer_ix = 0; outer_ix < OUTER_GRID_WIDTH; ++outer_ix) {
                                if (abs(outer_iy - atom_outer_y) < GRID_WIDTH_THRESHOLD &&
                                    abs(outer_ix - atom_outer_x) < GRID_WIDTH_THRESHOLD) {
                                        continue;
                                }
                                MassApprox *approx = outer_grid + outer_iy * OUTER_GRID_WIDTH + outer_ix;
                                Point acceleration_gravity = get_gravity_of_a(approx->total_mass, &approx->center_of_mass, &atom->position);
                                atom->acceleration.x += acceleration_gravity.x;
                                atom->acceleration.y += acceleration_gravity.y;
                        }
                }

                // go through the appropriate inner grids, checking if they are far enough and close enough
                int outer_to_inner_factor = INNER_GRID_WIDTH / OUTER_GRID_WIDTH;
                for (int inner_iy = 0; inner_iy < INNER_GRID_WIDTH; ++inner_iy) {
                        int outer_iy_equivalent = inner_iy / outer_to_inner_factor;
                        for (int inner_ix = 0; inner_ix < INNER_GRID_WIDTH; ++inner_ix) {
                                int outer_ix_equivalent = inner_ix / outer_to_inner_factor;

                                // check if the atom lies in the outer part (already checked)
                                if (abs(outer_iy_equivalent - atom_outer_y) >= GRID_WIDTH_THRESHOLD ||
                                    abs(outer_ix_equivalent - atom_outer_x) >= GRID_WIDTH_THRESHOLD) {
                                        continue;
                                }

                                // check if the atom lies inside the part where it needs to check perpoint 
                                if (abs(inner_iy - atom_inner_y) < GRID_WIDTH_THRESHOLD && 
                                    abs(inner_ix - atom_inner_x) < GRID_WIDTH_THRESHOLD) {
                                        AtomNode *head = *(atomnode_grid + inner_iy * INNER_GRID_WIDTH + inner_ix);
                                        while (head->data != NULL) {
                                                Atom *other = head->data;
                                                if (other != atom) {
                                                        Point acceleration_gravity = get_gravity_of_a(other->mass, &other->position, &atom->position);
                                                        atom->acceleration.x += acceleration_gravity.x;
                                                        atom->acceleration.y += acceleration_gravity.y;
                                                }
                                                head = head->next;
                                        }
                                        continue;
                                }

                                // otherwise, this point is in the zone to be checked against the inner grid
                                MassApprox *approx = inner_grid + inner_iy * INNER_GRID_WIDTH + inner_ix;
                                Point acceleration_gravity = get_gravity_of_a(approx->total_mass, &approx->center_of_mass, &atom->position);
                                atom->acceleration.x += acceleration_gravity.x;
                                atom->acceleration.y += acceleration_gravity.y;
                        }
                }
        }
        
        // free all the allocated things
        for (int i = 0; i < INNER_GRID_LEN; ++i) {
                free_atomnode_list(atomnode_grid[i]);
        }
        free(atomnode_grid);
        free(inner_grid);
        free(outer_grid);
}

// calculates the gravity of each point by considering the gravity of each other point - O(n^2) always
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                for (int j = 0; j < n_atoms; ++j) {
                        if (i == j) {
                                continue;
                        }
                        Atom *other = &atoms[j];
                        Point acceleration_gravity = get_gravity_of_a(other->mass, &other->position, &atom->position);
                        atom->acceleration.x += acceleration_gravity.x;
                        atom->acceleration.y += acceleration_gravity.y;
                }
        }
}

// gravitate points towards given coordinates
void apply_gravity_to_mouse(Atom *atoms, int n_atoms, int mouse_x, int mouse_y, int mouse_mass) {
        Point mouse_coords = (Point) {.x = mouse_x, .y = mouse_y};
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                Point acceleration_gravity = get_gravity_of_a(mouse_mass, &mouse_coords, &atom->position);
                atom->acceleration.x += acceleration_gravity.x;
                atom->acceleration.y += acceleration_gravity.y;
        }
}

// approximation of gravity pointing to the center of mass of everything but each point
void apply_gravity_to_center(Atom *atoms, int n_atoms) {

        // center of mass = \frac{\sum{m_i * r_i}}{\sum{m_i}}
        Point center_of_mass;
        int total_mass = 0; 
        for (int i = 0; i < n_atoms; ++i) {
                int mass = atoms[i].mass;
                center_of_mass.x += atoms[i].position.x * mass;
                center_of_mass.y += atoms[i].position.y * mass;
                total_mass += mass; 
        }
        center_of_mass.x /= total_mass;
        center_of_mass.y /= total_mass;
        
        // to find the distance, first undo the center of mass formula to adjust for the "missing" point 
        // then, use that new center of mass to calculate gravitational acceleration
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                double adjusted_com_x = (center_of_mass.x * total_mass - atom->position.x * atom->mass) / (total_mass - atom->mass);
                double adjusted_com_y = (center_of_mass.y * total_mass - atom->position.y * atom->mass) / (total_mass - atom->mass);

                Point acceleration_gravity = get_gravity_of_a(total_mass - atom->mass, 
                                                              &(Point){.x = adjusted_com_x, .y = adjusted_com_y}, 
                                                              &atom->position);
                atom->acceleration.x += acceleration_gravity.x;
                atom->acceleration.y += acceleration_gravity.y;
        }
}

// apply the following: 
// g = \frac{G * m_i}{r_i^2}
// where r_i is the distance between com_a and com_b
Point get_gravity_of_a(int mass_a, Point *com_a, Point *com_b) {
        Point distance = subtract_a_minus_b(com_a, com_b);
        double overall_distance = abs_point(&distance);
        double acceleration_gravity = mass_a * GRAVITATIONAL_CONSTANT / (overall_distance * overall_distance + GRAVITATIONAL_DISTANCE_GUARD);
        return (Point) {.x = acceleration_gravity * fast_cos_atan2(&distance),
                        .y = acceleration_gravity * fast_sin_atan2(&distance)};
}

#endif
