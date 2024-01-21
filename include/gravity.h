#include "atom.h"
#include "point.h"

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
#define GRAVITATIONAL_CONSTANT 2
#define GRAVITATIONAL_DISTANCE_GUARD 20000
#define GRAVITATIONAL_FUNCTION apply_gravity_barnes_hut

// this is theta as defined here: https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation#Algorithm
#define BARNES_HUT_THRESHOLD 0.6

typedef struct {
        Point max_values, min_values;
} Box;

typedef struct {
        Point center_of_mass;
        int total_mass;
} MassApprox;

typedef struct quad_tree_node {
        struct quad_tree_node *top_left, *top_right, *bottom_left, *bottom_right;
        int total_mass;
        Box bounds;
        Point center_of_mass;
} QuadTreeNode;

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

// create a quad tree for the Barnes-Hut Algorithm 
QuadTreeNode *build_barnes_hut_tree(Atom *atoms, int n_atoms);

// determine the gravity using the barnes-hut algorithm 
void apply_gravity_barnes_hut(Atom *atoms, int n_atoms);

// free a quad tree 
void free_quad_tree(QuadTreeNode *root);

#endif
