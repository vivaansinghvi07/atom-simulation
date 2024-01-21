#include "point.h"

#ifndef SIM_ATOM
#define SIM_ATOM

// min number of std devs to consider for removing atoms
#define MIN_STD_DEV_FOR_REMOVAL 7

// constants for repulsion force between particles
#define REPULSION_FUNC_A 0.000002
#define REPULSION_FUNC_B 5

// min distance between two centers to not be considered a collision, represents 2 * radius
#define COLLISION_ATOM_WIDTH 8

typedef struct {
        int mass;
        Point position;
        Point velocity;
        Point acceleration;
} Atom; 

typedef struct atom_node {
        Atom *data;
        struct atom_node *next;
} AtomNode;

// absolute value of a double
double absd(double d);

// prepends an atomnode to a linked list 
void prepend_to_atomnode_list(AtomNode **head, Atom *atom);

// free a linked list of atomnodes
void free_atomnode_list(AtomNode *head);

// generates a list of atoms that have posisions (0 -> x_max / 2, 0 -> y_max / 2)
Atom *init_atoms(int n_atoms, int x_max, int y_max);

// prints information for each atom on the screen
void debug_atoms(Atom *atoms, int n_atoms);

// remove atoms that are too far from the rest
void remove_faraway_atoms(Atom **atoms_pointer, int *n_atoms_pointer);

// adds the velocity to the position, and the acceleration to the velocity 
void iterate_kinematics(Atom *atoms, int n_atoms);

// applies collision detection between each atom
// despite the O(n^2) time complexity, this shouldn't matter (up to a point) because operations are very simple
void apply_collision_detection_naive(Atom *atoms, int n_atoms);

// performs a collision between the two atoms, adjusting their velocity
void apply_collision(Atom *atom, Atom *other, Point *distance);

// apply a force of \frac{B}{x^3} - \frac{A}{x^4} between each atom 
// less cool of an effect as the repulsive force tends to make eveyrthing go away
void apply_attraction_repulsion_function_naive(Atom *atoms, int n_atoms);  // O(n^2)

// computes the above force model between two masses 
Point get_attractive_repulsive_force_of_a(int mass_a, Point *com_a, Point *com_b);

#endif
