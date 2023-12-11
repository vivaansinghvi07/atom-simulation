#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "point.c"

#ifndef SIM_ATOM
#define SIM_ATOM

// fine tune these parameters to control gravitation 
#define GRAVITATIONAL_CONSTANT 0.5
#define GRAVITATIONAL_DISTANCE_GUARD 20000

typedef struct atom {
        int mass;
        Point position;
        Point velocity;
        Point acceleration;
} Atom; 

Atom *init_atoms(int n_atoms, int x_max, int y_max);
void debug_atoms(Atom *atoms, int n_atoms);
void apply_gravity_to_center(Atom *atoms, int n_atoms);
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms);
Point get_gravity_of_a(int mass_a, Point *com_a, Point *com_b);
void step_simulation(Atom *atoms, int n_atoms);

Atom *init_atoms(int n_atoms, int x_max, int y_max) {
        srand(100);
        Atom *atoms = malloc(sizeof(Atom) * n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
                        atoms[i].mass = 1; 
                atoms[i].position = (Point) {.x = rand() % (x_max / 2) + ((double) x_max / 4), 
                                             .y = rand() % (y_max / 2) + ((double) y_max / 4) };
        }
        return atoms;
}

void debug_atoms(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                printf("Mass: %d\n", atom->mass);
                printf("Position: %.2g %.2g\n", atom->position.x, atom->position.y);
                printf("Velocity: %.2g %.2g\n", atom->velocity.x, atom->velocity.y);
                printf("Acceleration: %.2g %.2g\n", atom->acceleration.x, atom->acceleration.y);
        }
}

void step_simulation(Atom *atoms, int n_atoms) {

        apply_gravity_to_center(atoms, n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                atom->position.x += atom->velocity.x;
                atom->position.y += atom->velocity.y;
                atom->velocity.x += atom->acceleration.x;
                atom->velocity.y += atom->acceleration.y;
        }
}

// calculates the gravity of each point by considering the gravity of each other point - O(n^2) always
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
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
                Atom *atom = &atoms[i];
                double adjusted_com_x = (center_of_mass.x * total_mass - atom->position.x * atom->mass) / (total_mass - atom->mass);
                double adjusted_com_y = (center_of_mass.y * total_mass - atom->position.y * atom->mass) / (total_mass - atom->mass);

                Point acceleration_gravity = get_gravity_of_a(total_mass - atom->mass, 
                                                              &(Point){.x = adjusted_com_x, .y = adjusted_com_y}, 
                                                              &atom->position);
                atom->acceleration.x = acceleration_gravity.x;
                atom->acceleration.y = acceleration_gravity.y;
        }
}

// apply the following: 
// g = \frac{G * m_i}{r_i^2}
// where r_i is the distance between com_a and com_b
Point get_gravity_of_a(int mass_a, Point *com_a, Point *com_b) {
        
        Point distance = {.x = com_a->x - com_b->x, .y = com_a->y - com_b->y};
        double angle_to_a = atan2(distance.y, distance.x); 
        double overall_distance = sqrt(distance.y * distance.y + distance.x * distance.x);
        double acceleration_gravity = mass_a * GRAVITATIONAL_CONSTANT / (overall_distance * overall_distance + GRAVITATIONAL_DISTANCE_GUARD);
        return (Point) { .x = acceleration_gravity * cos(angle_to_a), .y = acceleration_gravity * sin(angle_to_a) };
}

#endif
