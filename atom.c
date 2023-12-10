#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// fine tune these parameters to control gravitation 
#define GRAVITATIONAL_CONSTANT 0.0005
#define GRAVITATIONAL_DISTANCE_GUARD 20000

typedef struct point { 
        double x;
        double y;
} Point;

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
void step_simulation(Atom *atoms, int n_atoms);

Atom *init_atoms(int n_atoms, int x_max, int y_max) {
        srand(100);
        Atom *atoms = malloc(sizeof(Atom) * n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
                atoms[i].mass = 1; 
                atoms[i].position = (Point) {rand() % x_max, rand() % y_max};
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

// google "pauli repulsive force"
void apply_repulsive_force(Atom *atoms, int n_atoms) {

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
                        Point distance = {.x = other->position.x - atom->position.x, .y = other ->position.y - atom->position.y};
                        double angle_to_other = atan2(distance.y, distance.x);
                        double overall_distance = sqrt(distance.y * distance.y + distance.x * distance.x);
                        double acceleration_gravity = other->mass * GRAVITATIONAL_CONSTANT / (overall_distance * overall_distance + GRAVITATIONAL_DISTANCE_GUARD);
                        atom->acceleration.x += acceleration_gravity * cos(angle_to_other);
                        atom->acceleration.y += acceleration_gravity * sin(angle_to_other);
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
        
        // g = \frac{G * m_i}{r_i^2} // in this case, a constant is being added
        // to find the distance, first undo the center of mass formula to adjust for the "missing" point 
        // then, use that new center of mass to calculate 
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                double adjusted_com_x = (center_of_mass.x * total_mass - atom->position.x * atom->mass) / (total_mass - atom->mass);
                double adjusted_com_y = (center_of_mass.y * total_mass - atom->position.y * atom->mass) / (total_mass - atom->mass);
                Point distance = {.x = adjusted_com_x - atom->position.x, .y = adjusted_com_y - atom->position.y};
                double angle_to_com = atan2(distance.y, distance.x);
                double overall_distance = sqrt(distance.y * distance.y + distance.x * distance.x);
                double acceleration_gravity = (total_mass - atom->mass) * GRAVITATIONAL_CONSTANT / (overall_distance * overall_distance + GRAVITATIONAL_DISTANCE_GUARD);
                atom->acceleration.x = acceleration_gravity * cos(angle_to_com);
                atom->acceleration.y = acceleration_gravity * sin(angle_to_com);
        }
}
