#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "point.c"

#ifndef SIM_ATOM
#define SIM_ATOM

// constants for repulsion force between particles
#define REPULSION_FUNC_A 0.000002
#define REPULSION_FUNC_B 5

// min distance between two centers to not be considered a collision
#define COLLISION_ATOM_WIDTH 15 

// modifications for gravitational functions
#define GRAVITATIONAL_CONSTANT 0.8
#define GRAVITATIONAL_DISTANCE_GUARD 20000
#define GRAVITATIONAL_FUNCTION apply_gravity_each_point_naive

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

int SIMULATION_STEPS = 0;

// generates a list of atoms that have posisions (0 -> x_max / 2, 0 -> y_max / 2)
Atom *init_atoms(int n_atoms, int x_max, int y_max);

// prints information for each atom on the screen
void debug_atoms(Atom *atoms, int n_atoms);

// goes to the next simulation step, modifying acceleration and velocity when applicable
void step_simulation(Atom *atoms, int n_atoms);

// adds the velocity to the position, and the acceleration to the velocity 
void iterate_kinematics(Atom *atoms, int n_atoms);

// applies gravity towards the mouse for each point 
void apply_gravity_to_mouse(Atom *atoms, int n_atoms, int mouse_x, int mouse_y, int mouse_mass);

// apply gravity of each atom towards the center of mass of every other one
void apply_gravity_to_center(Atom *atoms, int n_atoms);  // O(n)

// apply gravity on each point by considering each other point
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms);  // O(n^2)

// helper funciton to determine gravity from one object to another
Point get_gravity_of_a(int mass_a, Point *com_a, Point *com_b);

// applies collision detection between each atom
void apply_collision_detection_naive(Atom *atoms, int n_atoms); // O(n^2) 

// performs a collision between the two atoms, adjusting their velocity
void apply_collision(Atom *atom, Atom *other, Point *distance);

// apply a force of \frac{B}{x^3} - \frac{A}{x^4} between each atom 
// less cool of an effect as the repulsive force tends to make eveyrthing go away
void apply_attraction_repulsion_function_naive(Atom *atoms, int n_atoms);  // O(n^2)

// computes the above force model between two masses 
Point get_attractive_repulsive_force_of_a(int mass_a, Point *com_a, Point *com_b);

/*
 * FUNCTION IMPLEMENTATIONS BELOW
 * END OF PSEUDO HEADER FILE
 */

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

void iterate_kinematics(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                atom->position.x += atom->velocity.x;
                atom->position.y += atom->velocity.y;
                atom->velocity.x += atom->acceleration.x;
                atom->velocity.y += atom->acceleration.y;
        }
}

void step_simulation(Atom *atoms, int n_atoms) {
        SIMULATION_STEPS++;
        GRAVITATIONAL_FUNCTION(atoms, n_atoms);
        if (SIMULATION_STEPS % 7 == 0) {  // collision detection step
                apply_collision_detection_naive(atoms, n_atoms);
        }
}

// helper function for the below function, changes an atom's velocity for a collision
void _apply_new_velocity(Atom *atom, Point *velocity_com, double angle_of_collision) {

        // determine details about the velocity relative to center of mass
        Point velocity_atom_com = (Point) {.x = atom->velocity.x - velocity_com->x, 
                                           .y = atom->velocity.y - velocity_com->y};
        double angle_of_velocity = atan2(velocity_atom_com.y, velocity_atom_com.x);

        // determine the new angle at which the point bounces off of 
        angle_of_velocity = fmod(angle_of_velocity + M_PI, M_PI * 2); 
        angle_of_collision -= M_PI;
        while (angle_of_velocity - angle_of_collision > M_PI_2) {
                angle_of_collision += M_PI;
        }
        double new_angle_of_velocity = angle_of_velocity - 2 * (angle_of_velocity - angle_of_collision);

        // adjust the atom's velocity with the new information
        atom->velocity = (Point) {.x = abs_point(&velocity_atom_com) * cos(new_angle_of_velocity) + velocity_com->x,
                                  .y = abs_point(&velocity_atom_com) * sin(new_angle_of_velocity) + velocity_com->y};
}

// use conservation of momentum on the two atoms to simulate an elastic collision
void apply_collision(Atom *atom, Atom *other, Point *distance) {

        // determine the velocity of the center of mass and the line of collision between the two objects
        Point velocity_com = (Point) {.x = (atom->velocity.x * atom->mass + other->velocity.x * other->mass) / (atom->mass + other->mass), 
                                      .y = (atom->velocity.y * atom->mass + other->velocity.y * other->mass) / (atom->mass + other->mass)};
        double angle_of_collision = fmod(atan2(distance->y, distance->x) + M_PI, M_PI);
        
        // change the velocities of both of the atom 
        _apply_new_velocity(atom, &velocity_com, angle_of_collision);
        _apply_new_velocity(other, &velocity_com, angle_of_collision);
}

// adds collision detection to stop atoms from going into each other
void apply_collision_detection_naive(Atom *atoms, int n_atoms) {
        Point distance;
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                for (int j = i + 1; j < n_atoms; ++j) {
                        Atom *other = &atoms[j];
                        distance = subtract_a_minus_b(&other->position, &atom->position);
                        if (abs_point(&distance) < COLLISION_ATOM_WIDTH) {
                                apply_collision(atom, other, &distance);
                        }
                } 
        }
}

// using a model of \frac{B}{x^3} - \frac{A}{x^4} to determine force between particles
Point get_attractive_repulsive_force_of_a(int mass_a, Point *com_a, Point *com_b) {
        Point distance = subtract_a_minus_b(com_a, com_b);
        double overall_distance = abs_point(&distance);
        double force = - REPULSION_FUNC_A / pow(overall_distance, 4)  
                       + REPULSION_FUNC_B / pow(overall_distance, 3);
        return (Point) {.x = fast_cos_atan2(&distance) * force,
                        .y = fast_sin_atan2(&distance) * force};
}

void apply_attraction_repulsion_function_naive(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                atom->acceleration = (Point) {0};
                for (int j = 0; j < n_atoms; ++j) {
                        if  (i == j) { 
                                continue;
                        }
                        Atom *other = &atoms[j];
                        Point force = get_attractive_repulsive_force_of_a(other->mass, &other->position, &atom->position);
                        atom->acceleration.x += force.x / atom->mass;
                        atom->acceleration.y += force.y / atom->mass;
                }
        } 
}

// calculates the gravity of each point by considering the gravity of each other point - O(n^2) always
void apply_gravity_each_point_naive(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = &atoms[i];
                atom->acceleration = (Point) {0};
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
                Atom *atom = &atoms[i];
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
        
        Point distance = subtract_a_minus_b(com_a, com_b);
        double overall_distance = abs_point(&distance);
        double acceleration_gravity = mass_a * GRAVITATIONAL_CONSTANT / (overall_distance * overall_distance + GRAVITATIONAL_DISTANCE_GUARD);
        return (Point) {.x = acceleration_gravity * fast_cos_atan2(&distance),
                        .y = acceleration_gravity * fast_sin_atan2(&distance)};
}

#endif
