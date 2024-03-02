#include "../include/atom.h"
#include "../include/point.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// absolute value of a double
double absd(double d) {
        return d > 0 ? d : -d;
}

// prepend an atomnode to a linked list of atomnodes
void prepend_to_atomnode_list(AtomNode **head, Atom *data) {
        AtomNode *new = malloc(sizeof(AtomNode)); 
        new->next = *head;
        new->data = data;
        *head = new;
}

void free_atomnode_list(AtomNode *head) {
        while (head != NULL) {
                AtomNode *temp = head->next;
                free(head);
                head = temp; 
        }
}

Atom *init_atoms(int n_atoms, int x_max, int y_max) {
        srand(100);
        Atom *atoms = malloc(sizeof(Atom) * n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
                atoms[i].mass = MIN_MASS + rand() % (MAX_MASS - MIN_MASS + 1);
                atoms[i].position = (Point) {.x = rand() % (x_max / 2) + ((double) x_max / 4), 
                                             .y = rand() % (y_max / 2) + ((double) y_max / 4) };
        }
        return atoms;
}

void debug_atoms(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                printf("Mass: %d\n", atom->mass);
                printf("Position: %.2g %.2g\n", atom->position.x, atom->position.y);
                printf("Velocity: %.2g %.2g\n", atom->velocity.x, atom->velocity.y);
                printf("Acceleration: %.2g %.2g\n", atom->acceleration.x, atom->acceleration.y);
        }
}

void iterate_kinematics(Atom *atoms, int n_atoms) {
        for (int i = 0; i < n_atoms; ++i) {
                Atom *atom = atoms + i;
                atom->position.x += atom->velocity.x;
                atom->position.y += atom->velocity.y;
                atom->velocity.x += atom->acceleration.x;
                atom->velocity.y += atom->acceleration.y;
                atom->acceleration = (Point) {0};
        }
}

// removes atoms that are a significant length away from the mean
void remove_faraway_atoms(Atom **atoms_pointer, int *n_atoms_pointer) {

        // determine average of all the atoms
        Point center;
        for (int i = 0; i < *n_atoms_pointer; ++i) {
                Atom *atom = &(*atoms_pointer)[i];
                center.x += atom->position.x;
                center.y += atom->position.y;
        }
        center.x /= *n_atoms_pointer;
        center.y /= *n_atoms_pointer;

        // determine standard deviation of the atoms
        Point std_dev;
        for (int i = 0; i < *n_atoms_pointer; ++i) {
                Atom *atom = &(*atoms_pointer)[i];
                std_dev.x += pow(center.x - atom->position.x, 2);
                std_dev.y += pow(center.y - atom->position.y, 2);
        }
        std_dev.x = sqrt(std_dev.x / *n_atoms_pointer);
        std_dev.y = sqrt(std_dev.y / *n_atoms_pointer);
                
        // add atoms that are too far from the mean
        AtomNode *keep_atoms = malloc(sizeof(AtomNode));
        keep_atoms->data = NULL;  // this will be the last element in the list
        int n_atoms_kept = 0;
        for (int i = 0; i < *n_atoms_pointer; ++i) {
                Atom *atom = &(*atoms_pointer)[i];
                if (absd(atom->position.x - center.x) <= std_dev.x * MIN_STD_DEV_FOR_REMOVAL && 
                    absd(atom->position.y - center.y) <= std_dev.y * MIN_STD_DEV_FOR_REMOVAL) {
                        prepend_to_atomnode_list(&keep_atoms, atom);
                        n_atoms_kept++;
                }
        }

        // reallocate the memory and fill the new list of atoms in reverse to preserve the randomness
        Atom *new_atoms = malloc(sizeof(Atom) * n_atoms_kept);
        AtomNode *temp = keep_atoms;  // go to the next one because the first is blank
        int curr_index = n_atoms_kept - 1;
        while (temp->data != NULL) {
                new_atoms[curr_index--] = *temp->data;  // copy the atom to this
                temp = temp->next;
        }

        // free old memory and set new values
        free(*atoms_pointer); 
        free_atomnode_list(keep_atoms);
        *atoms_pointer = new_atoms;
        *n_atoms_pointer = n_atoms_kept;
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
                Atom *atom = atoms + i;
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
                Atom *atom = atoms + i;
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

