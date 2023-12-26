#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef SIM_POINT
#define SIM_POINT

typedef struct { 
        double x;
        double y;
} Point;

typedef struct point_node {
        Point data;
        struct point_node *next;
} PointNode;

// computes \cos{(\tan^{-1}{y/x})} as \frac{x}{\sqrt{x^2 + y^2}}
double fast_cos_atan2(Point *p);

// computes \sin{(\tan^{-1}{y/x})} as \frac{y}{\sqrt{x^2 + y^2}}
double fast_sin_atan2(Point *p);

// gets the absolute value of a point relative to the origin
double abs_point(Point *p);

// returns a new point that is the difference of a and b
Point subtract_a_minus_b(Point *a, Point *b);

// prepends a Point to a linked list of PointNodes
void prepend_to_pointnode_list(PointNode **head, Point *data);

// frees a list of PointNodes
void free_pointnode_list(PointNode *head);

/*
 * FUNCTION IMPLEMENTATIONS BELOW
 * END OF PSEUDO HEADER FILE
 */

// does \cos{(\tan^{-1}{x})} for point quick
double fast_cos_atan2(Point *p) {
        double dist = abs_point(p);
        if (dist == 0) {
                return 0;
        } 
        return p->x / dist;
}

// does \sin{(\tan^{-1}{x})} for point quick
double fast_sin_atan2(Point *p) {
        double dist = abs_point(p);
        if (dist == 0) {
                return 0;
        }
        return p->y / dist;
}

// get the absolute value of a point
double abs_point(Point *p) {
        return sqrt(p->x * p->x + p->y * p->y);
}

// denotes the subtraction operation of two points
Point subtract_a_minus_b(Point *a, Point *b) {
        return (Point) {.x = a->x - b->x, .y = a->y - b->y};
}

// creates a new node and adds it to the beginning of the list, replacing the head
void prepend_to_pointnode_list(PointNode **head, Point *data) {
        PointNode *new = malloc(sizeof(PointNode)); 
        new->next = *head;
        new->data = *data;
        *head = new;
}

// free a linked list
void free_pointnode_list(PointNode *head) {
        while (head != NULL) {
                PointNode *temp = head->next;
                free(head);
                head = temp; 
        }
}

#endif
