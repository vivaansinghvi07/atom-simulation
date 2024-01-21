#include "../include/point.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

