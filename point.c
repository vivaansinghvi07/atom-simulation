#include <math.h>
#include <stdlib.h>

#ifndef SIM_POINT
#define SIM_POINT

typedef struct point { 
        double x;
        double y;
} Point;

typedef struct point_node {
        Point data;
        struct point_node *next;
} PointNode;

double abs_point(Point *p);
void prepend_to_list(PointNode **head, Point *data);
void free_list(PointNode *head);

// get the absolute value of a point
double abs_point(Point *p) {
        return sqrt(p->x * p->x + p->y * p->y);
}

// creates a new node and adds it to the beginning of the list, replacing the head
void prepend_to_list(PointNode **head, Point *data) {
        PointNode *new = malloc(sizeof(PointNode)); 
        new->next = *head;
        new->data = *data;
        *head = new;
}

// free a linked list
void free_list(PointNode *head) {
        while (head != NULL) {
                PointNode *temp = head->next;
                free(head);
                head = temp; 
        }
}

#endif
