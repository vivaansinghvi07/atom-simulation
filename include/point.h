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

#endif
