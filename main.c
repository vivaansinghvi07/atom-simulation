#include <SDL2/SDL.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "include/point.c"
#include "include/atom.c"

#define SCREEN_X 1600
#define SCREEN_Y 900
#define MOUSE_MASS 200
#define ATOM_DISPLAY_WIDTH 1
#define DISPLAY_COLOR COLOR_VELOCITY
#define CLICK_PLACE_WIDTH 25
#define CLICK_PLACE_GAP 5

int N_ATOMS = 20000;

typedef struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
} RGB_color;

enum ColorMode {
        COLOR_NONE, 
        COLOR_RANDOM, 
        COLOR_VELOCITY
};

void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color);
void clear_screen(SDL_Surface *surface);
void display_atoms(SDL_Surface *surface, Atom *atoms, enum ColorMode color_mode);
void add_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y);
int main(void);

// https://github.com/MikeShah/SDL2_Tutorials/blob/main/8_ModifyingSurface/main.cpp
void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color) {
        if (x < 0 || x >= SCREEN_X || y < 0 || y >= SCREEN_Y) {
                return;
        }
        uint8_t *pixels = (uint8_t *)surface->pixels;
        int pixel_address = y * surface->pitch + x * surface->format->BytesPerPixel;
        pixels[pixel_address + 0] = color.b;
        pixels[pixel_address + 1] = color.g;
        pixels[pixel_address + 2] = color.r;
}

void clear_screen(SDL_Surface *surface) {
        uint8_t *pixels = (uint8_t *)surface->pixels; 
        memset(pixels, 0, SCREEN_X * SCREEN_Y * 4);  // stored as gbra
}

void display_atoms(SDL_Surface *surface, Atom *atoms, enum ColorMode color_mode) {

        // by setting the seed here, we keep the same random colors for each atom
        srand(100);

        // determine minimum and maximum value of velocity for the relative coloring
        double min_velocity, max_velocity;
        if (color_mode == COLOR_VELOCITY) {
                for (int i = 0; i < N_ATOMS; ++i) {
                        double velocity = abs_point(&atoms[i].velocity);
                        if (velocity < min_velocity) {
                                min_velocity = velocity;
                        } else if (velocity > max_velocity) { 
                                max_velocity = velocity;
                        }
                }
        }

        for (int i = 0; i < N_ATOMS; ++i) { 
                // store atom position more conveniently
                int atom_x = atoms[i].position.x;
                int atom_y = atoms[i].position.y;

                // determine what the color will be=
                int r = 255, g = 255, b = 255;
                switch (color_mode) {
                        case COLOR_NONE:
                                break;
                        case COLOR_RANDOM:
                                r = rand() % 255;
                                g = rand() % 255;
                                b = rand() % 255;
                                break;
                        case COLOR_VELOCITY:
                                g = (int) ((abs_point(&atoms[i].velocity) + min_velocity) / (max_velocity - min_velocity) * 255);
                                break;
                }

                // apply coloring to each pixel for the atom
                for (int x = -(ATOM_DISPLAY_WIDTH / 2); x < ATOM_DISPLAY_WIDTH / 2 + 1; ++x) {
                        for (int y = -(ATOM_DISPLAY_WIDTH / 2); y < ATOM_DISPLAY_WIDTH / 2 + 1; ++y) {
                                set_pixel(surface, atom_x + x, atom_y + y, (RGB_color){.r = r, .g = g, .b = b});
                        }
                }
        }
}

void add_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y) { 

        // fill a list of points in which new atoms are going to go
        PointNode *head = malloc(sizeof(PointNode));
        head->data = (Point) {.x = mouse_x, .y = mouse_y};
        int n_new_atoms = 1;
        for (int x = -CLICK_PLACE_WIDTH; x < CLICK_PLACE_WIDTH; x += CLICK_PLACE_GAP) {
                for (int y = -CLICK_PLACE_WIDTH; y < CLICK_PLACE_WIDTH; y += CLICK_PLACE_GAP) {
                        Point atom_location = { .x = mouse_x + x, .y = mouse_y + y};
                        if (sqrt(x*x + y*y) < CLICK_PLACE_WIDTH) {
                                prepend_to_pointnode_list(&head, &atom_location);
                                n_new_atoms++;
                        }
                }
        }

        // https://stackoverflow.com/a/63861885
        // extend the array to contain these new atoms
        Atom *reallocated_pointer = realloc(*atoms_pointer, (*n_atoms_pointer + n_new_atoms) * sizeof(Atom));
        if (reallocated_pointer == NULL) {
                goto quit;
        }
        *atoms_pointer = reallocated_pointer;
        memset(*atoms_pointer + *n_atoms_pointer, 0, n_new_atoms * sizeof(Atom));

        // create each new atom and add it to the array
        int current_atom = *n_atoms_pointer;
        PointNode *temp = head;
        while (temp != NULL) {
                (*atoms_pointer)[current_atom] = (Atom) {
                        .mass = 1,
                        .position = (Point) { .x = temp->data.x, .y = temp->data.y },
                        .velocity = (Point) {0},
                        .acceleration = (Point) {0}
                };
                temp = temp->next;
                current_atom++;
        }
        *n_atoms_pointer += n_new_atoms;
quit:
        free_pointnode_list(head);
}

int main(void) {
        Atom *atoms = init_atoms(N_ATOMS, SCREEN_X, SCREEN_Y);
        SDL_Window *window = SDL_CreateWindow(
                "Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                SCREEN_X, SCREEN_Y, 0
        );
        SDL_Surface *surface = SDL_GetWindowSurface(window);
        SDL_Event event;
        int mouse_x, mouse_y;
        bool quit = false, gravitation_to_mouse = false;
        while (!quit) {

                SDL_GetMouseState(&mouse_x, &mouse_y);

                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                quit = true;
                        } else if (event.button.button == SDL_BUTTON_LEFT) {
                                add_atoms_upon_click(&atoms, &N_ATOMS, mouse_x, mouse_y);
                        } else if (event.button.button == SDL_BUTTON_RIGHT) {
                                // create a toggle because pressing and letting go each send events
                                gravitation_to_mouse = !gravitation_to_mouse; 
                        }
                }
        
                step_simulation(&atoms, &N_ATOMS);
                if (gravitation_to_mouse) {
                        apply_gravity_to_mouse(atoms, N_ATOMS, mouse_x, mouse_y, MOUSE_MASS);
                }
                iterate_kinematics(atoms, N_ATOMS);
                clear_screen(surface);
                display_atoms(surface, atoms, DISPLAY_COLOR);
                SDL_UpdateWindowSurface(window);
        }
        return 0;
}
