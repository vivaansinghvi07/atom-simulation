#include <SDL2/SDL.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "include/point.h"
#include "include/atom.h"
#include "include/gravity.h"

// dimensions for the SDL window
#define SCREEN_X 1600
#define SCREEN_Y 900

// effective mass of the mouse when right click is pressed
#define MOUSE_MASS 2000

// color mode of coloring the atoms
#define DISPLAY_COLOR COLOR_RANDOM

// controls how atoms are placed to the screen
#define CLICK_PLACE_FUNC add_rotating_atoms_upon_click
#define CLICK_PLACE_WIDTH 75
#define CLICK_PLACE_GAP 15

// controls how atom count is displayed
#define TEXT_BLOCK_WIDTH 10
#define TEXT_OFFSET 20

// detecting collisions or not
#define COLLISION_DETECTION_ON false

int N_ATOMS = 1000;
int SIMULATION_STEPS = 0;

typedef struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
} RGB_Color;

typedef struct {
        bool ready;
        bool on;
} KeySwitch;

enum ColorMode {
        COLOR_NONE, 
        COLOR_RANDOM, 
        COLOR_VELOCITY
};

// change a single pixel to a given color on the sdl surface
void set_pixel(SDL_Surface *surface, int x, int y, RGB_Color color);

// make an sdl surface entirely black
void clear_screen(SDL_Surface *surface);

// display a given array of atoms
void display_atoms(SDL_Surface *surface, Atom *atoms, enum ColorMode color_mode);

// display the number of atoms to the screen and the fps to the screen 
void display_atom_count(SDL_Surface *surface, int n_atoms);
void display_fps(SDL_Surface *surface, int fps);

// add atoms when the mouse is clicked
void add_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y);

// adds atoms to the screen but initializes them with some rotation
void add_rotating_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y);

// display a barnes hut tree around a group of atoms 
void display_barnes_hut_tree(SDL_Surface *surface, QuadTreeNode *root, int depth);

// step the simulation
void step_simulation(Atom **atoms_pointer, int *n_atoms_pointer);

// toggles a key switch, changing its state to true only if it is ready to change 
void toggle_switch(KeySwitch *key_switch);

// main loop
int main(void);

// https://github.com/MikeShah/SDL2_Tutorials/blob/main/8_ModifyingSurface/main.cpp
void set_pixel(SDL_Surface *surface, int x, int y, RGB_Color color) {
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
                int radius = atoms[i].mass;
                for (int x = -radius; x < radius + 1; ++x) {
                        for (int y = -radius; y < radius + 1; ++y) {
                                if (sqrt(x*x + y*y) < radius) {
                                        set_pixel(surface, atom_x + x, atom_y + y, (RGB_Color){.r = r, .g = g, .b = b});
                                }
                        }
                }
        }
}

/* 
 * 0  1  2
 * 3  4  6
 * 6  7  8
 * 9  10 11
 * 12 13 14
 *
 * Returns an integer where each bit in its binary representation corresponds to a place on this grid
 * 0b111001001001001  would correspond to bits 0, 1, 2, 6, 8, 11, 14 (from left ro right)
  */
int _char_to_disp_map(char c) {
        switch (c) {
                case '0': return 0b111101101101111;
                case '1': return 0b001001001001001; 
                case '2': return 0b111001111100111;
                case '3': return 0b111001111001111;
                case '4': return 0b101101111001001;
                case '5': return 0b111100111001111;
                case '6': return 0b111100111101111;
                case '7': return 0b111001001001001;
                case '8': return 0b111101111101111;
                case '9': return 0b111101111001111;
                default:  return 0b000000000000000;
        }
}

// x and y are at the top left corner of the 3x3 grid
void _draw_block_at(SDL_Surface *surface, int x, int y) {
        for (int i = 0; i < TEXT_BLOCK_WIDTH; ++i) {
                for (int j = 0; j < TEXT_BLOCK_WIDTH; ++j) {
                        set_pixel(surface, x+j, y+i, (RGB_Color) {255, 255, 255});
                }
        }
}

// x and y are at the top left corner of the 3x3 grid
void _draw_digit(SDL_Surface *surface, char digit, int x, int y) {
        int digit_map = _char_to_disp_map(digit);
        for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 3; ++j) {
                        if (digit_map & 0b100000000000000) {
                                _draw_block_at(surface, x + j * TEXT_BLOCK_WIDTH, y + i * TEXT_BLOCK_WIDTH);
                        }
                        digit_map <<= 1;
                }
        }
}

void display_atom_count(SDL_Surface *surface, int n_atoms) {

        // create a buffer storing the string
        int len = snprintf(NULL, 0, "%d", n_atoms);
        char buf[len + 1];
        sprintf(buf, "%d", n_atoms);

        // draw each number to the screen 
        for (int i = 0; i < len; ++i) {
                _draw_digit(surface, buf[i], TEXT_OFFSET + i * 4 * TEXT_BLOCK_WIDTH, TEXT_OFFSET);
        }
}

void display_fps(SDL_Surface *surface, int fps) {

        // create a buffer storing the string
        int len = snprintf(NULL, 0, "%d", fps);
        char buf[len + 1];
        sprintf(buf, "%d", fps);

        // go through each number in reverse order to print from the back
        for (int i = len - 1; i >= 0; --i) {
                _draw_digit(surface, buf[i], SCREEN_X - TEXT_OFFSET - (len - i) * 4 * TEXT_BLOCK_WIDTH, TEXT_OFFSET);
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
                        .mass = MIN_MASS + rand() % (MAX_MASS - MIN_MASS + 1),
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

void add_rotating_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y) { 

        // add the atoms normally
        int old_n_atoms = *n_atoms_pointer;
        add_atoms_upon_click(atoms_pointer, n_atoms_pointer, mouse_x, mouse_y);

        // go through each other atom and determine the resulting velocity
        // i tried to calculate this but it wasn't working so this is what we are gonna run with
        double acceleration_gravity = GRAVITATIONAL_CONSTANT * M_PI / (CLICK_PLACE_GAP * CLICK_PLACE_GAP * 400);
        Point mouse_pos = (Point) {.x = mouse_x, .y = mouse_y};
        for (int i = old_n_atoms; i < *n_atoms_pointer; ++i) {

                // calculate the resulting velocity
                Atom *atom = *atoms_pointer + i;
                Point distance = subtract_a_minus_b(&mouse_pos, &atom->position);
                double velocity = sqrt(acceleration_gravity * pow(abs_point(&distance), 2)); 
                
                // apply the velocity to be tangent 
                atom->velocity = (Point) {.x = velocity * -fast_sin_atan2(&distance),
                                          .y = velocity * fast_cos_atan2(&distance)};
        }
}

RGB_Color _get_node_color(int depth) {
        switch (depth % 8) {
                case 0:  return (RGB_Color) {255,   0,   0};
                case 1:  return (RGB_Color) {255, 127,   0};
                case 2:  return (RGB_Color) {255, 255,   0};
                case 3:  return (RGB_Color) {  0, 255,   0};
                case 4:  return (RGB_Color) {  0, 255, 255};
                case 5:  return (RGB_Color) {  0, 127, 255};
                case 6:  return (RGB_Color) {127,   0, 255};
                case 7:  return (RGB_Color) {255,   0, 255};
                default: return (RGB_Color) {  0,   0,   0};  // this will never happen
        }
}

void draw_vertical_line(SDL_Surface *surface, RGB_Color color, int x, int y_min, int y_max) {
        for (int y = y_min; y <= y_max; ++y) {
                set_pixel(surface, x, y, color);
        }
}

void draw_horizontal_line(SDL_Surface *surface, RGB_Color color, int y, int x_min, int x_max) {
        for (int x = x_min; x <= x_max; ++x) {
                set_pixel(surface, x, y, color);
        }
}

void display_barnes_hut_tree(SDL_Surface *surface, QuadTreeNode *root, int depth) {
        
        RGB_Color color = _get_node_color(depth);
        int min_x = root->bounds.min_values.x,
            max_x = root->bounds.max_values.x,
            min_y = root->bounds.min_values.y,
            max_y = root->bounds.max_values.y;
        draw_vertical_line(surface, color, min_x, min_y, max_y);
        draw_vertical_line(surface, color, max_x, min_y, max_y);
        draw_horizontal_line(surface, color, min_y, min_x, max_x);
        draw_horizontal_line(surface, color, max_y, min_x, max_x);

        if (root->top_left) {
                display_barnes_hut_tree(surface, root->top_left, depth + 1);
                display_barnes_hut_tree(surface, root->top_right, depth + 1);
                display_barnes_hut_tree(surface, root->bottom_left, depth + 1);
                display_barnes_hut_tree(surface, root->bottom_right, depth + 1);
        }
}

void step_simulation(Atom **atoms_pointer, int *n_atoms_pointer) {
        SIMULATION_STEPS++;
        GRAVITATIONAL_FUNCTION(*atoms_pointer, *n_atoms_pointer);
        if (COLLISION_DETECTION_ON) {
                apply_collision_detection_naive(*atoms_pointer, *n_atoms_pointer);
        }
        if (SIMULATION_STEPS % 14 == 0) {
                remove_faraway_atoms(atoms_pointer, n_atoms_pointer);
        }
}

void toggle_switch(KeySwitch *key_switch) {
        if (key_switch->ready) {
                key_switch->on = !key_switch->on;
                key_switch->ready = false;
        }
}

int main(void) {

        // initialize things used throughout the program
        Atom *atoms = init_atoms(N_ATOMS, SCREEN_X, SCREEN_Y);
        SDL_Window *window = SDL_CreateWindow(
                "Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                SCREEN_X, SCREEN_Y, 0
        );
        SDL_Surface *surface = SDL_GetWindowSurface(window);

        // initialize "helper variables" that control what is displayed
        SDL_Event event;
        int mouse_x, mouse_y, frames = 0;
        bool quit = false, 
             gravitation_to_mouse = false;
        KeySwitch showing_n_atoms = { .ready = true, .on = false },
                  showing_barnes_hut_tree = { .ready = true, .on = false },
                  showing_fps = { .ready = true, .on = false },
                  reset_atoms = { .ready = true, .on = false };
        double fps_start = 0;
        while (!quit) {

                SDL_GetMouseState(&mouse_x, &mouse_y);

                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                quit = true;
                        } else if (event.button.button == SDL_BUTTON_LEFT && event.button.state) {
                                CLICK_PLACE_FUNC(&atoms, &N_ATOMS, mouse_x, mouse_y);
                        } else if (event.button.button == SDL_BUTTON_RIGHT) {
                                gravitation_to_mouse = !gravitation_to_mouse; 
                        } else if (event.type == SDL_KEYUP) {
                                switch (event.key.keysym.sym) {
                                        case SDLK_n: showing_n_atoms.ready = true; break;
                                        case SDLK_d: showing_barnes_hut_tree.ready = true; break;
                                        case SDLK_f: showing_fps.ready = true; break;
                                        case SDLK_r: reset_atoms.ready = true; break;
                                }
                        } else if (event.type == SDL_KEYDOWN) {
                                switch (event.key.keysym.sym) {
                                        case SDLK_n: toggle_switch(&showing_n_atoms); break;
                                        case SDLK_d: toggle_switch(&showing_barnes_hut_tree); break;
                                        case SDLK_f: 
                                                if (showing_fps.ready) {
                                                        fps_start = (double) clock(); 
                                                        frames = 0; 
                                                }
                                                toggle_switch(&showing_fps);
                                                break;
                                        case SDLK_r:
                                                if (reset_atoms.ready) {
                                                        N_ATOMS = 0;
                                                        free(atoms);
                                                        atoms = malloc(sizeof(Atom) * 0);
                                                }
                                                toggle_switch(&reset_atoms);
                                                break;
                                }
                        } 
                }
        
                step_simulation(&atoms, &N_ATOMS);
                if (gravitation_to_mouse) {
                        apply_gravity_to_mouse(atoms, N_ATOMS, mouse_x, mouse_y, MOUSE_MASS);
                }
                iterate_kinematics(atoms, N_ATOMS);
                clear_screen(surface);
                if (showing_barnes_hut_tree.on) {
                        QuadTreeNode *root = build_barnes_hut_tree(atoms, N_ATOMS);
                        display_barnes_hut_tree(surface, root, 0);
                        free_quad_tree(root);
                }
                display_atoms(surface, atoms, DISPLAY_COLOR);
                if (showing_n_atoms.on) {
                        display_atom_count(surface, N_ATOMS);
                }
                if (showing_fps.on) {
                        display_fps(surface, frames++ * CLOCKS_PER_SEC / ((double) clock() - fps_start));
                }
                SDL_UpdateWindowSurface(window);
        }
        return 0;
}
