#include <SDL2/SDL.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "include/point.c"
#include "include/atom.c"
#include "include/gravity.c"

// dimensions for the SDL window
#define SCREEN_X 1600
#define SCREEN_Y 900

// effective mass of the mouse when right click is pressed
#define MOUSE_MASS 200

// width of each atom being displayed to the screen
#define ATOM_DISPLAY_WIDTH 3

// color mode of coloring the atoms
#define DISPLAY_COLOR COLOR_VELOCITY

// controls how atoms are placed to the screen
#define CLICK_PLACE_WIDTH 50
#define CLICK_PLACE_GAP 5

// controls how atom count is displayed
#define TEXT_BLOCK_WIDTH 10
#define TEXT_OFFSET 20

int N_ATOMS = 0;
int SIMULATION_STEPS = 0;

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

// change a single pixel to a given color on the sdl surface
void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color);

// make an sdl surface entirely black
void clear_screen(SDL_Surface *surface);

// display a given array of atoms
void display_atoms(SDL_Surface *surface, Atom *atoms, enum ColorMode color_mode);

// add atoms when the mouse is clicked
void add_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y);

// step the simulation
void step_simulation(Atom **atoms_pointer, int *n_atoms_pointer);

// main loop
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

void _draw_block_at(SDL_Surface *surface, int x, int y) {
        for (int i = 0; i < TEXT_BLOCK_WIDTH; ++i) {
                for (int j = 0; j < TEXT_BLOCK_WIDTH; ++j) {
                        set_pixel(surface, x+j, y+i, (RGB_color) {255, 255, 255});
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
                int digit_map = _char_to_disp_map(buf[i]);
                for (int y = 0; y < 5; ++y) {
                        for (int x = 0; x < 3; ++x) {
                                if (digit_map & 0b100000000000000) {
                                        _draw_block_at(surface, TEXT_OFFSET + (i * 4 + x) * TEXT_BLOCK_WIDTH, 
                                                                TEXT_OFFSET + y * TEXT_BLOCK_WIDTH);
                                }
                                digit_map <<= 1;
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

void add_rotating_atoms_upon_click(Atom **atoms_pointer, int *n_atoms_pointer, int mouse_x, int mouse_y) { 

        // add the atoms normally
        int old_n_atoms = *n_atoms_pointer;
        add_atoms_upon_click(atoms_pointer, n_atoms_pointer, mouse_x, mouse_y);

        // go through each other atom and determine the resulting velocity
        double acceleration_gravity = 0.2 / (CLICK_PLACE_GAP * CLICK_PLACE_GAP);  // TODO: MODIFY THIS
        Point mouse_pos = (Point) {.x = mouse_x, .y = mouse_y};
        for (int i = old_n_atoms; i < *n_atoms_pointer; ++i) {

                // calculate the resulting velocity
                Atom *atom = *atoms_pointer + i;
                Point distance = subtract_a_minus_b(&mouse_pos, &atom->position);
                double velocity = sqrt(acceleration_gravity * abs_point(&distance));  // \sqrt{a * r} 
                
                // apply the velocity to be tangent 
                atom->velocity = (Point) {.x = velocity * -fast_sin_atan2(&distance),
                                          .y = velocity * fast_cos_atan2(&distance)};
        }

}

void step_simulation(Atom **atoms_pointer, int *n_atoms_pointer) {
        SIMULATION_STEPS++;
        GRAVITATIONAL_FUNCTION(*atoms_pointer, *n_atoms_pointer);
        if (SIMULATION_STEPS % 7 == 0) {
                apply_collision_detection_naive(*atoms_pointer, *n_atoms_pointer);
        }
        if (SIMULATION_STEPS % 14 == 0) {
                remove_faraway_atoms(atoms_pointer, n_atoms_pointer);
        }
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
        bool quit = false, gravitation_to_mouse = false, showing_n_atoms = false;
        while (!quit) {

                SDL_GetMouseState(&mouse_x, &mouse_y);

                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                quit = true;
                        } else if (event.button.button == SDL_BUTTON_LEFT && event.button.state) {
                                add_rotating_atoms_upon_click(&atoms, &N_ATOMS, mouse_x, mouse_y);
                        } else if (event.button.button == SDL_BUTTON_RIGHT) {
                                gravitation_to_mouse = !gravitation_to_mouse; 
                        } else if (event.type == SDL_KEYUP) {
                                switch (event.key.keysym.sym) {
                                        case SDLK_n: showing_n_atoms = false; break;
                                }
                        } else if (event.type == SDL_KEYDOWN) {
                                switch (event.key.keysym.sym) {
                                        case SDLK_n: showing_n_atoms = true; break;
                                }
                        } 
                }
        
                step_simulation(&atoms, &N_ATOMS);
                if (gravitation_to_mouse) {
                        apply_gravity_to_mouse(atoms, N_ATOMS, mouse_x, mouse_y, MOUSE_MASS);
                }
                iterate_kinematics(atoms, N_ATOMS);
                clear_screen(surface);
                display_atoms(surface, atoms, DISPLAY_COLOR);
                if (showing_n_atoms) {
                        display_atom_count(surface, N_ATOMS);
                }
                SDL_UpdateWindowSurface(window);
        }
        return 0;
}
