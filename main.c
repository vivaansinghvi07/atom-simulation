#include <SDL2/SDL.h>
#include <math.h>
#include <stdlib.h>
#include "atom.c"

#define SCREEN_X 800
#define SCREEN_Y 800
#define N_ATOMS 2000
#define ATOM_WIDTH 5

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

// https://github.com/MikeShah/SDL2_Tutorials/blob/main/8_ModifyingSurface/main.cpp
void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color) {
        if (x < 0 || x >= SCREEN_X || y < 0 || y >= SCREEN_Y) {
                return;
        }
        uint8_t *pixels = (uint8_t *)surface->pixels;
        int pixel_address = y * surface->pitch + x * surface->format->BytesPerPixel;
        pixels[pixel_address + 0] = color.g;
        pixels[pixel_address + 1] = color.b;
        pixels[pixel_address + 2] = color.r;
}

void clear_screen(SDL_Surface *surface) {
        uint8_t *pixels = (uint8_t *)surface->pixels; 
        memset(pixels, 0, SCREEN_X * SCREEN_Y * 4);  // stored as rgba
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
                                r = g = (int) ((abs_point(&atoms[i].velocity) + min_velocity) / (max_velocity - min_velocity) * 255);
                                break;
                }

                // apply coloring to each pixel for the atom
                for (int x = -(ATOM_WIDTH / 2); x < ATOM_WIDTH / 2 + 1; ++x) {
                        for (int y = -(ATOM_WIDTH / 2); y < ATOM_WIDTH / 2 + 1; ++y) {
                                set_pixel(surface, atom_x + x, atom_y + y, (RGB_color){.r = r, .g = g, .b = b});
                        }
                }
        }
}

int main() {
        Atom *atoms = init_atoms(N_ATOMS, SCREEN_X, SCREEN_Y);
        SDL_Window *window = SDL_CreateWindow(
                "Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                SCREEN_X, SCREEN_Y, 0
        );
        SDL_Surface *surface = SDL_GetWindowSurface(window);
        SDL_Event event;
        int quit = 0;
        while (!quit) {

                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                quit = 1;
                        }
                }

                step_simulation(atoms, N_ATOMS);
                clear_screen(surface);
                display_atoms(surface, atoms, COLOR_RANDOM);
                SDL_UpdateWindowSurface(window);
        }
        return 0;
}
