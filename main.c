#include "stdio.h"
#include <SDL2/SDL.h>
#include <stdint.h>
#include <string.h>
#include "include/test.c"

#define SCREEN_X 400
#define SCREEN_Y 400
#define N_ATOMS 100

typedef struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
} RGB_color;

void debug_atoms(Atom *atoms);

// https://github.com/MikeShah/SDL2_Tutorials/blob/main/8_ModifyingSurface/main.cpp
void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color) {
        if (x < 0 || x >= SCREEN_X || y < 0 || y >= SCREEN_Y) {
                return;
        }
        uint8_t *pixels = (uint8_t *)surface->pixels;
        uint pixel_address = y * surface->pitch + x * surface->format->BytesPerPixel;
        pixels[pixel_address + 0] = color.r;
        pixels[pixel_address + 1] = color.g;
        pixels[pixel_address + 2] = color.b;
}

void clear_screen(SDL_Surface *surface) {
        uint8_t *pixels = (uint8_t *)surface->pixels; 
        memset(pixels, 0, SCREEN_X * SCREEN_Y * 4);  // stored as rgba
}

void init_atoms(Atom *atoms) {
        for (int i = 0; i < N_ATOMS; ++i) {
                atoms[i].mass = 1; 
                atoms[i].position = (Point) {rand() % SCREEN_X, rand() % SCREEN_Y};
        }
}

void debug_atoms(Atom *atoms) {
        for (int i = 0; i < N_ATOMS; ++i) {
                Atom atom = atoms[i];
                printf("Mass: %d\n", atom.mass);
                printf("Position: %.2g %.2g\n", atom.position.x, atom.position.y);
                printf("Velocity: %.2g %.2g\n", atom.velocity.x, atom.velocity.y);
                printf("Acceleration: %.2g %.2g\n", atom.acceleration.x, atom.acceleration.y);
        }
}

void display_atoms(SDL_Surface *surface, Atom *atoms) {
        for (int i = 0; i < N_ATOMS; ++i) { 
                int atom_x = atoms[i].position.x;
                int atom_y = atoms[i].position.y;
                for (int x = -1; x < 2; ++x) {
                        for (int y = -1; y < 2; ++y) {
                                set_pixel(surface, atom_x + x, atom_y + y, (RGB_color){255, 255, 255});
                        }
                }
        }
}

int main() {
         
        // determine what the dimensions for the simulation will be
        Atom *atoms = malloc(sizeof(Atom) * N_ATOMS);
        init_atoms(atoms);
        
        // initialize the window and run the window control sequence
        SDL_Window *window = SDL_CreateWindow(
                "Simulation", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                SCREEN_X, SCREEN_Y, 0
        );
        SDL_Surface *surface = SDL_GetWindowSurface(window);
        SDL_Event event;
        int quit = 0, running = 0;
        while (!quit) {

                while (SDL_PollEvent(&event)) {
                        if (event.type == SDL_QUIT) {
                                quit = 1;
                        } else if (event.button.button == SDL_BUTTON_LEFT) { 
                                running = !running;
                        }
                }
                if (!running) {
                        continue;
                }

                clear_screen(surface);
                display_atoms(surface, atoms);
                SDL_UpdateWindowSurface(window);
        }
        free(atoms);
        return 0;
}
