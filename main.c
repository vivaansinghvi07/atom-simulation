#include <SDL2/SDL.h>
#include "atom.c"

#define SCREEN_X 1600
#define SCREEN_Y 800
#define N_ATOMS 100
#define ATOM_WIDTH 5

typedef struct {
        uint8_t r;
        uint8_t g;
        uint8_t b;
} RGB_color;

// https://github.com/MikeShah/SDL2_Tutorials/blob/main/8_ModifyingSurface/main.cpp
void set_pixel(SDL_Surface *surface, int x, int y, RGB_color color) {
        if (x < 0 || x >= SCREEN_X || y < 0 || y >= SCREEN_Y) {
                return;
        }
        uint8_t *pixels = (uint8_t *)surface->pixels;
        int pixel_address = y * surface->pitch + x * surface->format->BytesPerPixel;
        pixels[pixel_address + 0] = color.r;
        pixels[pixel_address + 1] = color.g;
        pixels[pixel_address + 2] = color.b;
}

void clear_screen(SDL_Surface *surface) {
        uint8_t *pixels = (uint8_t *)surface->pixels; 
        memset(pixels, 0, SCREEN_X * SCREEN_Y * 4);  // stored as rgba
}

void display_atoms(SDL_Surface *surface, Atom *atoms) {
        for (int i = 0; i < N_ATOMS; ++i) { 
                int atom_x = atoms[i].position.x;
                int atom_y = atoms[i].position.y;
                for (int x = -(ATOM_WIDTH / 2); x < ATOM_WIDTH / 2 + 1; ++x) {
                        for (int y = -(ATOM_WIDTH / 2); y < ATOM_WIDTH / 2 + 1; ++y) {
                                set_pixel(surface, atom_x + x, atom_y + y, (RGB_color){255, 255, 255});
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
                display_atoms(surface, atoms);
                SDL_UpdateWindowSurface(window);
        }
        return 0;
}
