# Particle Simulation

This is my first project in C, and simulates particles in a 2D space with many adjustible parameters and an interactive display. Note: throughout this file and the project, "particles" are referred to as "atoms".

<video src="./recording.mov" width="320" height="240" controls></video>

## Features

- Interactive display of atoms in the simulation 
- Easy customizability of how the program works, including several functions for applying gravity
- Removal of atoms that have "escaped" the simulation 
- Optional collision detection for more realism

## Approximating Gravity

The majority of this project was exploring different methods of simulating gravity. First, gravity was applied to the center of mass of the body, which was both the most efficient and the most inaccurate. I then implemented a naive approach which checks each atom against each other atom, which is the slowest due to its quadratic time complexity.

I then started looking at approximations for gravity. I came up with my own method, which involved the following: create one grid (called the "outer grid") that surrounds the entire body of atoms. Each square of this grid is relatively large, and stores the center of mass and total mass of all the atoms in that square. Another "inner grid" is also made, which does the same thing, but with smaller (and therefore more) inner squares. Then, a final grid is made, but instead of storing center of mass and total mass, it stores the particles contained in each square of the grid as linked lists. To calculate gravity for a particle, squares in the outer grid that are far enough away are used, then squares in the inner grid that are closer but still far enough to avoid direct particle-to-particle comparison are used, and finally, for the closest regions, the linked lists are traversed to guarantee a certain accuracy percent. This algorithm ended up being significantly faster than the naive approach, but still struggled handling atom counts greater than 10,000.

I then tried to implement what is widely regarded as one of the most efficient existing algorithms for calculating gravity, the [Barnes-Hut algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation#). The Wikipedia page contains a great explanation of the process, so I won't cover it in detail here. Additionally, [this repository](https://github.com/ntta/barnes-hut-algorithm) was a great reference for a basic idea of how to implement the algorithm optimally. This ended up being significantly faster than the previous two approaches, and this project also includes a feature that allows the user to visualize the Barnes-Hut tree.

I also implemented a function that follows the (a/x\*\*n - b/x\*\*m) model for attraction/repulsion between particles. However, this model ends up creating an extremely high repulsive force at close distances that eventually causes all atoms to bounce off each other at extreme speeds.

## Usage

To use this project, simply compile the `main.c` file. Note: all source files are imported from directly, mainly because I was having linker issues when creating seperate header files. The project has one dependency, which is the SDL library.

### Configuration

The simulation is controlled by many parameters (like a gravitational constant) that can be adjusted to change the behavior of the simulation. In each file, there are parameters present that are `#define`'d; they are as follows:

<details>
<summary> In `main.c` </summary>

- `SCREEN_Y` and `SCREEN_X`: These simply control the dimensions of the SDL window that displays the simulation.
- `MOUSE_MASS`: This controls how much effective mass the mouse, which is used when applying gravity relative to the mouse.
- `ATOM_DISPLAY_WIDTH`: This is the side length of atoms that are displayed to the screen. They are displayed as boxes of pixels (this may be subject to change).
- `DISPLAY_COLOR`: This is the method of coloring atoms on the screen. Available options are `COLOR_NONE` (plain white), `COLOR_VELOCITY` (relative to atom speed), and `COLOR_RANDOM`.
- `CLICK_PLACE_WIDTH`: This is the radius of the circle in which atoms are placed.
- `CLICK_PLACE_FUNC`: This is the function that is applied when placing atoms on the window. Available options are `add_rotating_atoms_upon_click` and `add_atoms_upon_click`.
- `CLICK_PLACE_GAP`: This controls the gap between atoms when placed in the area. A higher gap makes atoms more sparse.
- `TEXT_BLOCK_WIDTH`: This controls how big each "block" is when printing text (text is displayed in a 3-by-5 grid). 
- `TEXT_OFFSET`: This is how far text is from the borders of the screen.
- `COLLISION_DETECTION_ON`: Controls if collision detection is applied or not.

</details>

<details>
<summary> In `include/gravity.c` </summary>

- `OUTER_GRID_WIDTH` and `INNER_GRID_WIDTH`: These are how many squares build up each side of the outer grid and inner grid respectively in my approach to optimizing simulating gravity.
- `OUTER_GRID_LEN` and `INNER_GRID_LEN`: These are the previous two values squared to store the length of each array used to store the grids.
- `GRID_WIDTH_THRESHOLD`: For each grid, this is the minimum number of squares you can be away from the atom to use approximations in the grid to guarantee a certain value of accuracy. This value is calculatied in the `gravity_approx_opt.py` file.
- `GRAVITATIONAL_CONSTANT`: This is simply the constant used in the equation for gravity. A higher value of this constant increases the strength of gravity.
- `GRAVITATIONAL_DISTANCE_GUARD`: In the equation for gravity, this value is added to `distance^2` in order to prevent the acceleration from being too high for particles that are too close, allowing a smoother-looking simulation. This value can be set to 0 to make gravity more accurate.
- `GRAVITATIONAL_FUNCTION`: This simply sets the function used to simulate gravity, and can be set to one of the following values: `apply_gravity_to_center` (to the center of mass), `apply_gravity_each_point_naive` (calculating gravity for each point by considering every other point directly), `apply_gravity_approx` (my algorithm for approximating gravity), and `apply_gravity_barnes_hut` (using the Barnes-Hut algorithm to approximate gravity).
- `BARNES_HUT_THRESHOLD`: The threshold that controls how accurate a Barnes-Hut simulation is, representing θ (the quotient `width of rectangle / distance of center of mass`) in the Wikipedia page.

</details>

<details>
<summary> In `include/atom.c` </summary>

- `MIN_STD_DEV_FOR_REMOVAL`: This is the minimum number of standard deviations an atom has to be away from the center of mass of all the atoms to consider removing it from the simulation, as that atom is pretty much gone.
- `REPULSION_FUNC_A` and `REPULSION_FUNC_B`: This is the value for `a` and `b` respectively in the repulsion function.
- `COLLISION_ATOM_WIDTH`: This is the effective diameter of an atom for collision detection. Two atoms with positions closer than this value are in a collision (if collision detection is turned on).

</details>

### Controls

As aforementioned, the program offers a lof of control over the display. Users can:
- Left click to add atoms according to the `CLICK_PLACE_FUNC`.
- Right click to add gravity towards the mouse.
- Press `n` to toggle displaying the atom count in the top-left corner.
- Press `f` to toggle displaying the FPS count in the top-right corner.
- Press `d` to display the Barnes-Hut tree on top of the atoms (currently unoptimized).
- Press `r` to reset the atoms and remove all of them from the simulation.
