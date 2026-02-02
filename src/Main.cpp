#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>

#include <glm/glm.hpp>

#include <vector>

#ifdef __EMSCRIPTEN__
    #include <emscripten.h>
    #include <emscripten/html5.h>
#endif

static constexpr int WINDOW_SIZE = 800;

// SDL
SDL_Window* window     = nullptr;
SDL_Renderer* renderer = nullptr;
bool isRunning         = true;

// MLS-MPM
struct Particle
{
    glm::vec2 x;  // position
    glm::vec2 v;  // velocity
    float mass;
};

struct Cell
{
    glm::vec2 v;  // velocity
    float mass;
};

std::vector<Particle> particles;
std::vector<Cell> grid;

// SDL
void InitializeSDL();
void Shutdown();
void Render();

// MLS-MPM
void InitializeMlsMpm();
void EachSimulationStep();

int main(int argc, char* argv[])
{
    InitializeSDL();
    InitializeMlsMpm();

    auto mainLoop = []()
    {
#ifdef __EMSCRIPTEN__
        if (!isRunning)
        {
            emscripten_cancel_main_loop();
            Shutdown();
            return;
        }
#endif

        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            switch (event.type)
            {
                case SDL_QUIT:
                    isRunning = false;
                    break;

                default:
                    break;
            }
        }

        auto state = SDL_GetKeyboardState(NULL);
        if (state[SDL_SCANCODE_ESCAPE])
        {
            isRunning = false;
        }

        EachSimulationStep();
        Render();
    };

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(mainLoop, 0, 1);
#else
    while (isRunning)
    {
        mainLoop();
    }
    Shutdown();
#endif

    return 0;
}

void InitializeSDL()
{
    int sdlResult = SDL_Init(SDL_INIT_VIDEO);
    if (sdlResult != 0)
    {
        SDL_Log("Failed to initialize SDL: %s", SDL_GetError());
        exit(1);
    }

    window = SDL_CreateWindow("MLS-MPM Sample",
                              SDL_WINDOWPOS_UNDEFINED,
                              SDL_WINDOWPOS_UNDEFINED,
                              WINDOW_SIZE,
                              WINDOW_SIZE,
                              0);
    if (!window)
    {
        SDL_Log("Failed to create window: %s", SDL_GetError());
        exit(1);
    }

    renderer = SDL_CreateRenderer(window, -1, 0);
    if (!renderer)
    {
        SDL_Log("Failed to create renderer: %s", SDL_GetError());
        exit(1);
    }
}

void Shutdown()
{
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
}

void Render()
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    filledCircleRGBA(renderer, WINDOW_SIZE / 2, WINDOW_SIZE / 2, 100, 0, 0, 255, 255);
    SDL_RenderPresent(renderer);
}

void InitializeMlsMpm()
{
    // 1.  initialise your grid - fill your grid array with (grid_res * grid_res) cells.

    // 2. create a bunch of particles. set their positions somewhere in your simulation domain.
    // initialise their deformation gradients to the identity matrix, as they're in their undeformed state.

    // 3. optionally precompute state variables e.g. particle initial volume, if your model calls for it
}

void EachSimulationStep()
{
    // 1. reset our scratch-pad grid completely
    for (auto& cell : grid)
    {
        // zero out mass and velocity for this cell
    }

    // 2. particle-to-grid (P2G).
    // goal: transfers data from particles to our grid
    for (auto& p : particles)
    {
        // 2.1: calculate weights for the 3x3 neighbouring cells surrounding the particle's position
        // on the grid using an interpolation function

        // 2.2: calculate quantities like e.g. stress based on constitutive equation

        // 2.3:
        for (auto& cell : particleNeighborhood)
        {
            // scatter our particle's momentum to the grid, using the cell's interpolation weight calculated in 2.1
        }
    }

    // 3. calculate grid velocities
    for (auto& cell : grid)
    {
        // 3.1: calculate grid velocity based on momentum found in the P2G stage

        // 3.2: enforce boundary conditions
    }

    // 4. grid-to-particle (G2P).
    // goal: report our grid's findings back to our particles, and integrate their position + velocity forward
    for (auto& p : particles)
    {
        // 4.1: update particle's deformation gradient using MLS-MPM's velocity gradient estimate
        // Reference: MLS-MPM paper, Eq. 17

        // 4.2: calculate neighbouring cell weights as in step 2.1.
        // note: our particle's haven't moved on the grid at all by this point, so the weights will be identical

        // 4.3: calculate our new particle velocities
        for (auto& cell : particleNeighborhood)
        {
            // 4.3.1:
            // get this cell's weighted contribution to our particle's new velocity
        }

        // 4.4: advect particle positions by their velocity
    }
}
