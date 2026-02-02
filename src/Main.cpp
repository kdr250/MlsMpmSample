#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>

#include <glm/glm.hpp>

#include <vector>
#include <random>

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
// Grid resolution (cells)
const int n = 80;

const float dt       = 1e-4f;
const float frame_dt = 1e-3f;
const float dx       = 1.0f / n;
const float inv_dx   = 1.0f / dx;

// Snow material properties
const auto particle_mass = 1.0f;
const auto vol           = 1.0f;   // Particle Volume
const auto hardening     = 10.0f;  // Snow hardening factor
const auto E             = 1e4f;   // Young's Modulus
const auto nu            = 0.2f;   // Poisson ratio
const bool plastic       = true;

// Initial Lamé parameters
const float mu_0     = E / (2 * (1 + nu));
const float lambda_0 = E * nu / ((1 + nu) * (1 - 2 * nu));

struct Particle
{
    // Position and velocity
    glm::vec2 x, v;
    // Deformation gradient
    glm::mat2 F;
    // Affine momentum from APIC
    glm::mat2 C;
    // Determinant of the deformation gradient (i.e. volume)
    float Jp;
    // Color
    uint32_t c;

    Particle(glm::vec2 x, uint32_t c, glm::vec2 v = glm::vec2(0)) :
        x(x), v(v), F(1), C(0), Jp(1), c(c)
    {
    }
};

std::vector<Particle> particles;
glm::vec3 grid[n + 1][n + 1];

// SDL
void InitializeSDL();
void Shutdown();
void Render();
glm::vec2 ZOToSDLPosition(const glm::vec2& zeroToOnePosition);
SDL_Color ToSDLColor(uint32_t color);

float Random();
glm::vec2 RandomVec();

// MLS-MPM
void InitializeMlsMpm();
void AddObject(const glm::vec2& center, uint32_t color);
void Advance(float dt);

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

        Advance(dt);
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
    for (auto& particle : particles)
    {
        auto screenPos = ZOToSDLPosition(particle.x);
        auto color     = ToSDLColor(particle.c);
        filledCircleRGBA(renderer, screenPos.x, screenPos.y, 1, color.r, color.g, color.b, color.a);
    }
    SDL_RenderPresent(renderer);
}

glm::vec2 ZOToSDLPosition(const glm::vec2& zeroToOnePosition)
{
    glm::vec2 result = zeroToOnePosition;
    result.y         = 1.0f - result.y;

    return result * (float)WINDOW_SIZE;
}

SDL_Color ToSDLColor(uint32_t color)
{
    return SDL_Color {(Uint8)((color >> 24) & 0xFF),
                      (Uint8)((color >> 16) & 0xFF),
                      (Uint8)((color >> 8) & 0xFF),
                      (Uint8)(color & 0xFF)};
}

float Random()
{
    static std::random_device device;
    static std::mt19937_64 generator(device());
    static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    return distribution(generator);
}

glm::vec2 RandomVec()
{
    return glm::vec2(Random(), Random());
}

void InitializeMlsMpm()
{
    AddObject(glm::vec2(0.55, 0.45), 0xFF0000FF);
    AddObject(glm::vec2(0.45, 0.65), 0x00FF00FF);
    AddObject(glm::vec2(0.55, 0.85), 0x0000FFFF);
}

void AddObject(const glm::vec2& center, uint32_t color)
{
    // Randomly sample 1000 particles in the square
    for (int i = 0; i < 1000; ++i)
    {
        particles.push_back(
            Particle {(RandomVec() * 2.0f - glm::vec2(1.0f)) * 0.08f + center, color});
    }
}

void Advance(float dt)
{
    // Reset grid
    std::memset(grid, 0, sizeof(grid));

    // P2G
    for (auto& p : particles)
    {
        // element-wise floor
        glm::vec2 base_coord = p.x * inv_dx - glm::vec2(0.5f);
        glm::ivec2 ibase_coord((int)base_coord.x, (int)base_coord.y);

        glm::vec2 fx = p.x * inv_dx - base_coord;

        // Quadratic kernels [http://mpm.graphics Eqn. 123, with x=fx, fx-1,fx-2]
        glm::vec2 w[3] = {
            glm::vec2(0.5f) * glm::sqrt(glm::vec2(1.5f) - fx),
            glm::vec2(0.75f) - glm::sqrt(fx - glm::vec2(1.0f)),
            glm::vec2(0.5f) * glm::sqrt(fx - glm::vec2(0.5f)),
        };

        // Compute current Lamé parameters [http://mpm.graphics Eqn. 86]
        auto e      = std::exp(hardening * (1.0f - p.Jp));
        auto mu     = mu_0 * e;
        auto lambda = lambda_0 * e;

        // current volume
        float J = glm::determinant(p.F);

        // TODO: Polar decomposition for fixed corotated model
        glm::mat2 r, s;

        // [http://mpm.graphics Paragraph after Eqn. 176]
        float Dinv = 4 * inv_dx * inv_dx;
        // [http://mpm.graphics Eqn. 52]
        auto PF = (2 * mu * (p.F - r) * glm::transpose(p.F) + lambda * (J - 1) * J);

        // Cauchy stress times dt and inv_dx
        auto stress = -(dt * vol) * (Dinv * PF);

        // Fused APIC momentum + MLS-MPM stress contribution
        // See http://taichi.graphics/wp-content/uploads/2019/03/mls-mpm-cpic.pdf
        // Eqn 29
        auto affine = stress + particle_mass * p.C;

        // P2G
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                auto dpos = (glm::vec2(i, j) - fx) * dx;
                // Translational momentum
                glm::vec3 mass_x_velocity(p.v * particle_mass, particle_mass);
                grid[ibase_coord.x + i][ibase_coord.y + j] +=
                    (w[i].x * w[j].y * (mass_x_velocity + glm::vec3(affine * dpos, 0)));
            }
        }
    }

    // For all grid nodes
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            auto& g = grid[i][j];
            // No need for epsilon here
            if (g[2] > 0)
            {
                // Normalize by mass
                g /= g[2];
                // Gravity
                g += dt * glm::vec3(0, -200, 0);

                // boundary thickness
                float boundary = 0.05;
                // Node coordinates
                float x = (float)i / n;
                float y = float(j) / n;

                // Sticky boundary
                if (x < boundary || x > 1 - boundary || y > 1 - boundary)
                {
                    g = glm::vec3(0);
                }
                // Separate boundary
                if (y < boundary)
                {
                    g[1] = std::max(0.0f, g[1]);
                }
            }
        }
    }

    // G2P
    for (auto& p : particles)
    {
        // element-wise floor
        glm::vec2 base_coord = (p.x * inv_dx - glm::vec2(0.5f));
        glm::ivec2 ibase_coord((int)base_coord.x, (int)base_coord.y);
        glm::vec2 fx   = p.x * inv_dx - base_coord;
        glm::vec2 w[3] = {glm::vec2(0.5) * glm::sqrt(glm::vec2(1.5) - fx),
                          glm::vec2(0.75) - glm::sqrt(fx - glm::vec2(1.0)),
                          glm::vec2(0.5) * glm::sqrt(fx - glm::vec2(0.5))};

        p.C = glm::mat2(0);
        p.v = glm::vec2(0);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                glm::vec2 dpos   = glm::vec2(i, j) - fx;
                glm::vec2 grid_v = glm::vec2(grid[ibase_coord.x + i][ibase_coord.y + j]);
                auto weight      = w[i].x * w[j].y;
                // Velocity
                p.v += weight * grid_v;
                // APIC C
                p.C += 4 * inv_dx
                       * glm::cross(glm::vec3(weight * grid_v, 0.0f), glm::vec3(dpos, 0.0f)).z;
            }
        }

        // Advection
        p.x += dt * p.v;

        // TODO: MLS-MPM F-update
    }
}
