#include <GLFW/glfw3.h>

#include <butter/butter.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

using namespace butter;
using namespace butter::math;

namespace {

struct Camera {
    float yaw = 0.0f;
    float pitch = 0.6f;
    float dist = 22.0f;
};

struct Crate {
    Body* body = nullptr;
    Vec3 half{0.5f, 0.5f, 0.5f};
};

struct Debris {
    Body* body = nullptr;
    Vec3 half{0.12f, 0.12f, 0.12f};
};

struct Flash {
    Vec3 position;
    float life = 0.0f;
};

class ExplosionScene {
public:
    ExplosionScene() { reset(); }

    World& world() { return world_; }
    const std::vector<Crate>& crates() const { return crates_; }
    const std::vector<Debris>& debris() const { return debris_; }
    const std::vector<Flash>& flashes() const { return flashes_; }
    Body* ground_body() const { return ground_; }
    bool exploded() const { return exploded_; }
    float shake() const { return shake_; }

    void reset() {
        world_ = World{};
        crates_.clear();
        debris_.clear();
        flashes_.clear();
        exploded_ = false;
        shake_ = 0.0f;
        blast_time_ = 0.0f;
        build();
    }

    void tick(float dt) {
        for (auto& flash : flashes_) {
            flash.life -= dt;
        }
        flashes_.erase(std::remove_if(flashes_.begin(), flashes_.end(),
                                      [](const Flash& flash) { return flash.life <= 0.0f; }),
                       flashes_.end());

        shake_ = std::max(0.0f, shake_ - dt * 2.6f);
    }

    void detonate() {
        if (exploded_) return;
        exploded_ = true;
        shake_ = 1.0f;

        const Vec3 center{0.0f, 2.0f, 0.0f};
        constexpr float power = 120.0f;
        std::uniform_real_distribution<float> unit(0.0f, 1.0f);

        blast_center_ = center;
        blast_time_ = 0.35f;
        flashes_.push_back({center, 0.55f});

        // Compute the explosion impulse for a position. Near bodies receive a
        // much larger impulse and therefore break into more fragments.
        auto explosion_impulse = [&](const Vec3& position) {
            Vec3 delta = position - center;
            const float distance = std::max(0.2f, delta.length());
            const Vec3 direction = delta / distance;
            const float strength = power / (1.0f + distance * distance);
            return direction * strength + Vec3{0, strength * 0.45f, 0};
        };

        constexpr float fracture_threshold = 28.0f;

        for (auto it = crates_.begin(); it != crates_.end();) {
            Body& body = *it->body;
            if (body.is_destroyed()) {
                it = crates_.erase(it);
                continue;
            }

            const Vec3 impulse = explosion_impulse(body.position);
            if (impulse.length() > fracture_threshold) {
                fracture_crate(*it, impulse);
                it = crates_.erase(it);
            } else {
                const Vec3 contact = body.position +
                                     Vec3{unit(rng_) - 0.5f, unit(rng_) - 0.5f, unit(rng_) - 0.5f} * 0.6f;
                body.apply_impulse_at_point(impulse, contact);
                ++it;
            }
        }

        for (auto& piece : debris_) {
            if (piece.body && !piece.body->is_destroyed()) {
                const Vec3 impulse = explosion_impulse(piece.body->position);
                const Vec3 contact = piece.body->position +
                                     Vec3{unit(rng_) - 0.5f, unit(rng_) - 0.5f, unit(rng_) - 0.5f} * 0.4f;
                piece.body->apply_impulse_at_point(impulse, contact);
            }
        }
    }

    // Sustained blast wave. A single impulse is only one frame of motion; this
    // keeps pushing bodies for a short time so they visibly accelerate along
    // different radial/tangential trajectories.
    void apply_blast_force(float dt) {
        if (blast_time_ <= 0.0f) return;
        blast_time_ -= dt;

        std::uniform_real_distribution<float> unit(0.0f, 1.0f);
        constexpr float radius = 9.0f;

        auto push = [&](Body& body) {
            Vec3 delta = body.position - blast_center_;
            const float distance = std::max(0.3f, delta.length());
            if (distance > radius) return;

            const Vec3 radial = delta / distance;
            const float falloff = 120.0f / (1.0f + distance * distance);

            // A random tangential component gives each object its own spin and
            // curved trajectory instead of everyone flying straight outward.
            Vec3 random_dir{unit(rng_) - 0.5f, unit(rng_) - 0.5f, unit(rng_) - 0.5f};
            random_dir = random_dir.normalized();
            Vec3 tangential = random_dir - radial * random_dir.dot(radial);
            const float tangential_length = tangential.length();
            if (tangential_length > 1.0e-3f) {
                tangential = tangential / tangential_length;
            } else {
                tangential = {0, 1, 0};
            }

            const Vec3 force = (radial * falloff + tangential * (falloff * 0.4f)) * body.mass();
            body.apply_force(force);
        };

        for (auto& crate : crates_) {
            if (crate.body && !crate.body->is_destroyed()) push(*crate.body);
        }
        for (auto& piece : debris_) {
            if (piece.body && !piece.body->is_destroyed()) push(*piece.body);
        }
    }

private:
    void build() {
        world_.gravity = {0, -9.81f, 0};
        world_.config.enable_sleeping = true;
        world_.config.sleep_threshold = 0.15f;

        // Smaller physics steps dramatically reduce tunneling for fast debris.
        world_.config.fixed_timestep = 1.0f / 240.0f;
        world_.config.velocity_iterations = 6;
        world_.config.position_iterations = 4;

        ground_ = &world_.create_body()
            .static_body()
            .at(0, 0, 0)
            .box(20, 0.5f, 20)
            .friction(0.95f)
            .build();

        constexpr int layers = 3;
        constexpr int side = 4;
        const Vec3 half{0.5f, 0.5f, 0.5f};
        std::uniform_real_distribution<float> unit(0.0f, 1.0f);

        for (int layer = 0; layer < layers; ++layer) {
            const float base_y = 5.0f + static_cast<float>(layer) * 2.0f;
            for (int x = 0; x < side; ++x) {
                for (int z = 0; z < side; ++z) {
                    const float base_x = static_cast<float>(x - side / 2) + 0.5f;
                    const float base_z = static_cast<float>(z - side / 2) + 0.5f;
                    const float px = base_x + (unit(rng_) - 0.5f) * 0.4f;
                    const float pz = base_z + (unit(rng_) - 0.5f) * 0.4f;
                    const float y = base_y + (unit(rng_) - 0.5f) * 0.5f;

                    Body& body = world_.create_body()
                        .dynamic()
                        .at(px, y, pz)
                        .mass(1.0f)
                        .box(0.5f, 0.5f, 0.5f)
                        .friction(0.85f)
                        .bounciness(0.05f)
                        .build();
                    body.linear_damping = 0.15f;
                    body.angular_damping = 0.6f;
                    body.rotation = Quat::from_euler((unit(rng_) - 0.5f) * 0.6f,
                                                     (unit(rng_) - 0.5f) * 0.6f,
                                                     (unit(rng_) - 0.5f) * 0.6f);
                    crates_.push_back({&body, half});
                }
            }
        }
    }

    void fracture_crate(const Crate& crate, const Vec3& impulse) {
        Body& original = *crate.body;
        const Vec3 position = original.position;
        const Vec3 velocity = original.velocity;
        const Vec3 angular_velocity = original.angular_velocity;
        const Quat rotation = original.rotation;
        const float mass = original.mass();

        world_.destroy(original);

        std::uniform_real_distribution<float> unit(0.0f, 1.0f);
        const int fragment_count = std::clamp(static_cast<int>(impulse.length() / 12.0f), 6, 14);
        const float fragment_mass = mass / static_cast<float>(fragment_count);

        for (int i = 0; i < fragment_count; ++i) {
            const Vec3 offset{unit(rng_) * 0.7f - 0.35f,
                              unit(rng_) * 0.7f - 0.35f,
                              unit(rng_) * 0.7f - 0.35f};
            const Vec3 half{0.06f + unit(rng_) * 0.10f,
                            0.06f + unit(rng_) * 0.10f,
                            0.06f + unit(rng_) * 0.10f};

            Body& fragment = world_.create_body()
                .dynamic()
                .at(position + offset)
                .mass(fragment_mass)
                .box(half.x, half.y, half.z)
                .friction(0.5f)
                .bounciness(0.25f)
                .build();
            fragment.linear_damping = 0.3f;
            fragment.angular_damping = 0.9f;
            fragment.rotation = rotation *
                Quat::from_euler((unit(rng_) - 0.5f) * 1.2f,
                                 (unit(rng_) - 0.5f) * 1.2f,
                                 (unit(rng_) - 0.5f) * 1.2f);
            fragment.angular_velocity = angular_velocity +
                Vec3{(unit(rng_) - 0.5f) * 12.0f,
                     (unit(rng_) - 0.5f) * 12.0f,
                     (unit(rng_) - 0.5f) * 12.0f};

            Vec3 random_dir{unit(rng_) - 0.5f, unit(rng_) - 0.5f, unit(rng_) - 0.5f};
            random_dir = random_dir.normalized();
            fragment.velocity = velocity + impulse + random_dir * (impulse.length() * 0.3f);
            debris_.push_back({&fragment, half});
        }
    }

    World world_;
    Body* ground_ = nullptr;
    std::vector<Crate> crates_;
    std::vector<Debris> debris_;
    std::vector<Flash> flashes_;
    std::mt19937 rng_{1234};
    bool exploded_ = false;
    float shake_ = 0.0f;
    float blast_time_ = 0.0f;
    Vec3 blast_center_{0.0f, 2.0f, 0.0f};
};

struct App {
    GLFWwindow* window = nullptr;
    Camera camera;
    ExplosionScene scene;
    bool paused = false;
    bool dragging = false;
    double last_x = 0.0;
    double last_y = 0.0;
    float time_scale = 1.0f;
};

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

void key_callback(GLFWwindow* window, int key, int /*scancode*/, int action, int /*mods*/) {
    App* app = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (action != GLFW_PRESS) return;

    if (key == GLFW_KEY_ESCAPE) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    } else if (key == GLFW_KEY_SPACE) {
        app->scene.detonate();
    } else if (key == GLFW_KEY_P) {
        app->paused = !app->paused;
    } else if (key == GLFW_KEY_R) {
        app->scene.reset();
    } else if (key == GLFW_KEY_UP) {
        app->camera.dist = std::max(10.0f, app->camera.dist - 1.0f);
    } else if (key == GLFW_KEY_DOWN) {
        app->camera.dist = std::min(80.0f, app->camera.dist + 1.0f);
    } else if (key == GLFW_KEY_LEFT) {
        app->time_scale = std::max(0.1f, app->time_scale * 0.8f);
    } else if (key == GLFW_KEY_RIGHT) {
        app->time_scale = std::min(3.0f, app->time_scale * 1.25f);
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int /*mods*/) {
    App* app = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        app->dragging = (action == GLFW_PRESS);
        glfwGetCursorPos(window, &app->last_x, &app->last_y);
    }
}

void cursor_position_callback(GLFWwindow* window, double x, double y) {
    App* app = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (!app->dragging) return;

    const double dx = x - app->last_x;
    const double dy = y - app->last_y;
    app->camera.yaw += static_cast<float>(dx) * 0.005f;
    app->camera.pitch += static_cast<float>(dy) * 0.005f;
    app->camera.pitch = std::clamp(app->camera.pitch, -1.45f, 1.45f);
    app->last_x = x;
    app->last_y = y;
}

void scroll_callback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    App* app = static_cast<App*>(glfwGetWindowUserPointer(window));
    if (yoffset > 0) {
        app->camera.dist = std::max(10.0f, app->camera.dist * 0.9f);
    } else {
        app->camera.dist = std::min(80.0f, app->camera.dist * 1.1f);
    }
}

void draw_cube(const Vec3& half, const Vec3& color, float alpha = 1.0f) {
    const float hx = half.x;
    const float hy = half.y;
    const float hz = half.z;

    auto shade = [&](float f) {
        return Vec3{std::min(1.0f, color.x * f),
                    std::min(1.0f, color.y * f),
                    std::min(1.0f, color.z * f)};
    };

    auto face = [&](const Vec3& tint, float ax, float ay, float az,
                    float bx, float by, float bz,
                    float cx, float cy, float cz,
                    float dx, float dy, float dz) {
        glColor4f(tint.x, tint.y, tint.z, alpha);
        glVertex3f(ax, ay, az);
        glVertex3f(bx, by, bz);
        glVertex3f(cx, cy, cz);
        glVertex3f(dx, dy, dz);
    };

    glBegin(GL_QUADS);
    face(shade(0.95f), -hx, hy, -hz, hx, hy, -hz, hx, hy, hz, -hx, hy, hz); // +Y
    face(shade(0.30f), -hx, -hy, -hz, -hx, -hy, hz, hx, -hy, hz, hx, -hy, -hz); // -Y
    face(shade(0.80f), hx, -hy, -hz, hx, hy, -hz, hx, hy, hz, hx, -hy, hz); // +X
    face(shade(0.55f), -hx, -hy, -hz, -hx, -hy, hz, -hx, hy, hz, -hx, hy, -hz); // -X
    face(shade(0.85f), -hx, -hy, hz, hx, -hy, hz, hx, hy, hz, -hx, hy, hz); // +Z
    face(shade(0.60f), -hx, -hy, -hz, -hx, hy, -hz, hx, hy, -hz, hx, -hy, -hz); // -Z
    glEnd();
}

void set_model_transform(const Vec3& position, const Quat& rotation) {
    glPushMatrix();
    glTranslatef(position.x, position.y, position.z);

    const Mat3 rot = rotation.to_mat3();
    const float matrix[16] = {
        rot.m[0][0], rot.m[1][0], rot.m[2][0], 0.0f,
        rot.m[0][1], rot.m[1][1], rot.m[2][1], 0.0f,
        rot.m[0][2], rot.m[1][2], rot.m[2][2], 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f,
    };
    glMultMatrixf(matrix);
}

void render(const ExplosionScene& scene, const Camera& camera) {
    int width = 0;
    int height = 0;
    GLFWwindow* window = glfwGetCurrentContext();
    glfwGetFramebufferSize(window, &width, &height);
    if (width == 0 || height == 0) return;

    glViewport(0, 0, width, height);
    glClearColor(0.06f, 0.06f, 0.09f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const float aspect = static_cast<float>(width) / static_cast<float>(height);
    const float fov = 45.0f * (pi / 180.0f);
    const float near_plane = 0.1f;
    const float top = near_plane * std::tan(fov * 0.5f);
    glFrustum(-top * aspect, top * aspect, -top, top, near_plane, 500.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -camera.dist);
    glRotatef(camera.pitch * (180.0f / pi), 1.0f, 0.0f, 0.0f);
    glRotatef(camera.yaw * (180.0f / pi), 0.0f, 1.0f, 0.0f);

    // Camera shake from the blast.
    const float shake = scene.shake();
    if (shake > 0.0f) {
        const double t = glfwGetTime() * 80.0;
        glTranslatef(static_cast<float>(std::sin(t)) * 0.14f * shake,
                     static_cast<float>(std::cos(t * 1.3)) * 0.14f * shake,
                     0.0f);
    }

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    // Ground.
    const Vec3 ground_pos = scene.ground_body()->position;
    glPushMatrix();
    glTranslatef(ground_pos.x, ground_pos.y, ground_pos.z);
    draw_cube({20.0f, 0.5f, 20.0f}, {0.28f, 0.28f, 0.32f});
    glPopMatrix();

    // Crates.
    for (const auto& crate : scene.crates()) {
        if (!crate.body || crate.body->is_destroyed()) continue;
        const float speed = crate.body->velocity.length();
        const Vec3 color = speed > 3.0f ? Vec3{1.0f, 0.45f, 0.10f}
                                        : Vec3{0.78f, 0.42f, 0.12f};
        set_model_transform(crate.body->position, crate.body->rotation);
        draw_cube(crate.half, color);
        glPopMatrix();
    }

    // Debris fragments.
    for (const auto& piece : scene.debris()) {
        if (!piece.body || piece.body->is_destroyed()) continue;
        const Vec3 color{0.95f, 0.55f, 0.15f};
        set_model_transform(piece.body->position, piece.body->rotation);
        draw_cube(piece.half, color);
        glPopMatrix();
    }

    // Explosion flash: hot core plus an expanding shockwave ring.
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glDepthMask(GL_FALSE);
    for (const auto& flash : scene.flashes()) {
        constexpr float max_life = 0.55f;
        const float age = max_life - flash.life;
        const float progress = std::clamp(age / max_life, 0.0f, 1.0f);
        const float fade = 1.0f - progress;

        glPointSize(10.0f + 60.0f * fade);
        glColor4f(1.0f, 0.9f, 0.6f, fade * 0.9f);
        glBegin(GL_POINTS);
        glVertex3f(flash.position.x, flash.position.y, flash.position.z);
        glEnd();

        const float radius = 1.2f + progress * 8.0f;
        glLineWidth(2.5f * fade + 0.5f);
        glColor4f(1.0f, 0.85f, 0.55f, fade * 0.8f);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 140; ++i) {
            const float a = static_cast<float>(i) / 140.0f * 2.0f * pi;
            glVertex3f(flash.position.x + std::cos(a) * radius,
                       flash.position.y,
                       flash.position.z + std::sin(a) * radius);
        }
        glEnd();

        glLineWidth(1.5f * fade + 0.5f);
        glColor4f(1.0f, 1.0f, 0.9f, fade * 0.5f);
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 140; ++i) {
            const float a = static_cast<float>(i) / 140.0f * 2.0f * pi;
            glVertex3f(flash.position.x + std::cos(a) * radius * 0.82f,
                       flash.position.y + 0.1f,
                       flash.position.z + std::sin(a) * radius * 0.82f);
        }
        glEnd();
    }
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
}

} // namespace

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1280, 720,
                                          "Butter Physics: Exploding Crate Stack",
                                          nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    App app;
    app.window = window;

    glfwSetWindowUserPointer(window, &app);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetScrollCallback(window, scroll_callback);

    auto previous_time = std::chrono::steady_clock::now();
    while (!glfwWindowShouldClose(window)) {
        const auto current_time = std::chrono::steady_clock::now();
        float frame_dt = std::chrono::duration<float>(current_time - previous_time).count();
        previous_time = current_time;
        frame_dt = std::clamp(frame_dt, 0.0f, 0.05f);

        if (!app.paused) {
            app.scene.apply_blast_force(frame_dt);
            app.scene.world().step(frame_dt * app.time_scale);
        }
        app.scene.tick(frame_dt);

        render(app.scene, app.camera);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
