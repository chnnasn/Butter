#include <butter/physics2d/butter2d.h>
#include <iostream>
#include <stdexcept>
using namespace butter::physics2d;
static void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
static Fixture &circle(World &w, Body &b, Vec2 offset = {}, bool sensor = false) {
    Fixture f;
    f.shape = Circle{0.5f};
    f.local.position = offset;
    f.trigger = sensor;
    return w.add_fixture(b, f);
}
int main() try {
    World w({{0, 0}, 8});
    auto &a = w.create_empty_body();
    a.type = BodyType::Static;
    a.inverse_mass = 0;
    a.inverse_inertia = 0;
    auto &b = w.create_empty_body();
    b.transform.position = {0.5f, 0};
    b.linear_damping = 0;
    b.angular_damping = 0;
    auto &fa = circle(w, a, {}, true);
    auto &fb = circle(w, b);
    auto &extra = circle(w, b, {4, 0});
    int enters = 0, exits = 0;
    bool lock_checked = false;
    w.on_contact = [&](Fixture &, Fixture &, bool enter) {
        enter ? ++enters : ++exits;
        if (w.locked()) {
            try {
                w.create_empty_body();
            } catch (const std::logic_error &) {
                lock_checked = true;
            }
        }
    };
    w.step();
    w.step();
    check(enters == 1 && exits == 0 && lock_checked, "contact lifecycle and step lock");
    b.sleeping = true;
    w.step();
    check(exits == 0, "sleep must preserve trigger contacts");
    w.contact_filter = [](const Fixture &, const Fixture &) { return false; };
    w.step();
    check(exits == 1, "custom filter exit");
    w.contact_filter = {};
    w.step();
    check(enters == 2, "custom filter enter");
    fb.collision_mask = 0;
    w.step();
    check(exits == 2, "fixture mask must filter triggers");
    fb.collision_mask = 0xffffffffu;
    w.step();
    check(enters == 3, "fixture mask restored");
    w.destroy_fixture(fb);
    check(exits == 3 && b.fixtures.size() == 1 && b.fixtures.front().get() == &extra,
          "targeted fixture destruction");
    w.destroy_body(b);
    check(w.body_count() == 1, "body destruction");
    (void)fa;

    World motion({{0, -10}, 8});
    auto &k = motion.create_body().kinematic().velocity(2, 0).build();
    motion.step(0.5f);
    check(std::abs(k.transform.position.x - 1) < 1e-5f && k.transform.position.y == 0,
          "kinematic integration");
    auto &d = motion.create_empty_body();
    d.linear_damping = 0;
    d.force = {4, 0};
    d.inverse_mass = 0.5f;
    motion.step(0.5f);
    check(std::abs(d.velocity.x - 1) < 1e-5f && d.force == Vec2{},
          "force integration and clearing");

    World joints({{0, 0}, 8});
    auto &ja = joints.create_body().static_body().at(0, 0).build();
    auto &jb = joints.create_body().dynamic().at(4, 0).build();
    auto &joint = joints.add_distance_joint(ja, jb, 1);
    joint.anchor_a = {1, 0};
    joint.anchor_b = {-1, 0};
    joint.collide_connected = false;
    joints.step();
    check(std::abs(jb.transform.position.x - 3) < 1e-4f, "anchored distance joint converges");
    joints.destroy_joint(joint);
    joints.destroy_body(ja);
    joints.step();
    check(joints.body_count() == 1, "joint/body cleanup");
    auto &jc = joints.create_body().dynamic().at(6, 0).build();
    joints.add_spring_joint(jb, jc, 1, 20);
    joints.destroy_body(jc);
    joints.step();

    check(!ray_shape({-2, 0.9f}, {1, 0}, 1.2f, Shape{Circle{1}}, {}),
          "circle ray must not hit AABB corner");
    auto hit = ray_shape({-2, 0}, {1, 0}, 5, Shape{Circle{1}}, {});
    check(hit && std::abs(hit->distance - 1) < 1e-5f && hit->normal.x < -0.99f, "circle exact ray");
    check(!ray_shape({0, 0}, {1, 0}, 5, Shape{Circle{1}}, {}), "ignore rays starting inside");
    auto box_hit = ray_shape({-3, 0}, {1, 0}, 5, Shape{Box{{1, 1}}}, {{0, 0}, 0.785398163f});
    check(box_hit && std::abs(box_hit->distance - (3 - std::sqrt(2.0f))) < 1e-4f,
          "rotated box exact ray");
    Contact contact;
    check(test(Shape{Circle{0.1f}}, {}, Shape{Box{{2, 2}}}, {}, contact) && contact.penetration > 2,
          "circle contained in box");
    World compound({{0, -9.8f}, 6});
    auto &fixed = compound.create_empty_body();
    fixed.type = BodyType::Static;
    fixed.inverse_mass = 0;
    fixed.inverse_inertia = 0;
    auto &moving = compound.create_empty_body();
    moving.linear_damping = 0;
    for (Body *body : {&fixed, &moving}) {
        Fixture box;
        box.shape = Box{{5, 5}};
        compound.add_fixture(*body, box);
        Fixture disk;
        disk.shape = Circle{4};
        compound.add_fixture(*body, disk);
    }
    int active = 0;
    compound.on_contact = [&](Fixture &, Fixture &, bool enter) { active += enter ? 1 : -1; };
    compound.step();
    compound.step();
    check(active > 0, "deep compound contact must persist across bounded projection steps");
    World queries({{0, 0}, 8});
    auto &queryBody = queries.create_empty_body();
    circle(queries, queryBody, {5, 0});
    auto offsetHit = queries.raycast({2, 0}, {1, 0}, 10);
    check(offsetHit && std::abs(offsetHit->distance - 2.5f) < 1e-5f,
          "world query respects fixture offset");
    std::cout << "2D engine integration passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
}
