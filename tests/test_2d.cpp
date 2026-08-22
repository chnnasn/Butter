#include <butter/physics2d/butter2d.h>
#include <cassert>
#include <cmath>

using namespace butter::physics2d;

int main() {
    World world({{0, -10}, 8, 1.0f / 60.0f});
    world.create_body().static_body().at(0, -1).box(10, 1).friction(0.8f).build();
    auto& ball = world.create_body().dynamic().at(0, 5).circle(0.5f).restitution(0).build();
    for (int i = 0; i < 240; ++i) world.step();
    assert(ball.transform.position.y > -0.51f && ball.transform.position.y < 0.2f);
    assert(std::abs(ball.velocity.y) < 0.5f);

    Contact contact;
    assert(test(Shape{Circle{1}}, {{0, 0}, 0}, Shape{Circle{1}}, {{1.5f, 0}, 0}, contact));
    assert(contact.penetration > 0.4f);
    assert(!test(Shape{Circle{1}}, {{0, 0}, 0}, Shape{Circle{1}}, {{3, 0}, 0}, contact));
    assert(test(Shape{Circle{0.5f}}, {{0, 0}, 0}, Shape{Box{{1, 1}}}, {{0.9f, 0}, 0}, contact));
    assert(test(Shape{Capsule{0.25f, 0.75f}}, {{0, 0}, 0}, Shape{Circle{0.5f}}, {{0, 1.1f}, 0}, contact));
    Shape mesh = Mesh{{Polygon{{{-1, 0}, {1, 0}, {0, 1}}}}};
    assert(test(Shape{Circle{0.25f}}, {{0, 0.1f}, 0}, mesh, {{0, 0}, 0}, contact));

    bool entered = false, exited = false;
    World trigger_world({{0, 0}, 2, 1.0f / 60.0f});
    trigger_world.on_trigger = [&](Body&, Body&, bool enter) { entered = entered || enter; exited = exited || !enter; };
    trigger_world.create_body().static_body().at(0, 0).circle(2).trigger().build();
    auto& trigger_ball = trigger_world.create_body().dynamic().at(0, 0).circle(0.5f).build();
    trigger_world.step();
    assert(entered);
    trigger_ball.transform.position = {5, 0};
    trigger_world.step();
    assert(exited);

    World query_world({{0, 0}, 4, 1.0f / 60.0f});
    query_world.create_body().static_body().at(0, 0).box(1, 1).build();
    assert(query_world.query_aabb({{-2, -2}, {2, 2}}).size() == 1);
    const auto hit = query_world.raycast({-5, 0}, {1, 0}, 20);
    assert(hit && hit->distance > 3.9f && hit->distance < 4.1f);
    World joint_world({{0, 0}, 4, 1.0f / 60.0f});
    auto& ja = joint_world.create_body().dynamic().at(0, 0).circle(0.2f).build();
    auto& jb = joint_world.create_body().dynamic().at(3, 0).circle(0.2f).build();
    joint_world.add_distance_joint(ja, jb, 1.0f);
    joint_world.step();
    assert((jb.transform.position - ja.transform.position).length() < 2.5f);
    return 0;
}
