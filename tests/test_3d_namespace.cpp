#include <butter/physics3d/butter3d.h>
#include <cassert>

int main() {
    butter::physics3d::World world;
    auto& body = world.create_body().dynamic().at(0, 1, 0).sphere(0.5f).build();
    assert(body.is_dynamic());
    return 0;
}
