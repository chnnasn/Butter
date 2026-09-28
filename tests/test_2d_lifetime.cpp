// What destroying an object is allowed to disturb.
//
// Destroying a body used to mean `contact_cache_.clear(); constraints_.clear();
// broadphase_dirty_ = true;`, plus erasing the body from the middle of the
// array and renumbering every slot after it. So removing one object cost every
// other body in the world its contacts, its manifolds and its accumulated
// impulses -- a settled stack would re-solve from zero because an unrelated box
// a hundred metres away was deleted -- and the whole broad-phase index was
// rebuilt and re-paired to remove one set of proxies.
//
// The requirement is the opposite, and this file states it as four things that
// can be measured:
//
//   * A handle taken before a destroy does not go on naming whatever takes the
//     slot over. Slots are reused, so the index alone is not an identity.
//   * Creating and destroying bodies does not grow the world.
//   * The broad phase loses exactly the destroyed body's own candidate pairs.
//   * Only the contacts, the joints and the neighbours of the destroyed body
//     change. A settled stack's trajectory is the same whether an unrelated
//     body was deleted or left sitting where it was.
//
// Everything is read off the bodies' own transforms and the world's own counts,
// and every comparison is between two runs of the same scene over the same
// frames, so no tolerance is being spent to make a number look better.
#include <butter/physics2d/butter2d.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <ios>
#include <vector>
using namespace butter::physics2d;
static int checks = 0;
static void check(bool value, const char *message) {
    ++checks;
    if (!value)
        throw std::runtime_error(message);
}
namespace {
constexpr float kDt = 1.0f / 60.0f;
// The scene everything else is measured against: a settled ten-box stack on a
// floor, plus a loose body parked far enough away that it shares no candidate
// pair with anything. Every body is woken each frame, so the stack is genuinely
// being solved and not merely asleep.
struct Scene {
    std::unique_ptr<World> world;
    std::vector<Body *> stack;
    Body *loose{};
};
Scene scene(bool with_loose, int settle_frames = 300) {
    Scene out{};
    out.world = std::make_unique<World>();
    World &w = *out.world;
    w.create_body().static_body().at(0, -0.05f).box(40.0f, 0.05f).friction(0.6f).build();
    for (int i = 0; i < 10; ++i)
        out.stack.push_back(&w.create_body()
                                 .dynamic()
                                 .at(0, 0.5f + 1.01f * float(i))
                                 .box(0.5f, 0.5f)
                                 .friction(0.6f)
                                 .restitution(0)
                                 .build());
    if (with_loose)
        out.loose = &w.create_body().dynamic().at(60.0f, 10.0f).box(0.2f, 0.2f).build();
    for (int frame = 0; frame < settle_frames; ++frame) {
        for (Body *b : out.stack)
            b->wake();
        if (out.loose)
            out.loose->wake();
        w.step(kDt);
    }
    return out;
}
// Position, angle and both velocities of every body in `bodies`, flattened.
std::vector<float> snapshot(const std::vector<Body *> &bodies) {
    std::vector<float> out;
    out.reserve(bodies.size() * 5);
    for (const Body *b : bodies) {
        out.push_back(b->transform.position.x);
        out.push_back(b->transform.position.y);
        out.push_back(b->transform.angle);
        out.push_back(b->velocity.x);
        out.push_back(b->velocity.y);
    }
    return out;
}
float max_delta(const std::vector<float> &a, const std::vector<float> &b) {
    float worst = 0;
    for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
        worst = std::max(worst, std::abs(a[i] - b[i]));
    return worst;
}
} // namespace
int main() try {
    std::cout.setf(std::ios_base::unitbuf); // So a crash still shows how far the checks got.
    {
        // A handle is an identity, not an index. The slot is reused, so the
        // thing that makes a stale handle detectable is the generation.
        World w;
        auto &doomed = w.create_body().dynamic().at(0, 5).box(0.5f, 0.5f).build();
        const World::BodyId stale = w.id(doomed);
        check(w.alive(stale), "a handle to a live body was not alive");
        check(w.body(stale) == &doomed, "a handle to a live body resolved to the wrong object");
        check(w.body(w.id(doomed)) == &doomed, "id() and body() disagree about a live body");
        w.destroy_body(doomed);
        check(!w.alive(stale), "a destroyed body's handle was still alive");
        check(w.body(stale) == nullptr, "a destroyed body's handle resolved to something");
        check(w.body_count() == 0, "destroying the only body left a body behind");
        // The slot the body occupied comes straight back, and the handle has to
        // stay dead through the reuse.
        auto &successor = w.create_body().dynamic().at(0, 7).box(0.5f, 0.5f).build();
        check(w.body(stale) == nullptr, "a reused slot revived a stale handle");
        check(!w.alive(stale), "a reused slot made a stale handle alive again");
        check(w.id(successor) != stale, "the successor of a slot got the same handle");
        check(w.body(w.id(successor)) == &successor, "the successor's own handle was not alive");
        // Destroying the successor too leaves every handle dead.
        const World::BodyId also = w.id(successor);
        w.destroy_body(successor);
        check(!w.alive(also) && w.body(also) == nullptr, "the second handle outlived its body");
    }
    {
        // Churn does not grow the world: a slot freed by a destroy is the slot
        // the next create fills, so the counts come back to exactly where they
        // started instead of the array growing and shifting.
        Scene s = scene(true);
        const std::size_t bodies = s.world->body_count();
        check(bodies == 12, "the scene did not start with the expected body count");
        std::vector<World::BodyId> handles;
        for (int round = 0; round < 2000; ++round) {
            Body &extra = s.world->create_body().dynamic().at(200.0f, 200.0f).box(0.1f, 0.1f).build();
            handles.push_back(s.world->id(extra));
            s.world->destroy_body(extra);
        }
        check(s.world->body_count() == bodies, "churn grew the body count");
        for (const World::BodyId handle : handles)
            check(!s.world->alive(handle) && s.world->body(handle) == nullptr,
                  "a churned body's handle outlived it");
        check(s.world->joint_count() == 0, "churn invented a joint");
    }
    {
        // The broad phase loses exactly the destroyed body's own candidate
        // pairs. The count of candidate pairs is the index's own answer to "what
        // can touch what", so a delete that re-files the world shows up here as
        // a count that does not come back.
        Scene s = scene(true);
        s.world->step(kDt);
        const std::size_t before = s.world->broadphase_candidate_count();
        check(before >= 9, "the settled stack produced no candidate pairs");
        // Placed overlapping the top box, so it owns exactly one candidate pair
        // and the delete has something of its own to take away.
        auto &extra = s.world->create_body().dynamic().at(0, 10.4f).box(0.5f, 0.5f).build();
        s.world->step(kDt);
        const std::size_t grown = s.world->broadphase_candidate_count();
        check(grown == before + 1, "the added body did not own exactly one candidate pair");
        s.world->destroy_body(extra);
        s.world->step(kDt);
        std::cout << "candidate pairs: settled " << before << ", +1 body " << grown
                  << ", back to " << s.world->broadphase_candidate_count() << '\n';
        check(s.world->body_count() == 12, "the deleted body was not removed from the world");
        check(s.world->broadphase_candidate_count() == before,
              "destroying a body did not give the candidate count back");
    }
    {
        // The one that matters: the rest of the world keeps its contacts and its
        // accumulated impulses. Two identical scenes, one deleted body between
        // them, and the stack has to end up in the same place either way.
        Scene control = scene(true);
        Scene twin = scene(true);
        Scene edited = scene(true);
        check(edited.loose != nullptr, "the scene has no unrelated body to delete");
        // The comparison itself has to be exact before anything is concluded
        // from it: two identical scenes that were never edited must agree bit for
        // bit, or the tolerance below would be hiding a difference that has
        // nothing to do with the delete.
        const float identical = max_delta(snapshot(control.stack), snapshot(twin.stack));
        check(identical == 0.0f, "two identical scenes did not agree bit for bit");
        edited.world->destroy_body(*edited.loose);
        edited.loose = nullptr;
        for (int frame = 0; frame < 120; ++frame) {
            for (Body *b : control.stack)
                b->wake();
            if (control.loose)
                control.loose->wake();
            control.world->step(kDt);
            for (Body *b : edited.stack)
                b->wake();
            edited.world->step(kDt);
        }
        // The stack's own motion over these frames is not the measurement: what
        // matters is that it is the *same* motion in a world where a body a
        // hundred metres away was deleted. `identical == 0` is what says the
        // difference below is the delete and nothing else.
        const float drift = max_delta(snapshot(control.stack), snapshot(edited.stack));
        std::cout << "with and without the unrelated body, the stack differs by " << drift << '\n';
        check(drift < 1.0e-6f,
              "the stack's trajectory depended on whether an unrelated body existed");
    }
    {
        // Joints go with the body, and only the joints that touch it. The two
        // hinges either side of the deleted link cannot survive it, the rest of
        // the chain has to be untouched, and the part that is still attached has
        // to still be hanging from the anchor.
        World w;
        w.create_body().static_body().at(0, -0.05f).box(40.0f, 0.05f).friction(0.6f).build();
        auto &anchor = w.create_body().static_body().at(0, 8).box(0.1f, 0.1f).build();
        std::vector<Body *> links;
        for (int i = 1; i <= 5; ++i)
            links.push_back(&w.create_body()
                                 .dynamic()
                                 .at(0, 8 - 0.5f * float(i))
                                 .box(0.04f, 0.25f)
                                 .friction(0.4f)
                                 .build());
        Body *previous = &anchor;
        for (Body *link : links) {
            w.add_hinge_joint_at(*previous, *link, {0, link->transform.position.y + 0.25f});
            previous = link;
        }
        check(w.joint_count() == 5, "the chain did not start with five joints");
        const std::size_t bodies = w.body_count();
        Body &middle = *links[2];
        w.destroy_body(middle);
        check(w.body_count() == bodies - 1, "destroying a link did not remove it");
        // The hinges on both sides of it are gone; the three at the far end are
        // not.
        check(w.joint_count() == 3, "destroying a link removed the wrong number of joints");
        for (int frame = 0; frame < 120; ++frame) {
            for (Body *link : links)
                if (link != &middle)
                    link->wake();
            w.step(kDt);
        }
        std::cout << "chain heights: " << links[0]->transform.position.y << ' '
                  << links[1]->transform.position.y << " hanging, "
                  << links[4]->transform.position.y << " detached\n";
        // The two links still joined to the anchor are still up there, and the
        // two that were cut loose have landed on the floor below.
        check(links[0]->transform.position.y > 7.0f,
              "the link still joined to the anchor fell away");
        check(links[1]->transform.position.y > 6.4f,
              "the chain below the surviving link fell away");
        check(links[4]->transform.position.y < 1.0f && links[4]->transform.position.y > 0.0f,
              "the links cut loose from the chain did not fall to the floor");
    }
    std::cout << checks << " lifetime checks passed\n";
    return 0;
} catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
}
