// Same scene definitions, two native engines. Timing never includes quality
// measurement, construction or forced wake-up. CSV is suitable for budget
// sweeps; it deliberately does not label unequal errors as a speedup.
#include <butter/physics2d/butter2d.h>
#include <box2d/box2d.h>
#include "task_pool.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

static thread_local bool count_new = false;
static thread_local std::size_t new_calls = 0, new_bytes = 0;
void *operator new(std::size_t n) {
    void *p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    if (count_new) { ++new_calls; new_bytes += n; }
    return p;
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }

namespace bp = butter::physics2d;
constexpr float dt = 1.0f / 60;
struct BodyDef {
    float x{}, y{}, hx{.5f}, hy{.5f}, mass{1}, vx{}, vy{};
    bool fixed{}, bullet{}, circle{};
    bool kinematic{};
};
struct JointDef { int a{}, b{}; float x{}, y{}; };
struct Scene {
    std::string name;
    std::vector<BodyDef> bodies;
    std::vector<JointDef> joints;
    float gravity{-9.81f};
    bool floor{true};
    int top{-1};
    float ideal_top{};
    int probe{-1};
};
Scene scene(const std::string &name, int count) {
    Scene s; s.name = name;
    if (name == "queries" || name == "free" || name == "fast") {
        s.gravity = 0; s.floor = false;
    }
    if (s.floor) s.bodies.push_back({0, -.5f, float(std::max(1000, count)), .5f, 0, 0, 0, true});
    if (name == "stack" || name == "mass_ratio" || name == "obstacles") {
        const int n = name == "mass_ratio" ? 11 : 20;
        for (int i = 0; i < n; ++i)
            s.bodies.push_back({0, .5f + float(i) * 1.01f, .5f, .5f,
                               name == "mass_ratio" && i == n-1 ? 50.0f : 1.0f});
        s.top = int(s.bodies.size()) - 1; s.ideal_top = float(n) - .5f;
        if (name == "obstacles")
            for (int i = 0; i < count; ++i)
                s.bodies.push_back({50 + float(i % 100) * 3, 10 + float(i / 100) * 3,
                                   .5f, .5f, 0, 0, 0, true});
    } else if (name == "chain") {
        s.floor = false;
        s.bodies.clear();
        s.bodies.push_back({0, 10, .1f, .1f, 0, 0, 0, true});
        for (int i = 1; i <= 24; ++i) {
            s.bodies.push_back({0, 10 - i * .5f, .04f, .25f});
            s.joints.push_back({i-1, i, 0, 10 - i * .5f + .25f});
        }
    } else if (name == "fast") {
        s.bodies.push_back({0, 0, .05f, float(count) * 2, 0, 0, 0, true});
        for (int i = 0; i < count; ++i)
            s.bodies.push_back({-5, float(i) * 2, .1f, .1f, 1, 100, 0, false, true});
    } else if (name == "circles") {
        for (int i = 0; i < count; ++i) {
            BodyDef d; d.x = (float(i) - count * .5f) * 1.5f; d.y = .5f; d.circle = true;
            s.bodies.push_back(d);
        }
    } else if (name == "queries" || name == "free") {
        for (int i = 0; i < count; ++i)
            s.bodies.push_back({5 + float(i) * 3, 0, .5f, .5f, 1,
                               name == "free" ? 1.0f : 0.0f, 0, name == "queries"});
    } else throw std::invalid_argument("unknown scenario");
    if (name != "queries") {
        s.probe = int(s.bodies.size());
        BodyDef probe; probe.y = 10000; probe.vx = 1; probe.kinematic = true;
        s.bodies.push_back(probe);
    }
    return s;
}
struct State { float x, y, angle, vx, vy, spin; bool awake; };
struct Butter {
    std::unique_ptr<TaskPool> tasks;
    bp::World world;
    std::vector<bp::Body *> bodies;
    std::vector<bp::HingeJoint *> joints;
    std::string mode;
    static bp::World::Config config(const Scene &s, int iterations, int substeps, const std::string &mode, TaskPool *pool) {
        bp::World::Config c;
        c.gravity = {0, s.gravity}; c.solver_iterations = iterations; c.substeps = substeps;
        if (mode == "baseline") {
            c.cache_solver_islands = c.buffered_velocities = false;
            c.ccd.tight_fixture_bounds = false;
            c.detailed_profiling = true;
        } else if (mode == "soft") c.contact_hertz = 20;
        else if (mode == "simd") c.simd_contacts = true;
        else if (mode == "anchors") c.translated_box_contacts = true;
        if (pool) c.parallel_for = [pool](std::size_t n, const std::function<void(std::size_t)> &job) { pool->run(n,job); };
        return c;
    }
    Butter(const Scene &s, int it, int sub, const std::string &m)
        : tasks(m == "parallel" ? std::make_unique<TaskPool>(3) : nullptr),
          world(config(s, it, sub, m, tasks.get())), mode(m) {
        for (const auto &d : s.bodies) {
            auto b = world.create_body().at(d.x, d.y).box(d.hx, d.hy).mass(d.mass)
                        .velocity(d.vx, d.vy).friction(.6f).restitution(0).bullet(d.bullet);
            if (d.fixed) b.static_body();
            if (d.kinematic) b.kinematic();
            if (d.circle) b.circle(d.hx);
            auto &body = b.build();
            body.linear_damping = body.angular_damping = 0;
            bodies.push_back(&body);
        }
        for (const auto &j : s.joints) {
            auto &joint = world.add_hinge_joint_at(*bodies[j.a], *bodies[j.b], {j.x, j.y});
            joint.collide_connected = false; joints.push_back(&joint);
        }
    }
    void wake() { for (auto *b : bodies) if (b->type == bp::BodyType::Dynamic) b->wake(); }
    void step() { world.step(dt); }
    State state(std::size_t i) const {
        const auto &b = *bodies[i];
        return {b.transform.position.x,b.transform.position.y,b.transform.angle,
                b.velocity.x,b.velocity.y,b.angular_velocity,!b.sleeping};
    }
    double gap() const {
        double result = 0;
        for (auto *j : joints) {
            auto a = j->a->transform.position + bp::rotate(j->anchor_a, bp::Rot(j->a->transform.angle));
            auto b = j->b->transform.position + bp::rotate(j->anchor_b, bp::Rot(j->b->transform.angle));
            result = std::max(result, double((a-b).length()));
        }
        return result;
    }
    double query() {
        double sum = 0;
        auto run = [&](const auto &q) {
            for (int i = 0; i < 100; ++i) {
                auto hit = q.raycast({0,0},{1,0},100000);
                if (hit) sum += hit->distance;
            }
        };
        if (mode == "baseline") run(world); else world.query_batch(run);
        return sum;
    }
    std::size_t limits() const {
        const auto &c = world.ccd_statistics();
        return c.budget_exhaustions + c.non_convergences + c.zero_time_repeats;
    }
};
struct Box2D {
    b2World world;
    std::vector<b2Body *> bodies;
    std::vector<b2Joint *> joints;
    int iterations, substeps;
    Box2D(const Scene &s, int it, int sub, const std::string &) : world({0,s.gravity}), iterations(it), substeps(sub) {
        for (const auto &d : s.bodies) {
            b2BodyDef bd; bd.type = d.fixed ? b2_staticBody : b2_dynamicBody;
            if (d.kinematic) bd.type = b2_kinematicBody;
            bd.position = {d.x,d.y}; bd.linearVelocity = {d.vx,d.vy}; bd.bullet = d.bullet;
            auto *b = world.CreateBody(&bd);
            b2PolygonShape shape; shape.SetAsBox(d.hx,d.hy);
            b2CircleShape circle; circle.m_radius = d.hx;
            b2FixtureDef fd; fd.shape = &shape; fd.density = d.mass / (4*d.hx*d.hy);
            if (d.circle) { fd.shape = &circle; fd.density = d.mass / (b2_pi*d.hx*d.hx); }
            fd.friction = .6f; fd.restitution = 0;
            b->CreateFixture(&fd); bodies.push_back(b);
        }
        for (const auto &j : s.joints) {
            b2RevoluteJointDef jd; jd.Initialize(bodies[j.a], bodies[j.b], {j.x,j.y});
            jd.collideConnected = false; joints.push_back(world.CreateJoint(&jd));
        }
    }
    void wake() { for (auto *b : bodies) if (b->GetType() == b2_dynamicBody) b->SetAwake(true); }
    void step() { for (int i=0;i<substeps;++i) world.Step(dt/substeps,iterations,3); }
    State state(std::size_t i) const {
        const auto *b = bodies[i]; auto p=b->GetPosition(),v=b->GetLinearVelocity();
        return {p.x,p.y,b->GetAngle(),v.x,v.y,b->GetAngularVelocity(),b->IsAwake()};
    }
    double gap() const {
        double result=0;
        for (auto *j:joints) result=std::max(result,double((j->GetAnchorA()-j->GetAnchorB()).Length()));
        return result;
    }
    double query() {
        struct Ray : b2RayCastCallback {
            float distance{};
            float ReportFixture(b2Fixture *, const b2Vec2 &, const b2Vec2 &, float fraction) override {
                distance=fraction*100000; return fraction;
            }
        };
        double sum=0;
        for(int i=0;i<100;++i) { Ray ray; world.RayCast(&ray,{0,0},{100000,0}); sum+=ray.distance; }
        return sum;
    }
    std::size_t limits() const { return 0; } // Not exposed by Box2D; printed as unavailable below.
};

template<class Engine> void run(const Scene &s, const char *engine_name, int frames, int it, int sub,
                                 bool awake, const std::string &mode, int round) {
    Engine engine(s,it,sub,mode);
    std::vector<double> times; times.reserve(frames);
    std::size_t active=0, allocations=0, bytes=0, limits=0;
    double penetration=0,gap=0,query_sum=0,residual=0;
    for(int frame=-60;frame<frames;++frame) {
        if(awake) engine.wake();
        if(frame>=0)
            for(std::size_t i=0;i<s.bodies.size();++i)
                active += !s.bodies[i].fixed && !s.bodies[i].kinematic && engine.state(i).awake;
        new_calls=new_bytes=0; count_new=frame>=0;
        const auto start=std::chrono::steady_clock::now();
        if(s.name=="queries") query_sum=engine.query(); else engine.step();
        const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        count_new=false;
        if(frame>=0) { times.push_back(ms); allocations+=new_calls; bytes+=new_bytes; limits+=engine.limits(); }
        if(s.probe >= 0) {
            const double expected = double(frame + 61) * dt;
            if(std::abs(engine.state(s.probe).x - expected) > 1e-4 + 1e-5 * expected)
                throw std::runtime_error("independent motion probe did not advance full time");
        } else if(std::abs(query_sum - 450) > .01) throw std::runtime_error("ray query checksum failed");
        for(std::size_t i=0;i<s.bodies.size();++i) if(!s.bodies[i].fixed && !s.bodies[i].kinematic) {
            auto b=engine.state(i);
            if(!std::isfinite(b.x+b.y+b.angle+b.vx+b.vy+b.spin))
                throw std::runtime_error("non-finite simulation state");
            const auto &d=s.bodies[i];
            const double ex=d.circle ? d.hx : std::abs(std::cos(b.angle))*d.hx+std::abs(std::sin(b.angle))*d.hy;
            const double ey=d.circle ? d.hx : std::abs(std::sin(b.angle))*d.hx+std::abs(std::cos(b.angle))*d.hy;
            if(s.floor) penetration=std::max(penetration,ey-b.y);
            if(s.name=="fast") penetration=std::max(penetration,double(b.x)+ex+.05);
            residual=std::max(residual,std::hypot(double(b.vx),double(b.vy)));
        }
        gap=std::max(gap,engine.gap());
    }
    double mean=0; for(double t:times) mean+=t; mean/=times.size();
    std::sort(times.begin(),times.end());
    auto quantile=[&](double q){return times[std::size_t(std::ceil(q*times.size()))-1];};
    const double top_error=s.top<0 ? 0 : std::abs(engine.state(s.top).y-s.ideal_top);
    std::printf("%s,%s,%s,%d,%d,%d,%d,%d,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%zu,%zu,%zu,",
        engine_name,mode.c_str(),s.name.c_str(),round,int(s.bodies.size()),frames,it,sub,awake?"awake":"natural",
        mean,quantile(.5),quantile(.95),quantile(.99),times.back(),penetration,gap,top_error,active,allocations,bytes);
    if constexpr(std::is_same_v<Engine,Butter>) std::printf("%zu",limits); else std::printf("unavailable");
    const bool quality = penetration <= (s.name == "fast" ? .001 : .025) && gap <= .05 && top_error <= .5;
    std::printf(",%.9f,%.9f,%.9f,%d\n",query_sum,residual,s.name=="queries"?0:double(frames)*dt,int(quality));
    std::fflush(stdout);
}
int main(int argc,char **argv) try {
    const std::string selected=argc>1?argv[1]:"all";
    const int frames=argc>2?std::stoi(argv[2]):300,count=argc>3?std::stoi(argv[3]):1000;
    const int iterations=argc>4?std::stoi(argv[4]):8,substeps=argc>5?std::stoi(argv[5]):1;
    const bool awake=argc>6 && std::string(argv[6])=="awake";
    const std::string mode=argc>7?argv[7]:"optimized";
    if(frames<1||frames>20000||count<1||count>10000||iterations<1||substeps<1||substeps>16 ||
       (mode!="optimized"&&mode!="baseline"&&mode!="soft"&&mode!="simd"&&mode!="parallel"&&mode!="anchors"))
        throw std::invalid_argument("invalid arguments");
    std::fprintf(stderr,"compiler=%s config=%s box2d=2.4.1 source_sha256=%s butter_sha256=%s dt=1/60 warmup=60\n",
        BUTTER_COMPILER,BUTTER_BUILD_CONFIG,BUTTER_BOX2D_HASH,BUTTER_SOURCE_HASH);
    std::fprintf(stderr,"Butter position iterations=velocity iterations; Box2D position iterations=3. Equal settings are not equal quality.\n");
    std::fprintf(stderr,"cpp_new excludes aligned allocation, malloc and Box2D block allocator; CCD limits unavailable for Box2D.\n");
    std::fprintf(stderr,"Box2D polygon skin is its default 0.01m; outline height error includes skin gaps. Quality gates: floor<=25mm, fast wall<=1mm, joint<=50mm, stack height error<=0.5m.\n");
    std::puts("engine,mode,scenario,round,bodies,frames,iterations,substeps,activity,mean_ms,p50_ms,p95_ms,p99_ms,max_ms,max_penetration_m,max_joint_gap_m,final_top_error_m,awake_body_steps,cpp_new_calls,cpp_new_bytes,ccd_limits,query_checksum,max_speed,measured_simulation_s,quality_pass");
    bool found = false;
    for(const std::string name:{"stack","mass_ratio","chain","fast","obstacles","queries","free","circles"}) {
        if(selected!="all"&&selected!=name) continue;
        found = true;
        const auto s=scene(name,count);
        for(int round=0;round<3;++round) {
            if(round%2==0) run<Butter>(s,"Butter",frames,iterations,substeps,awake,mode,round);
            run<Box2D>(s,"Box2D-2.4.1",frames,iterations,substeps,awake,mode,round);
            if(round%2!=0) run<Butter>(s,"Butter",frames,iterations,substeps,awake,mode,round);
        }
    }
    if (!found) throw std::invalid_argument("unknown scenario");
    return 0;
} catch(const std::exception &e) { std::fprintf(stderr,"FAIL %s\n",e.what()); return 1; }
