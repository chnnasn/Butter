#pragma once

#include <functional>
#include <utility>
#include <vector>

#include "butter/math/vec3.h"

namespace butter {

class Body;

struct CollisionEvent {
    Body& body_a;
    Body& body_b;
    math::Vec3 point;
    math::Vec3 normal;
    float impulse{0};
    float penetration{0};

    // Defined after Body is complete in body.h.
    float relative_velocity() const;
};

struct TriggerEvent {
    Body& trigger;
    Body& other;
    bool is_enter{true};
};

struct DestructionEvent {
    Body& body;
    math::Vec3 position;
};

template <typename... Args>
class Event {
public:
    using Callback = std::function<void(Args...)>;

    Event& operator+=(Callback callback) {
        callbacks_.push_back(std::move(callback));
        return *this;
    }

    // std::function does not support value comparison. Use clear() or keep a
    // handle outside this wrapper if a specific listener needs removal.
    Event& operator-=(const Callback&) { return *this; }

    void invoke(Args... args) {
        for (auto& callback : callbacks_) {
            callback(args...);
        }
    }

    void operator()(Args... args) {
        invoke(std::forward<Args>(args)...);
    }

    void clear() { callbacks_.clear(); }

    std::size_t size() const { return callbacks_.size(); }
    bool empty() const { return callbacks_.empty(); }

private:
    std::vector<Callback> callbacks_;
};

} // namespace butter
