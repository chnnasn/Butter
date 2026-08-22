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
    // `is_enter` is kept for backwards compatibility with the original
    // trigger API.  A newly-overlapping pair has is_enter=true and an ended
    // pair has is_exit=true.  Persistent overlaps are deliberately quiet:
    // trigger callbacks are enter/exit notifications, not a per-frame stream.
    // `is_stay` remains available for code that builds/forwards richer event
    // records, but World does not emit stay notifications by default.
    bool is_enter{true};
    bool is_stay{false};
    bool is_exit{false};

    bool entering() const { return is_enter; }
    bool staying() const { return is_stay; }
    bool exiting() const { return is_exit || (!is_enter && !is_stay); }
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
