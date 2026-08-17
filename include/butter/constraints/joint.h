#pragma once

#include "butter/core/body.h"

namespace butter {

class Joint {
public:
    Body* body_a{nullptr};
    Body* body_b{nullptr};
    bool enabled{true};

    virtual ~Joint() = default;

    bool is_valid() const { return body_a != nullptr && body_b != nullptr && enabled; }
    virtual void solve(float dt) = 0;
};

} // namespace butter
