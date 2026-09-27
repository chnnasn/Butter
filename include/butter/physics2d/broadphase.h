#pragma once

// Dynamic AABB tree for the 2D broad phase.
//
// Replaces the uniform grid that re-inserted every body (including static
// geometry) into cells on every step and sorted the resulting cell entries.
// A dynamic tree keeps one proxy per body and only re-inserts a proxy when its
// tight box leaves the fat box it was filed under, so resting and static bodies
// cost nothing per step. The same structure answers ray casts and AABB queries
// in O(log n) instead of scanning every body.
//
// The tree is a pure acceleration index: traversal uses fat boxes and callers
// always confirm a candidate with their own exact bounds, so results match a
// brute force scan exactly.

#include "butter/physics2d/shapes.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace butter::physics2d {

inline AABB merged(const AABB &a, const AABB &b) {
    return {{std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y)},
            {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y)}};
}

inline bool contains(const AABB &outer, const AABB &inner) {
    return outer.min.x <= inner.min.x && outer.min.y <= inner.min.y &&
           outer.max.x >= inner.max.x && outer.max.y >= inner.max.y;
}

inline bool valid(const AABB &box) {
    return box.min.x <= box.max.x && box.min.y <= box.max.y;
}

inline AABB fatten(const AABB &box, float margin) {
    return {box.min - Vec2{margin, margin}, box.max + Vec2{margin, margin}};
}

inline float perimeter(const AABB &box) {
    return 2.0f * ((box.max.x - box.min.x) + (box.max.y - box.min.y));
}

class DynamicTree {
  public:
    static constexpr std::uint32_t null_node = 0;

    DynamicTree() { reset(); }
    void reset() {
        nodes_.assign(1, Node{});
        nodes_[0].height = -1;
        root_ = null_node;
        free_ = null_node;
        count_ = 0;
    }
    std::size_t size() const { return count_; }
    std::uint32_t root() const { return root_; }
    const AABB &fat_box(std::uint32_t id) const { return nodes_[id].box; }
    const AABB &tight_box(std::uint32_t id) const { return nodes_[id].tight; }
    std::uint32_t payload(std::uint32_t id) const { return nodes_[id].payload; }
    int height() const { return root_ == null_node ? 0 : nodes_[root_].height; }

    std::uint32_t create_proxy(const AABB &tight, std::uint32_t payload, float margin) {
        const std::uint32_t id = allocate();
        nodes_[id].tight = tight;
        nodes_[id].box = fatten(tight, margin);
        nodes_[id].payload = payload;
        nodes_[id].height = 0;
        insert_leaf(id);
        ++count_;
        return id;
    }

    void destroy_proxy(std::uint32_t id) {
        remove_leaf(id);
        free_node(id);
        --count_;
    }

    // Refresh a proxy. The fat box is only rebuilt when the padded tight box
    // escapes it, which is what makes resting bodies free. The displacement
    // pads the new fat box so a steadily moving body is not re-filed every
    // step. Invariant: `box` always contains `tight` inflated by `margin`, so
    // a query can filter leaves with the same predicate the caller uses.
    bool move_proxy(std::uint32_t id, const AABB &tight, Vec2 displacement, float margin) {
        Node &node = nodes_[id];
        node.tight = tight;
        if (contains(node.box, fatten(tight, margin)))
            return false;
        remove_leaf(id);
        const float padding = margin + std::max(std::abs(displacement.x), std::abs(displacement.y));
        node.box = fatten(tight, padding);
        insert_leaf(id);
        return true;
    }

    // Visit every proxy whose bounds inflated by `margin` overlap `box`. This
    // is deliberately the same predicate callers apply themselves, so the tree
    // never hides a pair the caller would have accepted.
    template <class Visitor> void query(const AABB &box, float margin, Visitor &&visit) const {
        if (root_ == null_node)
            return;
        auto &stack = query_stack_;
        stack.clear();
        stack.push_back(root_);
        while (!stack.empty()) {
            const std::uint32_t index = stack.back();
            stack.pop_back();
            if (!nodes_[index].box.overlaps(box))
                continue;
            if (nodes_[index].leaf()) {
                if (fatten(nodes_[index].tight, margin).overlaps(box))
                    visit(nodes_[index].payload);
                continue;
            }
            stack.push_back(nodes_[index].child1);
            stack.push_back(nodes_[index].child2);
        }
    }

    // Visit the proxies whose fat box the ray crosses, nearest first, passing
    // the box entry distance. `visit` returns false to stop the traversal,
    // which lets a caller prune everything behind the closest exact hit.
    template <class Visitor>
    void raycast(Vec2 origin, Vec2 direction, float max_distance, Visitor &&visit) const {
        if (root_ == null_node)
            return;
        auto &stack = ray_stack_;
        stack.clear();
        stack.push_back(root_);
        while (!stack.empty()) {
            const std::uint32_t index = stack.back();
            stack.pop_back();
            const float distance = ray_box(origin, direction, nodes_[index].box);
            if (distance > max_distance)
                continue;
            if (nodes_[index].leaf()) {
                if (!visit(nodes_[index].payload, distance))
                    return;
                continue;
            }
            const std::uint32_t child1 = nodes_[index].child1, child2 = nodes_[index].child2;
            // Push the farther child first so the nearer one is popped first.
            const float d1 = ray_box(origin, direction, nodes_[child1].box);
            const float d2 = ray_box(origin, direction, nodes_[child2].box);
            if (d1 > d2) {
                stack.push_back(child1);
                stack.push_back(child2);
            } else {
                stack.push_back(child2);
                stack.push_back(child1);
            }
        }
    }

    // Entry distance of a ray into a box, or +inf when it misses.
    static float ray_box(Vec2 origin, Vec2 direction, const AABB &box) {
        float lower = 0, upper = std::numeric_limits<float>::max();
        for (int axis = 0; axis < 2; ++axis) {
            const float o = origin[axis], d = direction[axis];
            const float lo = box.min[axis], hi = box.max[axis];
            if (std::abs(d) < 1.0e-8f) {
                if (o < lo || o > hi)
                    return std::numeric_limits<float>::max();
                continue;
            }
            float t1 = (lo - o) / d, t2 = (hi - o) / d;
            if (t1 > t2)
                std::swap(t1, t2);
            lower = std::max(lower, t1);
            upper = std::min(upper, t2);
            if (lower > upper)
                return std::numeric_limits<float>::max();
        }
        return std::max(0.0f, lower);
    }

    // Structural sanity used by tests and diagnostics.
    bool validate() const {
        if (root_ == null_node)
            return true;
        if (nodes_[root_].parent != null_node)
            return false;
        return validate_node(root_);
    }
    std::size_t node_count() const { return nodes_.size(); }

  private:
    struct Node {
        AABB box{}, tight{};
        std::uint32_t parent{null_node}, child1{null_node}, child2{null_node};
        std::uint32_t payload{null_node};
        int height{-1};
        bool leaf() const { return child1 == null_node; }
    };

    bool validate_node(std::uint32_t index) const {
        if (nodes_[index].leaf())
            return nodes_[index].height == 0;
        const std::uint32_t child1 = nodes_[index].child1, child2 = nodes_[index].child2;
        if (nodes_[child1].parent != index || nodes_[child2].parent != index)
            return false;
        if (nodes_[index].height != 1 + std::max(nodes_[child1].height, nodes_[child2].height))
            return false;
        if (!contains(nodes_[index].box, nodes_[child1].box) ||
            !contains(nodes_[index].box, nodes_[child2].box))
            return false;
        return validate_node(child1) && validate_node(child2);
    }

    std::uint32_t allocate() {
        if (free_ == null_node) {
            nodes_.emplace_back();
            return std::uint32_t(nodes_.size() - 1);
        }
        const std::uint32_t id = free_;
        free_ = nodes_[id].parent;
        nodes_[id] = Node{};
        return id;
    }

    void free_node(std::uint32_t id) {
        nodes_[id] = Node{};
        nodes_[id].height = -1;
        nodes_[id].parent = free_;
        free_ = id;
    }

    void insert_leaf(std::uint32_t leaf) {
        if (root_ == null_node) {
            root_ = leaf;
            nodes_[leaf].parent = null_node;
            return;
        }
        const AABB leaf_box = nodes_[leaf].box;
        // Walk down towards the cheapest sibling using the surface-area
        // heuristic, exactly like Box2D's b2DynamicTree::InsertLeaf.
        std::uint32_t index = root_;
        while (!nodes_[index].leaf()) {
            const std::uint32_t child1 = nodes_[index].child1, child2 = nodes_[index].child2;
            const float area = perimeter(nodes_[index].box);
            const float combined = perimeter(merged(nodes_[index].box, leaf_box));
            const float inheritance = 2.0f * (combined - area);
            const float cost = 2.0f * combined;
            const float cost1 = 2.0f * (perimeter(merged(nodes_[child1].box, leaf_box)) -
                                        perimeter(nodes_[child1].box)) +
                                inheritance;
            const float cost2 = 2.0f * (perimeter(merged(nodes_[child2].box, leaf_box)) -
                                        perimeter(nodes_[child2].box)) +
                                inheritance;
            if (cost < cost1 && cost < cost2)
                break;
            index = cost1 < cost2 ? child1 : child2;
        }
        const std::uint32_t sibling = index;
        const std::uint32_t old_parent = nodes_[sibling].parent;
        const std::uint32_t new_parent = allocate();
        nodes_[new_parent].parent = old_parent;
        nodes_[new_parent].box = merged(leaf_box, nodes_[sibling].box);
        nodes_[new_parent].height = nodes_[sibling].height + 1;
        nodes_[new_parent].payload = null_node;
        nodes_[new_parent].child1 = sibling;
        nodes_[new_parent].child2 = leaf;
        nodes_[sibling].parent = new_parent;
        nodes_[leaf].parent = new_parent;
        if (old_parent != null_node) {
            if (nodes_[old_parent].child1 == sibling)
                nodes_[old_parent].child1 = new_parent;
            else
                nodes_[old_parent].child2 = new_parent;
        } else
            root_ = new_parent;
        refit_upward(nodes_[leaf].parent);
    }

    void remove_leaf(std::uint32_t leaf) {
        if (leaf == root_) {
            root_ = null_node;
            return;
        }
        const std::uint32_t parent = nodes_[leaf].parent;
        const std::uint32_t grand = nodes_[parent].parent;
        const std::uint32_t sibling =
            nodes_[parent].child1 == leaf ? nodes_[parent].child2 : nodes_[parent].child1;
        if (grand != null_node) {
            if (nodes_[grand].child1 == parent)
                nodes_[grand].child1 = sibling;
            else
                nodes_[grand].child2 = sibling;
            nodes_[sibling].parent = grand;
            free_node(parent);
            refit_upward(grand);
        } else {
            root_ = sibling;
            nodes_[sibling].parent = null_node;
            free_node(parent);
        }
        nodes_[leaf].parent = null_node;
    }

    void refit_upward(std::uint32_t index) {
        while (index != null_node) {
            index = balance(index);
            const std::uint32_t child1 = nodes_[index].child1, child2 = nodes_[index].child2;
            nodes_[index].box = merged(nodes_[child1].box, nodes_[child2].box);
            nodes_[index].height = 1 + std::max(nodes_[child1].height, nodes_[child2].height);
            index = nodes_[index].parent;
        }
    }

    // Rotate a subtree when one side is more than one level deeper.
    std::uint32_t balance(std::uint32_t a) {
        if (nodes_[a].leaf() || nodes_[a].height < 2)
            return a;
        const std::uint32_t b = nodes_[a].child1, c = nodes_[a].child2;
        const int difference = nodes_[c].height - nodes_[b].height;
        if (difference > 1) {
            const std::uint32_t f = nodes_[c].child1, g = nodes_[c].child2;
            const std::uint32_t parent = nodes_[a].parent;
            nodes_[c].child1 = a;
            nodes_[c].parent = parent;
            nodes_[a].parent = c;
            if (parent != null_node) {
                if (nodes_[parent].child1 == a)
                    nodes_[parent].child1 = c;
                else
                    nodes_[parent].child2 = c;
            } else
                root_ = c;
            if (nodes_[f].height > nodes_[g].height) {
                nodes_[c].child2 = f;
                nodes_[a].child2 = g;
                nodes_[g].parent = a;
                nodes_[a].box = merged(nodes_[b].box, nodes_[g].box);
                nodes_[c].box = merged(nodes_[a].box, nodes_[f].box);
                nodes_[a].height = 1 + std::max(nodes_[b].height, nodes_[g].height);
                nodes_[c].height = 1 + std::max(nodes_[a].height, nodes_[f].height);
            } else {
                nodes_[c].child2 = g;
                nodes_[a].child2 = f;
                nodes_[f].parent = a;
                nodes_[a].box = merged(nodes_[b].box, nodes_[f].box);
                nodes_[c].box = merged(nodes_[a].box, nodes_[g].box);
                nodes_[a].height = 1 + std::max(nodes_[b].height, nodes_[f].height);
                nodes_[c].height = 1 + std::max(nodes_[a].height, nodes_[g].height);
            }
            return c;
        }
        if (difference < -1) {
            const std::uint32_t f = nodes_[b].child1, g = nodes_[b].child2;
            const std::uint32_t parent = nodes_[a].parent;
            nodes_[b].child1 = a;
            nodes_[b].parent = parent;
            nodes_[a].parent = b;
            if (parent != null_node) {
                if (nodes_[parent].child1 == a)
                    nodes_[parent].child1 = b;
                else
                    nodes_[parent].child2 = b;
            } else
                root_ = b;
            if (nodes_[f].height > nodes_[g].height) {
                nodes_[b].child2 = f;
                nodes_[a].child1 = g;
                nodes_[g].parent = a;
                nodes_[a].box = merged(nodes_[c].box, nodes_[g].box);
                nodes_[b].box = merged(nodes_[a].box, nodes_[f].box);
                nodes_[a].height = 1 + std::max(nodes_[c].height, nodes_[g].height);
                nodes_[b].height = 1 + std::max(nodes_[a].height, nodes_[f].height);
            } else {
                nodes_[b].child2 = g;
                nodes_[a].child1 = f;
                nodes_[f].parent = a;
                nodes_[a].box = merged(nodes_[c].box, nodes_[f].box);
                nodes_[b].box = merged(nodes_[a].box, nodes_[g].box);
                nodes_[a].height = 1 + std::max(nodes_[c].height, nodes_[f].height);
                nodes_[b].height = 1 + std::max(nodes_[a].height, nodes_[g].height);
            }
            return b;
        }
        return a;
    }

    std::vector<Node> nodes_;
    std::uint32_t root_{null_node}, free_{null_node};
    std::size_t count_{};
    // Query scratch. Kept on the tree so repeated queries do not allocate.
    mutable std::vector<std::uint32_t> query_stack_, ray_stack_;
};

} // namespace butter::physics2d
