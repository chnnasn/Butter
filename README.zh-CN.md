# Butter

> **Butter** - Make Physics Smooth Again!
> 一个现代、灵活、符合直觉的 C++20 header-only 物理引擎。

[English](README.md)

## 简介

Butter 追求像 C#/TypeScript 一样流畅的 API，同时保持 C++ 的高性能和零开销抽象。

## 快速开始

```cpp
#include <butter/butter.h>
#include <iostream>

using namespace butter;
using namespace butter::math;

int main() {
    World world;
    world.gravity = {0, -9.81f, 0};

    world.create_body()
        .static_body()
        .at(0, 0, 0)
        .box(50, 0.5f, 50)
        .friction(0.8f)
        .build();

    auto& ball = world.create_body()
        .dynamic()
        .at(0, 10, 0)
        .sphere(0.5f)
        .bounciness(0.7f);

    for (int i = 0; i < 120; ++i) {
        world.step(1.0f / 60.0f);
        std::cout << ball.position.y << "\n";
    }
}
```

## 构建

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## 测试

测试使用轻量断言，不依赖第三方测试框架。

| 测试 | 覆盖内容 |
| --- | --- |
| `test_math` | 向量、矩阵、四元数、AABB |
| `test_body` | Builder、质量、冲量、碰撞体 |
| `test_collision` | 球-球、球-盒、盒-盒、分离检测 |
| `test_world` | 重力、raycast、球落地稳定 |
| `test_features` | 接触点冲量、碰撞过滤、阻尼 |

运行全部测试：

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

## 模块

- `butter/math` - 向量、矩阵、四元数、变换、包围盒
- `butter/core` - 世界、刚体、Builder、材质、事件
- `butter/shapes` - 球、盒、胶囊、凸包与碰撞检测
- `butter/constraints` - 距离、弹簧、铰链关节
- `butter/query` - 射线检测与重叠查询

## 示例

```bash
cmake -S . -B build -DBUTTER_BUILD_EXAMPLES=ON
cmake --build build
./build/examples/hello_butter
./build/examples/falling_boxes
./build/examples/chain
./build/examples/playground
```

基于 GLFW 的爆炸炸碎箱子堆示例有独立文档：
[examples/README.zh-CN.md](examples/README.zh-CN.md)。

```bash
cmake -S . -B build -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --target exploding_crates
./build/examples/exploding_crates
```

Windows 下若 `E:/Github/glfw` 存在，CMake 会自动检测；否则显式传入
`-DBUTTER_GLFW_DIR=<你的 GLFW 源码路径>`。GLFW 仅作为外部依赖链接，不复制进 Butter。

## 路线图

### 已完成

- C++20 header-only 核心库
- 数学库、刚体、Builder、世界、事件系统
- 球、盒、胶囊碰撞；盒-盒 SAT
- 冲量求解、摩擦、位置修正
- 接触点角速度计算，翻滚物体会被摩擦减速
- 碰撞过滤、线性/角阻尼、睡眠
- 爆炸 + 动态破碎示例

### 近期

- CCD 连续碰撞检测，彻底解决高速穿模
- 更高效的 broadphase（BVH / 网格）
- 堆叠稳定性和接触求解质量提升
- 更多关节：滑轨、固定、马达

### 中期

- 凸包与网格碰撞
- 多线程求解
- 简单软体、布料、弹簧骨骼
- 车辆与角色控制器

### 远期

- 流体与粒子
- 软体 FEM
- 可视化调试器与编辑器

## 许可证

本项目使用 [MIT License](LICENSE)。
