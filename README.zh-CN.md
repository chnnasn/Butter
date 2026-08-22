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
| `test_pbd` | PBD 重力、接触投影、摩擦与堆叠 |
| `test_triggers` | trigger enter/exit、方向、销毁清理 |
| `test_shapes` | 胶囊、凸包、索引三角网格窄相位 |
| `test_broadphase` | spatial-hash broadphase、AABB 查询、大代理 |

运行全部测试：

```bash
cmake --build build
ctest --test-dir build --output-on-failure
```

## 模块

- `butter/math` - 向量、矩阵、四元数、变换、包围盒
- `butter/core` - 世界、刚体、Builder、材质、事件
- `butter/shapes` - 球、盒、胶囊、凸包、索引网格与碰撞检测
- `butter/physics2d` - 2D 圆、盒、多边形、重力与冲量接触
- `butter/physics3d` - 现有 3D Butter API 的显式入口
- `butter/constraints` - 距离、弹簧、铰链关节
- `butter/query` - 射线检测与重叠查询

## Broadphase、PBD、Trigger 与碰撞形状

World 默认启用动态 spatial-hash broadphase。它按刚体聚合 AABB 入格，
再经过精确 AABB 和窄相位过滤；超大地面等代理会自动进入 large-proxy
回退列表。可以按场景调整：

```cpp
World::Config config;
config.enable_broadphase = true;
config.broadphase_cell_size = 2.0f;
config.broadphase_fat_margin = 0.05f;
World world(config);

auto& hull = world.create_body().dynamic()
    .convex({{-1, 0, -1}, {1, 0, -1}, {0, 1, 0}, {0, 0, 1}})
    .build();

std::vector<Vec3> vertices{{-10, 0, -10}, {10, 0, -10},
                           {10, 0, 10}, {-10, 0, 10}};
std::vector<MeshCollider::Triangle> triangles{{0, 1, 2}, {0, 2, 3}};
world.create_body().static_body().mesh(std::move(vertices), std::move(triangles)).build();
```

`ConvexCollider` 使用点集的凸包 support mapping（GJK/EPA），`MeshCollider`
使用索引三角形窄相位，支持球、胶囊、盒、凸包和网格之间的检测。
Trigger 只在状态转换时回调：`event.is_enter` 表示进入，`event.is_exit`
表示离开；持续重叠不会每帧重复触发。

`broadphase_candidate_count()` 会返回最近一次 broadphase 产生的刚体对
数量；`broadphase_max_cells_per_body` 限制单个刚体最多占用的网格单元，
避免超大地面代理淹没哈希表。

可以显式配置 PBD 和 Trigger 生命周期：

```cpp
World::Config config;
config.solver_mode = World::SolverMode::PBD;
config.position_iterations = 8;
World world(config);

world.on_trigger = [](const TriggerEvent& event) {
    if (event.is_enter) std::cout << "enter\n";
    if (event.is_exit)  std::cout << "exit\n";
};

world.create_body().static_body().at(0, 1, 0)
    .sphere(3.0f).trigger().build();
```

胶囊碰撞按“有限线段 + 半径”计算，包含胶囊-盒和胶囊-网格情况。凸包
使用 GJK/EPA；网格窄相位使用索引三角形，并用最近点、SAT、GJK 处理球、
胶囊、盒、凸包和网格之间的检测；无效或退化三角形会被忽略。碰撞体
offset 和旋转 OBB AABB 都在刚体局部变换下计算。持续重叠在 enter 后保持
安静，不会每帧重复回调。

当前 spatial-hash 是刚体级 broadphase；网格窄相位会遍历其索引三角形，
针对超大网格的内部 BVH 仍属于后续性能优化项。

## 实际验证

```bash
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

当前实现已通过全部 9 个测试程序。Windows 下查看 PBD 箱子堆和爆炸示例：

```powershell
cmake --build build --config Release --target exploding_crates
.\build\examples\Release\exploding_crates.exe
```

示例会先让受重力影响的箱子堆稳定，再自动施加爆炸冲量；按 `Space` 重复
爆炸，`R` 重置，`P` 暂停。

## 2D 模块

2D API 位于 `butter::physics2d`，与 3D API 并列，不会改变现有 3D 代码：

```cpp
#include <butter/physics2d/butter2d.h>
using namespace butter::physics2d;

World world;
world.create_body().static_body().at(0, -1).box(20, 1).build();
auto& ball = world.create_body().dynamic().at(0, 5).circle(0.5f).build();
for (int i = 0; i < 120; ++i) world.step();
```

当前 2D 模块包含圆、盒、胶囊、凸多边形和索引三角网格窄相位，
spatial-hash broadphase、PBD/冲量投影、摩擦、角运动、睡眠、碰撞过滤、
enter/exit trigger、raycast/AABB 查询、距离约束，以及
`physics2d_playground` 示例。`bench_2d` 会测量 100 刚体堆叠性能；CCD 和
专用网格内部 BVH 仍属于后续优化。

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
- PBD 求解模式、动态 spatial-hash broadphase、trigger enter/exit 状态
- 胶囊、凸包（GJK/EPA）和索引三角网格窄相位碰撞
- 爆炸 + 动态破碎示例

### 近期

- CCD 连续碰撞检测，彻底解决高速穿模
- 堆叠稳定性和接触求解质量提升
- 更多关节：滑轨、固定、马达

### 中期

- 多线程求解
- 简单软体、布料、弹簧骨骼
- 车辆与角色控制器

### 远期

- 流体与粒子
- 软体 FEM
- 可视化调试器与编辑器

## 许可证

本项目使用 [MIT License](LICENSE)。
