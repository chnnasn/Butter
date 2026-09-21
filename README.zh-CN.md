# Butter

> **Butter** - Make Physics Smooth Again!
> 一个现代、灵活、符合直觉的 C++20 header-only 物理引擎。

[English](README.md)

## 简介

Butter 追求像 C#/TypeScript 一样流畅的 API，同时保持 C++ 的高性能和零开销抽象。

项目同时提供 3D 与 2D API。下方快速开始使用 3D；最新优化和
[基准结果](#2d-基准测试报告2026-09-21)针对 2D。

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
        .bounciness(0.7f)
        .build();

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
| `test_2d` | 2D 刚体、碰撞、查询与约束 |
| `test_2d_ccd` | 薄墙、旋转/偏移扫掠、复合刚体位置修正、4–64 顶点多边形与失败隔离 |
| `test_2d_engine` | 显式 fixture、引擎接入与接触生命周期 |
| `test_3d_namespace` | 显式 3D 命名空间兼容性 |
| `test_2d_stability` | 大规模落地、箱堆长期静置与关联唤醒/休眠 |
| `test_2d_tomcat_stack` | TomCat 原始箱体参数、实际时间推进与地板检查 |
| `test_2d_tomcat_circles` | TomCat 1,000 圆形回归（运行 `test_2d_tomcat_stack 1000 circles`） |
| `test_2d_optimization` | 宽相暴力对照、缓存失效、全休眠快速路径、运动学运动与顶点/CCD 缓冲 |

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

`a8e4336` 的 Release、Debug 均通过全部 17 个 CTest 用例（16 个测试程序）。
这些是仓库回归测试，与下方外部基准测试分开记录。Windows 下查看 PBD 箱子堆和爆炸示例：

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
spatial-grid broadphase、PBD/冲量投影、摩擦、角运动、睡眠、碰撞过滤、
enter/exit trigger、raycast/AABB 查询、距离约束，以及
`physics2d_playground` 示例。`bench_2d [count] [circles|boxes] [awake]` 检查落地正确性，
分别报告活动与休眠阶段耗时；`awake` 显式保持刚体活动，以便对照。
这个内部落地基准与下方 TomCat 箱堆基准属于不同场景。
圆、盒和凸多边形已支持保守扫掠 CCD；专用网格内部 BVH 仍属于后续优化。

### 使用显式 fixture 接入引擎

`World::create_empty_body()` 创建没有 Builder 默认形状的动态刚体。
`add_fixture(body, definition)` 提供稳定且可独立移除的 fixture，支持局部变换、材质、trigger 和类别/掩码。
不要给 Builder 创建的单形状刚体再添加显式 fixture。

`contact_filter` 提供 fixture 级过滤；`on_contact(a, b, enter)` 报告显式 fixture 的实体接触与触发器状态转换。
休眠刚体保留接触；销毁 fixture 或刚体前会发送退出事件。回调期间 `locked() == true`，
世界结构修改应排队到回调或 `step()` 返回后。Builder 刚体继续使用 `on_trigger`；
Builder 与显式 fixture 的混合碰撞对参与求解，但不发布转换回调。

`destroy_fixture`、`destroy_joint`、`destroy_body` 分别释放对象，销毁刚体也会移除其关节。
`Body::user_data` 和 `Fixture::user_data` 为显式 64 位值。接入层可使用 `force`、`torque`、
`fixed_rotation` 和 `BodyType::Kinematic`；直接修改质量或类型时，必须同步质量和惯量的倒数。
显式复合刚体的 fixture 密度仅作为元数据，其质量与惯量由接入层提供。

距离关节支持旋转局部锚点与 `collide_connected`；软关节使用 `spring_stiffness` 和 `damping`。
`Config::max_position_correction` 限制深穿透修正。射线查询检测圆与凸多边形的实际形状，
包含旋转盒子，并忽略起点位于内部的形状。

### 2D 连续碰撞检测（CCD）

默认对动态刚体与静态/运动学刚体启用 CCD。设置 `body.bullet = true`（或 Builder 的 `.bullet()`）后，
只要动态碰撞对任一方为 bullet，也会进行连续检测。`config.ccd.enabled = false` 恢复离散推进；3D 求解器不受影响。

CCD 使用扫掠包围盒与分离平面的保守推进，支持圆、盒和凸多边形，计入双方平移、未折返的角度变化及旋转 fixture 偏移。
推进到最早接触时刻后求解弹性与摩擦冲量，再检测剩余运动；力与阻尼每个外层物理步只积分一次。
fixture 掩码、自定义过滤与关节禁碰设置均生效。显式 fixture 会为 CCD 撞击发送状态转换事件，
同一步内反弹离开的接触可产生一对进入/退出事件。

`config.ccd.max_impacts` 默认 32，按**运动刚体**限制撞击次数；`max_iterations` 默认 64，
`tolerance` 默认 0.0001 世界单位。同一时刻的接触成批求解；次数耗尽、未收敛或零时间重复命中时，
只保守钳制参与物体的剩余运动，其他物体继续经过 CCD 推进。
`world.ccd_statistics()` 分别记录上述失败原因、刚体/碰撞体索引、已推进时间和局部剩余时间。
`remaining_time` 表示最大的局部钳制区间，不是整个世界丢弃的时间；局部钳制仍是保守降级，不能视为等质量运动。

持续接触使用裁剪生成的多边形双接触点、特征/锚点匹配、累计法向与摩擦冲量和热启动。
速度求解与位置修正分开；接触和关节的位置修正也会检测静态及运动学障碍，包括原候选列表外的薄墙。
接触/关节相连的动态刚体整体唤醒与休眠，共用静态地板不会把无关物体并成一组。

胶囊和网格仍采用离散检测；传感器仅检测步末重叠，不阻止 CCD 运动。直接设置位置的瞬移不做扫掠，
初始穿透需要离散恢复。多边形必须凸且非退化，角度应保持连续（整圈是 `2*pi`，不是零）。
位置修正保护覆盖与 CCD 相同的凸形状。

原始 fixture 箱堆复现与独立静态环境时间线见 [TomCat 零时间 CCD 修复记录](docs/ccd-zero-time-followup.md)。
首次错误诊断、压力测试命令、计时边界和限制见 [2D 正确性与验证](docs/physics2d-correctness.md)。
历史耗时倍数不作为优化目标；先验证实际时间推进、穿透、休眠和活动状态是否可比。

`test_2d_ccd` 覆盖关闭 CCD 的穿透负对照、薄墙圆/盒/多边形撞击、运动学物体、bullet 对、
旋转与偏移扫掠、多次反弹、过滤、休眠、接触转换、预算回退和确定性的密集时间采样参照。
Release 构建保留这些检查。

## 2D 基准测试报告（2026-09-21）

### 最新结果：`a8e4336` 对照 `9baad05`

本轮减少 CCD 几何查询的临时分配，并在位置修正重建包围盒前排除远处障碍；可能命中时
仍执行完整 CCD 保护。此前已加入速度约束预计算、接触几何复用，以及全休眠世界跳过
重复 CCD/检测。没有降低求解次数、放宽休眠阈值或提高 CCD 预算，求解仍为串行。

测试环境与适配层：[TomCat Engine — `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter)。
本机使用该适配层的独立副本接入 Butter 头文件，不代表引擎分支更新了依赖。
Release、单线程、1/60 秒步长、预热 60 步、计时 300 步；六轮交替新旧版顺序，
舍弃首轮，取后五轮中位数。下表单位为 ms/步。

| 场景 | 刚体数 | `9baad05` | `a8e4336` |
| --- | ---: | ---: | ---: |
| 分散移动 | 100 | 0.019314 | 0.018888 |
| 分散移动 | 500 | 0.099302 | 0.098543 |
| 分散移动 | 1,000 | 0.212034 | 0.210657 |
| 圆形落地 | 100 | 0.026110 | 0.027397 |
| 圆形落地 | 500 | 0.356339 | 0.338405 |
| 圆形落地 | 1,000 | 1.139450 | 1.033050 |
| 箱体堆叠 | 100 | 0.068556 | 0.059709 |
| 箱体堆叠 | 500 | 5.405940 | 4.081260 |
| 箱体堆叠 | 1,000 | 12.166300 | 9.203220 |

500／1,000 箱体分别减少约 **25%／24%** 耗时，1,000 圆形减少约 **9%**。
分散移动基本持平；100 圆形增加约 1.29 微秒（5%），不能声称全部场景提速。
这是 Butter 版本对照，不是新的 Box2D 倍数；平均低于 16.67 ms 也不代表每步都能达到 60 Hz。

九个场景的新旧逐步状态哈希与活动刚体步数一致。全部计时运行完整推进 5 秒，
CCD 限制、预算耗尽、未收敛、零时间重复命中、非有限坐标和检测到的中心穿地均为零。
更严格的原生 500／1,000 箱体测试全程检查地面几何边界，分别约在 **9.57／12.17 秒**
全部休眠并保持到 60 秒。第六秒仍唤醒不能直接认定为永久无法休眠。
这些结果验证指定场景，不代表任意物理负载均已正确。

`world.step_statistics()` 提供 CCD、检测、速度、位置、缓存、休眠耗时及活动量，
并提供 `stationary_steps` 与 `projection_fast_rejections` 计数。
独立 1,000 箱体诊断中，位置求解约从 4.03 降至 1.49 ms/步，修正及候选数量不变。
接触检测仍是主要开销之一；后续比较仍需同时核对活动状态和轨迹。

详见[实现、失效规则与验证记录](docs/serial-optimization.zh-CN.md)和
[108 条原始样本](docs/benchmarks/projection-followup-20260921.csv)。原始数据中
`baseline9` 指 `9baad05`，`retest` 指后来提交为 `a8e4336` 的代码。
运行命令见[示例与诊断指南](examples/README.zh-CN.md#2d-回归与基准)。

### 外部 Box2D 参照：`52f15a5`

用户提供的后续 TomCat 报告在 i7-14650HX 上使用 Release、单线程及相同预热/计时步数，
比较 Butter `52f15a5` 与 Box2D 2.4.1：

| 场景（1,000 刚体） | Box2D | Butter `52f15a5` | Butter 耗时倍数 |
| --- | ---: | ---: | ---: |
| 分散移动 | 0.2118 ms | 0.2217 ms | 1.05 倍 |
| 圆形落地 | 0.5340 ms | 1.4868 ms | 2.78 倍 |
| 箱体堆叠 | 3.0651 ms | 13.6228 ms | 4.44 倍 |

无碰撞场景基本持平，接触密集场景仍有数倍差距。当时引擎分支仍固定旧依赖，
测试在独立目录中接入新版 Butter。两个引擎的轨迹和休眠进度不同，因此这是相同
初始场景的实际耗时，不是等质量吞吐量。不能用本轮本机耗时除以上述旧 Box2D 耗时，
宣称得到新的对比倍数。

<details>
<summary>历史报告：ea8ef90 与首轮串行优化</summary>

**串行优化后续复测：**同轮对照 `ea8ef90`，1,000 箱体从 **55.152 降至 13.968 ms／步**，
圆形为 **6.888 → 1.470 ms**，分散移动为 **1.351 → 0.221 ms**。
500／1,000 箱体约在 9.57／12.17 秒全部休眠，并保持到 60 秒。
缓存规则、匹配活动量的分配诊断、原始数据与限制见[中英文实现及验证记录](docs/serial-optimization.zh-CN.md)。
下文保留此前 `ea8ef90` 与 Box2D 的历史测量。

测试环境与接入代码：[TomCat Engine — `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter)。
本轮使用 **Butter `ea8ef90`** 和该分支的物理适配层，对照库为 TomCat `main` 分支使用的 **Box2D 2.4.1**。
**测试当时 `dev_butter` 仍固定旧版 Butter**：以下结果来自独立测试目录接入新版头文件，不代表该分支已经更新依赖。

硬件与构建：Intel Core i7-14650HX、Windows x64、MSVC 19.50.35724、C++20 Release（`/O2 /DNDEBUG`），单线程负载。
时间步长为 1/60 秒，每个场景先预热 60 步，再计时 300 步，默认开启 CCD 和休眠。
Box2D 使用 `Step(8,3)`，Butter 适配层映射为 8 次迭代。
六轮运行交替库的执行顺序，舍弃首轮，取后五轮中位数；计时期间没有并行测试或编译，未控制 CPU 亲和性与电源计划。
耗时不包含创建、渲染、脚本和引擎同步。

### 1,000 个动态刚体的物理步进性能

| 场景 | Box2D | Butter | Butter 耗时倍数 |
| --- | ---: | ---: | ---: |
| 分散移动、无碰撞 | 0.206 ms | 1.108 ms | 5.38 倍 |
| 圆形刚体落地 | 0.537 ms | 5.867 ms | 10.92 倍 |
| 箱体堆叠 | 2.967 ms | 51.237 ms | 17.27 倍 |

本次改善主要集中在箱体堆叠；圆形落地与无碰撞移动的性能差距基本保持原有量级。

| 箱体数量 | 上轮相对 Box2D 的耗时 | 本轮相对 Box2D 的耗时 |
| --- | ---: | ---: |
| 100 | 974.42 倍 | 35.80 倍 |
| 500 | 134.74 倍 | 28.43 倍 |
| 1,000 | 81.37 倍 | 17.27 倍 |

跨轮结果是观察值，不能直接换算成受控的版本提速倍数。两个引擎的轨迹与休眠状态仍有差异，
当前倍数代表相同初始场景下的实际耗时，不是完全相同模拟质量下的吞吐比较。

### 正确性与稳定性复查

对 100、500、1,000 刚体的圆形落地和箱体堆叠，共六个场景进行额外诊断：

| 检查项 | 本轮结果 |
| --- | --- |
| 300 个测量步的模拟时间推进 | 六个场景均完整推进 5 秒 |
| CCD 限制触发 | 0 |
| 碰撞预算耗尽 | 0 |
| CCD 未收敛 | 0 |
| 零时间重复碰撞 | 0 |
| 地面覆盖范围内检测到的穿地 | 0 |
| 非有限坐标 | 0 |

圆形场景和 100 箱体场景能够休眠；500、1,000 箱体在测试结束时仍保持唤醒。
本轮穿地诊断检查 `abs(x) < 99` 范围内的刚体中心，不是完整的几何穿透证明。
此前复现的整个世界停止推进、穿地和零时间重复碰撞问题，在本轮同类用例中均未再次出现。
这是限定场景的验证，不能替代完整物理正确性测试。

历史版本 `ea8ef90` 的 1,000 箱体单步约 **51 ms**，当时超过 60 Hz 的 **16.67 ms** 预算。
当前性能及休眠结论以上方最新测量为准。

本轮源码、原始数据和复现方法位于本地 TomCat 检出的 `.scratch/comparison-20260921-r2/README.md`，
并未通过上面的分支链接公开。仓库内的 [CCD 修复记录](docs/ccd-zero-time-followup.md) 提供实现诊断与回归覆盖；
该记录中的本地计时与本报告属于不同测试轮次。

</details>

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
- 2D 凸形状 CCD、持续接触流形、热启动与关联休眠
- TomCat 原始圆/箱体回归及区分活动状态的计时诊断
- 速度几何预计算、小型 CCD 顶点缓冲与全休眠快速路径
- 包含复合 fixture 偏移的位置修正保守预筛选

### 近期

- 将 CCD 扩展到 3D、胶囊和网格
- 堆叠稳定性和接触求解质量提升
- 在匹配活动状态的基准下减少密集接触检测和缓存开销
- 优化串行流程时保持薄墙 CCD 与长期静置稳定性
- 更多关节：滑轨、固定、马达

### 中期

- 串行正确性与性能分析完成后再引入多线程求解
- 简单软体、布料、弹簧骨骼
- 车辆与角色控制器

### 远期

- 流体与粒子
- 软体 FEM
- 可视化调试器与编辑器

## 许可证

本项目使用 [MIT License](LICENSE)。
