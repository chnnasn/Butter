# Butter 示例与诊断

[English](README.md)

以下命令均从仓库根目录运行。Butter 的 2D 与 3D 示例使用各自的 API；
最新 2D 优化没有修改 3D 求解器。

| 目标 | 用途 | 额外依赖 |
| --- | --- | --- |
| `hello_butter` | 最小 3D 模拟 | 无 |
| `falling_boxes` | 3D 刚体落地 | 无 |
| `chain` | 3D 关节链 | 无 |
| `playground` | 3D 控制台示例 | 无 |
| `physics2d_playground` | 控制台 2D 箱体模拟 | 无 |
| `exploding_crates` | 交互式 3D 爆炸与破碎 | GLFW 和 OpenGL |

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_EXAMPLES=ON
cmake --build build --config Release
```

单配置生成器运行 `./build/examples/physics2d_playground`；
Visual Studio 运行 `./build/examples/Release/physics2d_playground.exe`。
`physics2d_playground` 是控制台示例，不含交互渲染窗口。

## 爆炸炸碎箱子堆

一个使用 Butter 3D 物理引擎和 GLFW 渲染的示例。

## 展示内容

- 箱子从空中受重力掉落。
- 通过真实物理碰撞、翻滚，逐渐稳定成堆。
- 初始三秒用于让箱堆在重力下落地并稳定；之后会自动引爆，也可以提前按 `空格`。
- 受到的爆炸冲量超过阈值的箱子会破碎成更小的碎片。
- 碎片和未破碎的箱子沿不同轨迹飞散，随后减速并停下。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_EXAMPLES=ON -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --config Release --target exploding_crates
```

## 运行

```bash
./build/examples/exploding_crates
```

Windows 下运行：

```bash
build/examples/Release/exploding_crates.exe
```

## 操作

| 输入 | 功能 |
| --- | --- |
| `空格` | 引爆 |
| `P` | 暂停 / 继续 |
| `R` | 重置场景 |
| 鼠标拖动 | 旋转视角 |
| 鼠标滚轮 | 缩放 |
| `左` / `右` | 减速 / 加速模拟 |
| `Esc` | 退出 |

## 物理细节

- 持续爆炸力场在若干帧内持续推动物体。
- 偏离质心的冲量让箱子翻滚旋转。
- 接触点相对速度包含角速度，因此旋转的碎片会被摩擦减速。
- 线性和角阻尼让碎片失去能量并最终稳定。

## 2D 基准结果

本箱堆示例使用 3D API。独立的 TomCat 2D 适配层基准，包括 `a8e4336` 结果、历史 Box2D 参照与正确性边界，
见[主 README](../README.zh-CN.md#2d-基准测试报告2026-09-21)。这些耗时不代表本示例的性能。

最新本机六轮对照 `9baad05`，1,000 箱体为 **12.166 → 9.203 ms/步**（减少约 24%），
1,000 圆形为 **1.139 → 1.033 ms/步**（减少约 9%）；100 圆形增加约 1.29 微秒。
九个场景的新旧逐步状态哈希和活动刚体步数一致；这些不是新的 Box2D 倍数。
测试环境/适配层：[TomCat Engine `dev_butter`](https://github.com/chnnasn/TomCat_Engine/tree/dev_butter)。
测试在独立目录接入新版头文件，不代表该分支更新了依赖。
详见[中英文验证及原始数据](../docs/serial-optimization.zh-CN.md)。

## 2D 回归与基准

显式开启测试和内部基准：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUTTER_BUILD_TESTS=ON -DBUTTER_BUILD_BENCHMARKS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Visual Studio / Windows 命令：

```powershell
.\build\tests\Release\test_2d_ccd.exe
.\build\tests\Release\test_2d_tomcat_stack.exe 1000 circles 360
.\build\tests\Release\test_2d_tomcat_stack.exe 500 boxes 3600
.\build\tests\Release\test_2d_tomcat_stack.exe 1000 boxes 3600
.\build\benchmarks\Release\bench_2d.exe 1000 circles
.\build\benchmarks\Release\bench_2d.exe 1000 boxes
.\build\benchmarks\Release\bench_2d.exe 1000 boxes awake
```

单配置构建去掉路径中的 `Release/` 和 `.exe`。箱堆程序参数依次是刚体数、
`circles`/`boxes` 和总步数（包括前 60 个预热步）。3,600 步用于观察 60 秒，
检查地面几何边界、时间推进和 CCD 异常。`a8e4336` 的 500／1,000 箱体约在
9.57／12.17 秒全部休眠，并保持至结束。薄墙测试还覆盖旋转、复合偏移以及
4／16／17／64 顶点多边形。Release 和 Debug 均通过全部 17 个 CTest 用例。

`bench_2d` 是独立的内部落地负载，不是 TomCat 适配层表格的复现程序；它报告
活动/休眠阶段、活动刚体步数、分阶段耗时及分配量。`awake` 模式会主动逐步唤醒
刚体，只应与同等活动状态对照。适配层表格不包含渲染与场景创建。
不能用画面流畅或 FPS 代替实际物理时间、地面检查、CCD 诊断和休眠验证。

返回主 README：[../README.zh-CN.md](../README.zh-CN.md)。
