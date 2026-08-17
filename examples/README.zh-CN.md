# 爆炸炸碎箱子堆

[English](README.md)

一个使用 Butter 物理引擎和 GLFW 渲染的示例。

## 展示内容

- 箱子从空中受重力掉落。
- 通过真实物理碰撞、翻滚，逐渐稳定成堆。
- 按 `空格` 触发爆炸。
- 受到的爆炸冲量超过阈值的箱子会破碎成更小的碎片。
- 碎片和未破碎的箱子沿不同轨迹飞散，随后减速并停下。

## 构建

```bash
cmake -S . -B build -DBUTTER_GLFW_DIR=E:/Github/glfw
cmake --build build --target exploding_crates
```

## 运行

```bash
./build/examples/exploding_crates
```

Windows 下运行：

```bash
build/examples/Debug/exploding_crates.exe
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

返回主 README：[../README.zh-CN.md](../README.zh-CN.md)。
