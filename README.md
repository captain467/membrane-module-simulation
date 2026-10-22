# 平板膜组件模型数值模拟代码

本仓库包含平板膜组件的一维流体动力学与传质数值求解代码。

## 包含文件
- `membrane_schematic.py`: 绘制平板膜组件几何示意图。
- `membrane_solver.py`: 数值解（Euler, RK2）与解析解对比，零渗透极限验证。
- `membrane_sweep.py`: 参数扫描（Pi, Lambda, k），生成通量分布、非均匀性、回收率及反渗临界相图。

## 环境依赖
请安装 Python 3.8+ 及以下依赖：
`pip install numpy matplotlib`

## 运行方法
在终端中依次运行：
`python membrane_schematic.py`
`python membrane_solver.py`
`python membrane_sweep.py`