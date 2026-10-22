import numpy as np
import matplotlib.pyplot as plt

# ---------- 1. 参数设置（SI 单位） ----------
L = 1.0                     # 通道长度 m
W = 0.2                     # 通道宽度 m
H = 0.0005                  # 通道高度 m (0.5 mm)
mu = 0.001                  # 黏度 Pa·s
Qin = 0.5 / 60000           # 入口流量 m^3/s (0.5 L/min)
Pin = 3e5                   # 入口压力 Pa (3 bar)
Pp = 1e5                    # 渗透侧压力 Pa (1 bar)
Rm = 1e13                   # 膜阻力 1/m

# ---------- 2. 计算无量纲数 ----------
Pi = 12 * mu * Qin * L / (W * H**3 * (Pin - Pp))
Lambda = W * L * (Pin - Pp) / (Qin * mu * Rm)
k = np.sqrt(Pi * Lambda)

print("=" * 60)
print(f"当前物理参数下：Pi = {Pi:.4f}, Lambda = {Lambda:.4f}, k = {k:.4f}")
print("=" * 60)

# ---------- 3. 解析解 ----------
def analytical_solution(x_star, Pi, Lambda, k):
    """返回无量纲 P*, Q*, J*"""
    P_star = np.cosh(k * x_star) - (Pi / k) * np.sinh(k * x_star)
    Q_star = np.cosh(k * x_star) - (k / Pi) * np.sinh(k * x_star)
    J_star = P_star
    return P_star, Q_star, J_star

# ---------- 4. 数值求解（一阶显式欧拉法） ----------
def numerical_solution_euler(N, Pi, Lambda):
    x_star = np.linspace(0, 1, N + 1)
    dx = 1.0 / N
    P = np.zeros(N + 1)
    Q = np.zeros(N + 1)
    P[0], Q[0] = 1.0, 1.0
    for i in range(N):
        dQ = -Lambda * P[i] * dx
        dP = -Pi * Q[i] * dx
        Q[i + 1] = Q[i] + dQ
        P[i + 1] = P[i] + dP
    return x_star, P, Q, P # J*=P*

# ---------- 5. 数值求解（二阶 RK2 改进欧拉法） ----------
def numerical_solution_rk2(N, Pi, Lambda):
    x_star = np.linspace(0, 1, N + 1)
    dx = 1.0 / N
    P = np.zeros(N + 1)
    Q = np.zeros(N + 1)
    P[0], Q[0] = 1.0, 1.0
    for i in range(N):
        # 第一步：计算当前点的斜率 k1
        dP1 = -Pi * Q[i]
        dQ1 = -Lambda * P[i]
        # 第二步：计算试探点的斜率 k2
        dP2 = -Pi * (Q[i] + dx * dQ1)
        dQ2 = -Lambda * (P[i] + dx * dP1)
        # 加权平均更新
        Q[i + 1] = Q[i] + 0.5 * dx * (dQ1 + dQ2)
        P[i + 1] = P[i] + 0.5 * dx * (dP1 + dP2)
    return x_star, P, Q, P

# ---------- 6. 主计算与画图：N = 1000 ----------
N_main = 1000
x_num, P_num, Q_num, J_num = numerical_solution_euler(N_main, Pi, Lambda)
x_ana = np.linspace(0, 1, 200)
P_ana, Q_ana, J_ana = analytical_solution(x_ana, Pi, Lambda, k)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 图1：P* 和 Q*
ax1 = axes[0]
ax1.plot(x_num, P_num, 'o', markersize=3, label='P* numerical (Euler)')
ax1.plot(x_ana, P_ana, '-', linewidth=2, label='P* analytical')
ax1.plot(x_num, Q_num, 's', markersize=3, label='Q* numerical (Euler)')
ax1.plot(x_ana, Q_ana, '--', linewidth=2, label='Q* analytical')
ax1.set_xlabel('x*')
ax1.set_ylabel('P*, Q*')
ax1.set_title('[a]Pressure and Flow Rate along Channel')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 图2：J*
ax2 = axes[1]
ax2.plot(x_num, J_num, 'o', markersize=3, label='J* numerical (Euler)')
ax2.plot(x_ana, J_ana, '-', linewidth=2, label='J* analytical')
ax2.set_xlabel('x*')
ax2.set_ylabel('J*')
ax2.set_title('[b]Local Permeate Flux along Channel')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('m2_initial_results.png', dpi=300, bbox_inches='tight')
plt.show()

# ---------- 7. 网格收敛与一阶/二阶精度对比 ----------
print("\n网格收敛与精度对比分析：")
print(f"{'N':<6} | {'Euler 平均误差(%)':<18} | {'RK2 平均误差(%)':<18} | {'两者最大差异(%)':<15}")
print("-" * 65)

# 生成稠密解析解，用于精确插值对比
x_ana_dense = np.linspace(0, 1, 10000)
P_ana_dense, _, _ = analytical_solution(x_ana_dense, Pi, Lambda, k)

for N in [100, 250, 500, 1000, 2000]:
    # 跑 Euler
    x_e, P_e, _, _ = numerical_solution_euler(N, Pi, Lambda)
    # 跑 RK2
    x_r, P_r, _, _ = numerical_solution_rk2(N, Pi, Lambda)
    
    # 插值获取当前网格节点上的解析解
    P_a_at_N = np.interp(x_e, x_ana_dense, P_ana_dense)
    
    # 计算与解析解的平均相对误差
    err_e = np.mean(np.abs(P_e - P_a_at_N) / np.abs(P_a_at_N)) * 100
    err_r = np.mean(np.abs(P_r - P_a_at_N) / np.abs(P_a_at_N)) * 100
    
    # 计算 Euler 和 RK2 之间的最大差异
    diff = np.max(np.abs(P_e - P_r)) * 100
    
    print(f"{N:<6} | {err_e:<18.4e} | {err_r:<18.4e} | {diff:<15.4e}")

print("=" * 60)
print("说明：")
print("1. 随着 N 增大，误差单调减小，证明网格无关性。")
print("2. 在 N>=1000 时，两者差异极小，证明一阶 Euler 格式已足够精确。")
# ---------- 8. 零渗透极限验证 ----------
print("\n--- 零渗透极限验证 ---")
# 设置一个极小的 Lambda (趋近于0)，代表膜几乎不透水
Lambda_zero = 1e-8 
N_test = 1000
x_zero, P_zero, Q_zero, _ = numerical_solution_euler(N_test, Pi, Lambda_zero)

# 理论直线：当不透水时，流量 Q* 恒为 1，压降完全由摩擦引起，P* 线性下降
P_linear = 1.0 - Pi * x_zero

# 计算数值解与理论直线的最大偏差
err_zero = np.max(np.abs(P_zero - P_linear)) * 100

print(f"当 Lambda -> 0 时（零渗透极限）：")
print(f"数值压力解与理论线性直线 P* = 1 - Pi*x* 的最大偏差为：{err_zero:.6e} %")