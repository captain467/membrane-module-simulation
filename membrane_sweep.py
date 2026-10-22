import numpy as np
import matplotlib.pyplot as plt

# 设置全局字体
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'stix' 

# ============================================================
# 核心求解器 (所有工况统一调用此函数)
# ============================================================
def solve_membrane(Pi, Lambda, N=2000):
    x = np.linspace(0, 1, N + 1)
    dx = 1.0 / N
    P = np.zeros(N + 1)
    Q = np.zeros(N + 1)
    J = np.zeros(N + 1)
    
    P[0] = 1.0
    Q[0] = 1.0
    
    for i in range(N):
        J[i] = P[i] 
        dQ = -Lambda * P[i] * dx
        dP = -Pi * Q[i] * dx
        Q[i + 1] = Q[i] + dQ
        P[i + 1] = P[i] + dP
        
    J[N] = P[N]
    
    U = (J.max() - J.min()) / J.mean()
    recovery = np.sum(Lambda * J[:-1]) * dx 
    
    return x, P, Q, J, U, recovery

# ============================================================
# Part 1: 综合参数 k 的扫描 (原 M2 的 3 张图)
# ============================================================
print("="*50)
print("正在计算 Part 1: 综合参数 k 扫描...")
original_cases = [
    (0.02, 0.50), (0.13, 0.69), (0.35, 0.71), (0.70, 0.91),
    (0.95, 1.05), (1.45, 1.17), (1.70, 1.32)
]

results_k = []
for Pi, Lam in original_cases:
    x, P, Q, J, U, rec = solve_membrane(Pi, Lam)
    k = np.sqrt(Pi * Lam)
    results_k.append({'x': x, 'J': J, 'k': k, 'U': U, 'rec': rec})

# --- 图1 ---
fig1, ax1 = plt.subplots(figsize=(7, 5))
for r in results_k:
    ax1.plot(r['x'], r['J'], linewidth=2, label=f'k = {r["k"]:.2f}')
ax1.set_xlabel(r'Dimensionless Position $x^*$')
ax1.set_ylabel(r'Dimensionless Flux $J^*$')
ax1.set_title('Local Permeate Flux Distribution vs k')
ax1.legend()
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig1_orig_flux_vs_k.png', dpi=300)
plt.close(fig1)

# --- 图2 ---
fig2, ax2 = plt.subplots(figsize=(7, 5))
ks = [r['k'] for r in results_k]
Us = [r['U'] for r in results_k]
ax2.plot(ks, Us, 'o-', linewidth=2, markersize=8, color='C3')
ax2.set_xlabel(r'Coupling Strength $k$')
ax2.set_ylabel(r'Non-uniformity $U$')
ax2.set_title('Flux Non-uniformity vs k')
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig2_orig_U_vs_k.png', dpi=300)
plt.close(fig2)

# --- 图3 ---
fig3, ax3 = plt.subplots(figsize=(7, 5))
recs = [r['rec'] for r in results_k]
ax3.plot(ks, recs, 's-', linewidth=2, markersize=8, color='C4')
ax3.axvline(1.0, color='gray', linestyle='--', label='Peak approx. at k=1')
ax3.set_xlabel(r'Coupling Strength $k$')
ax3.set_ylabel(r'Recovery $Q_p/Q_{in}$')
ax3.set_title('Permeate Recovery vs k')
ax3.legend()
ax3.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig3_orig_recovery_vs_k.png', dpi=300)
plt.close(fig3)

# ============================================================
# Part 2: 独立因素扫描
# ============================================================
print("正在计算 Part 2: 独立因素扫描...")

# --- 2.1 固定 Lambda = 1.0, 扫描 Pi ---
Lambda_fixed = 1.0
Pi_list = [0.2, 0.5, 1.0, 1.5, 1.65, 2.0]
res_sweep_Pi = []
for Pi in Pi_list:
    _, _, _, J, U, rec = solve_membrane(Pi, Lambda_fixed)
    res_sweep_Pi.append({'Pi': Pi, 'J': J, 'U': U, 'rec': rec})

# 图4
fig4, ax4 = plt.subplots(figsize=(7, 5))
x = np.linspace(0, 1, 2001)
for r in res_sweep_Pi:
    k = np.sqrt(r['Pi'] * Lambda_fixed)
    ax4.plot(x, r['J'], linewidth=2, label=rf'$\Pi={r["Pi"]}, k={k:.2f}$')
ax4.axhline(0, color='black', linestyle='--', linewidth=1)
ax4.set_xlabel(r'Dimensionless Position $x^*$')
ax4.set_ylabel(r'Dimensionless Flux $J^*$')
ax4.set_title(r'Effect of $\Pi$ (Fixed $\Lambda=1.0$)')
ax4.legend()
ax4.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig4_sweep_Pi_flux.png', dpi=300)
plt.close(fig4)

# 图5
fig5, ax5 = plt.subplots(figsize=(7, 5))
U_Pi = [r['U'] for r in res_sweep_Pi]
rec_Pi = [r['rec'] for r in res_sweep_Pi]
ax5.plot(Pi_list, U_Pi, 'o-', color='C3', label=r'Non-uniformity $U$')
ax5.set_xlabel(r'$\Pi$ (Channel pressure drop parameter)')
ax5.set_ylabel(r'Non-uniformity $U$', color='C3')
ax5.tick_params(axis='y', labelcolor='C3')

ax5_twin = ax5.twinx()
ax5_twin.plot(Pi_list, rec_Pi, 's--', color='C4', label='Recovery')
ax5_twin.set_ylabel(r'Recovery $Q_p/Q_{in}$', color='C4')
ax5_twin.tick_params(axis='y', labelcolor='C4')
ax5.set_title(r'Metrics vs $\Pi$ (Fixed $\Lambda=1.0$)')
ax5.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig5_sweep_Pi_metrics.png', dpi=300)
plt.close(fig5)

# --- 2.2 固定 Pi = 1.0, 扫描 Lambda ---
Pi_fixed = 1.0
Lambda_list = [0.2, 0.5, 1.0, 1.5, 2.0, 3.0]
res_sweep_Lam = []
for Lam in Lambda_list:
    _, _, _, J, U, rec = solve_membrane(Pi_fixed, Lam)
    res_sweep_Lam.append({'Lam': Lam, 'J': J, 'U': U, 'rec': rec})

# 图6
fig6, ax6 = plt.subplots(figsize=(7, 5))
for r in res_sweep_Lam:
    k = np.sqrt(Pi_fixed * r['Lam'])
    ax6.plot(x, r['J'], linewidth=2, label=rf'$\Lambda={r["Lam"]}, k={k:.2f}$')
ax6.axhline(0, color='black', linestyle='--', linewidth=1)
ax6.set_xlabel(r'Dimensionless Position $x^*$')
ax6.set_ylabel(r'Dimensionless Flux $J^*$')
ax6.set_title(r'Effect of $\Lambda$ (Fixed $\Pi=1.0$)')
ax6.legend()
ax6.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig6_sweep_Lambda_flux.png', dpi=300)
plt.close(fig6)

# 图7
fig7, ax7 = plt.subplots(figsize=(7, 5))
U_Lam = [r['U'] for r in res_sweep_Lam]
rec_Lam = [r['rec'] for r in res_sweep_Lam]
ax7.plot(Lambda_list, U_Lam, 'o-', color='C3', label=r'Non-uniformity $U$')
ax7.set_xlabel(r'$\Lambda$ (Membrane permeability parameter)')
ax7.set_ylabel(r'Non-uniformity $U$', color='C3')
ax7.tick_params(axis='y', labelcolor='C3')

ax7_twin = ax7.twinx()
ax7_twin.plot(Lambda_list, rec_Lam, 's--', color='C4', label='Recovery')
ax7_twin.set_ylabel(r'Recovery $Q_p/Q_{in}$', color='C4')
ax7_twin.tick_params(axis='y', labelcolor='C4')
ax7.set_title(r'Metrics vs $\Lambda$ (Fixed $\Pi=1.0$)')
ax7.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig7_sweep_Lambda_metrics.png', dpi=300)
plt.close(fig7)

# ============================================================
# Part 3: 反渗临界相图
# ============================================================
print("正在计算 Part 3: 反渗临界边界...")

N_grid = 40
Pi_grid = np.linspace(0.01, 3.0, N_grid)
Lambda_grid = np.linspace(0.01, 3.0, N_grid)
Pi_mesh, Lambda_mesh = np.meshgrid(Pi_grid, Lambda_grid)
P_out_mesh = np.zeros_like(Pi_mesh)

for i in range(N_grid):
    for j in range(N_grid):
        _, P_out, _, _, _, _ = solve_membrane(Pi_mesh[i, j], Lambda_mesh[i, j], N=500)
        P_out_mesh[i, j] = P_out[-1]

fig8, ax8 = plt.subplots(figsize=(7, 5))
cf = ax8.contourf(Pi_mesh, Lambda_mesh, P_out_mesh, levels=[-1, 0, 1], colors=['#ff9999', '#99ff99'], alpha=0.4)
cs = ax8.contour(Pi_mesh, Lambda_mesh, P_out_mesh, levels=[0], colors='red', linewidths=2)
ax8.clabel(cs, inline=True, fontsize=10, fmt=r'$P^*(1)=0$')

ax8.text(2.0, 0.5, 'Reverse Osmosis\n(Backflow)', ha='center', color='darkred', fontweight='bold')
ax8.text(0.5, 2.5, 'Normal\nFiltration', ha='center', color='darkgreen', fontweight='bold')

for Pi, Lam in original_cases:
    ax8.plot(Pi, Lam, 'ko', markersize=5)
ax8.plot(1.70, 1.32, 'k*', markersize=12, label=r'M2 case $k=1.5$ (In RO region)')

ax8.set_xlabel(r'$\Pi$ (Channel pressure drop)')
ax8.set_ylabel(r'$\Lambda$ (Membrane permeability)')
ax8.set_title('Critical Boundary of Reverse Osmosis')
ax8.legend(loc='upper left')
ax8.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('fig8_critical_RO_boundary.png', dpi=300)
plt.close(fig8)

print("="*50)
print("全部计算完成！共生成 8 张图片（fig1 到 fig8），已保存至当前文件夹。")