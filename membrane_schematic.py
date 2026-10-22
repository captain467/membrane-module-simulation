import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.lines import Line2D

# ============================================================
# 平板膜组件模型示意图（英文版，避免乱码）
# ============================================================

fig, ax = plt.subplots(figsize=(12, 6))

# ---------- 几何参数 ----------
x0, x1 = 0.5, 11.5          # 通道左右边界
y_bottom = 3.0              # 通道下边界（膜面）
y_top = 4.5                 # 通道上边界
y_perm_bottom = 1.5         # 渗透侧下边界
y_perm_top = y_bottom       # 渗透侧上边界（膜面）

# ---------- 1. 料液通道 ----------
channel = patches.Rectangle(
    (x0, y_bottom), x1 - x0, y_top - y_bottom,
    linewidth=1.5, edgecolor='black', facecolor='#D6EAF8'
)
ax.add_patch(channel)

# ---------- 2. 膜（用粗线表示） ----------
ax.plot([x0, x1], [y_bottom, y_bottom],
        color='darkred', linewidth=4, label='Membrane')

# ---------- 3. 渗透侧 ----------
perm = patches.Rectangle(
    (x0, y_perm_bottom), x1 - x0, y_perm_top - y_perm_bottom,
    linewidth=1.5, edgecolor='black', facecolor='#FDEBD0'
)
ax.add_patch(perm)

# ---------- 4. 通道内主流箭头 ----------
n_arrows = 7
xs = np.linspace(x0 + 0.6, x1 - 0.6, n_arrows)
for xa in xs:
    ax.annotate(
        '', xy=(xa + 0.6, (y_bottom + y_top) / 2),
        xytext=(xa, (y_bottom + y_top) / 2),
        arrowprops=dict(arrowstyle='->', color='#1F618D', lw=2)
    )

# ---------- 5. 局部渗透通量箭头 ----------
n_flux = 9
xs_flux = np.linspace(x0 + 0.4, x1 - 0.4, n_flux)
for i, xf in enumerate(xs_flux):
    # 通量沿程递减：入口长，出口短
    length = 0.55 * (1 - 0.7 * i / (n_flux - 1))
    ax.annotate(
        '', xy=(xf, y_bottom - length),
        xytext=(xf, y_bottom),
        arrowprops=dict(arrowstyle='->', color='#C0392B', lw=1.8)
    )
    ax.text(xf, y_bottom - length - 0.15,
            f'$J_{i+1}$', ha='center', va='top',
            fontsize=9, color='#922B21')

# ---------- 6. 入口标注 ----------
ax.annotate(
    '', xy=(x0 + 0.6, (y_bottom + y_top) / 2),
    xytext=(x0 - 0.6, (y_bottom + y_top) / 2),
    arrowprops=dict(arrowstyle='->', color='#1F618D', lw=2.5)
)
ax.text(x0 - 0.7, (y_bottom + y_top) / 2 + 0.15,
        r'\(Q_{in},\ P_{in}\)', ha='right', va='bottom',
        fontsize=12, color='#1F618D')

# ---------- 7. 出口标注 ----------
ax.annotate(
    '', xy=(x1 + 0.6, (y_bottom + y_top) / 2),
    xytext=(x1 - 0.4, (y_bottom + y_top) / 2),
    arrowprops=dict(arrowstyle='->', color='#1F618D', lw=2.5)
)
ax.text(x1 + 0.7, (y_bottom + y_top) / 2 + 0.15,
        r'\(Q_{out},\ P_{out}\)', ha='left', va='bottom',
        fontsize=12, color='#1F618D')

# ---------- 8. 区域文字 ----------
ax.text((x0 + x1) / 2, (y_bottom + y_top) / 2 + 0.35,
        'Retentate Channel',
        ha='center', va='center', fontsize=12, color='#1B4F72')

ax.text((x0 + x1) / 2, (y_perm_bottom + y_perm_top) / 2,
        r'Permeate, \(P_p\) constant',
        ha='center', va='center', fontsize=12, color='#7E5109')

# ---------- 9. 尺寸标注 ----------
# 长度 L
ax.annotate(
    '', xy=(x0, y_top + 0.35), xytext=(x1, y_top + 0.35),
    arrowprops=dict(arrowstyle='<->', color='black', lw=1.2)
)
ax.text((x0 + x1) / 2, y_top + 0.45, r'\(L\)',
        ha='center', va='bottom', fontsize=12)

# 高度 H
ax.annotate(
    '', xy=(x1 + 0.25, y_bottom), xytext=(x1 + 0.25, y_top),
    arrowprops=dict(arrowstyle='<->', color='black', lw=1.2)
)
ax.text(x1 + 0.35, (y_bottom + y_top) / 2, r'\(H\)',
        ha='left', va='center', fontsize=12)

# 宽度 W（示意，标在通道上方）
ax.text(x0 + 0.1, y_top + 0.05, r'Width \(W\) (into page)',
        ha='left', va='bottom', fontsize=10, color='gray')

# ---------- 10. 坐标轴 ----------
# x 方向
ax.annotate(
    '', xy=(x1 + 0.3, y_perm_bottom - 0.8),
    xytext=(x0, y_perm_bottom - 0.8),
    arrowprops=dict(arrowstyle='->', color='black', lw=1.2)
)
ax.text(x1 + 0.4, y_perm_bottom - 0.8, r'\(x\)',
        ha='left', va='center', fontsize=12)

# ---------- 11. 图例 ----------
legend_elements = [
    Line2D([0], [0], color='darkred', lw=4, label='Membrane'),
    Line2D([0], [0], color='#1F618D', lw=2, label='Main flow'),
    Line2D([0], [0], color='#C0392B', lw=2, label='Local flux \(J(x)\)'),
]
ax.legend(handles=legend_elements, loc='upper center',
          bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=10)

# ---------- 12. 版面 ----------
ax.set_xlim(-1.5, 13.5)
ax.set_ylim(0.5, 6.0)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig('membrane_module_schematic.png', dpi=300, bbox_inches='tight')
plt.show()