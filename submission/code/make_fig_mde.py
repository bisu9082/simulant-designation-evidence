#!/usr/bin/env python3
"""Fig S5 — 최소검출효과. 음성 결과가 검정력 부족이 아님을 보인다.

본문·SI 다른 도판과 같은 서식 상수를 쓴다 (make_fig_main.py 와 동일).
"""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import fig_annot

FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT = 15, 14, 12, 12
C_REAL, C_SIM, C_INTF, C_REF, C_INK = '#C94F4A', '#3A6EA8', '#E8943A', '#B84040', '#333333'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.spines.top': False, 'axes.spines.right': False})

M = json.load(open('outputs/mde.json'))
P35 = M['panel35_incumbent3']
P28 = M['panel28_incumbent3']

# figsize 는 \textwidth(384 pt = 5.33 in) 기준으로 정한다. 24~27 in 로
# 그리면 축소율이 0.2 가 되어 15 pt 틱이 실효 3 pt 로 인쇄된다 (2026-10-20
# 외부 검토 지적, 전 도판 해당). 10 in 안쪽으로 그려 실효 7 pt 이상을 둔다.
fig, ax = plt.subplots(figsize=(8.6, 5.4))

for key, dat, col, lbl in [
        ('panel35_incumbent3', P35, C_SIM,
         'Thirty-five-compound panel (primary)'),
        ('panel28_incumbent3', P28, C_INTF,
         'Twenty-eight-compound panel')]:
    d = np.array([r['delta'] for r in dat['curve']])
    p = np.array([r['p'] for r in dat['curve']])
    ax.plot(d, p, color=col, lw=3.3, zorder=3, label=lbl)

ax.axhline(0.05, color=C_REF, ls=(0, (5, 3)), lw=2.3, zorder=2)
ax.text(1.46, 0.057, r'$\alpha = 0.05$', fontsize=FS_ANNOT, color=C_REF,
        ha='right', va='bottom')

# 최소검출효과 (주 패널)
mde = P35['mde_delta']
ax.plot([mde, mde], [0, 0.05], color=C_INK, ls=':', lw=2.3, zorder=2)
ax.scatter([mde], [0.05], s=85, facecolor='white', edgecolor=C_INK,
           lw=2.3, zorder=5)
ax.annotate(f'minimum detectable\neffect, {mde:.2f}',
            xy=(mde, 0.05), xytext=(mde + 0.10, 0.175),
            fontsize=FS_ANNOT, color=C_INK, va='center',
            arrowprops=dict(arrowstyle='-', color=C_INK, lw=1.5))

# 실제로 존재하는 결손. 이게 해석의 핵심이다.
gmed, gmin = P35['gap_median'], P35['gap_min']
ax.axvspan(mde, gmed, color=C_SIM, alpha=0.10, zorder=0)
# 라벨은 자기 선에 바로 붙여 세로로 둔다. 곡선은 이 구간에서 0.15 아래라
# 겹치지 않는다.
for x, lab in [(gmed, f'median deficit, {gmed:.2f}'),
               (gmin, f'smallest deficit (sarin), {gmin:.2f}')]:
    ax.axvline(x, color=C_REAL, lw=2.0, alpha=0.85, zorder=2)
    ax.text(x - 0.020, 0.205, lab, fontsize=FS_ANNOT - 1, color=C_REAL,
            va='bottom', ha='left', rotation=90)

ax.set_xlabel('Fidelity advantage added to the incumbent simulants, '
              r'$\delta$' '\n(typical nearest-neighbour separations)')
ax.set_ylabel('One-sided exact $p$ (mean difference)')
ax.set_xlim(0, 1.5)
ax.set_ylim(0, 0.62)
ax.legend(handles=[Line2D([], [], color=C_SIM, lw=3.3,
                          label='Thirty-five-compound panel (primary)'),
                   Line2D([], [], color=C_INTF, lw=3.3,
                          label='Twenty-eight-compound panel')],
          loc='upper left', frameon=False)

fig.tight_layout()
fig.savefig('outputs/Fig_mde.png', dpi=200, bbox_inches='tight')
print('[fig] outputs/Fig_mde.png')

for k, v in [('mde_delta35', f'{mde:.2f}'),
             ('gap_median35', f'{gmed:.3f}'),
             ('gap_min35', f'{gmin:.3f}'),
             ('p_at_zero35', f"{P35['p_at_zero']:.4f}"),
             ('p_at_zero28', f"{P28['p_at_zero']:.4f}")]:
    fig_annot.record('Fig_mde', k, v)
