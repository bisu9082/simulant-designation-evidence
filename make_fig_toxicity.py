#!/usr/bin/env python3
"""simulant — 독성 축 (EPA AEGL-3 10분) 사다리 + 결측 구조 명시"""
import json
import numpy as np
import pandas as pd

import fig_annot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 28, 16, 15, 13, 13, 11
C_REAL, C_SIM, C_INTF, C_REF, C_INK = '#C94F4A', '#3A6EA8', '#E8943A', '#B84040', '#333333'
CMAP = {'real_agent': C_REAL, 'simulant': C_SIM, 'interferent': C_INTF}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.facecolor': 'white', 'figure.facecolor': 'white'})
SHORT = {'sarin (GB)': 'GB  sarin', 'soman (GD)': 'GD  soman', 'tabun (GA)': 'GA  tabun',
         'cyclosarin (GF)': 'GF  cyclosarin', 'sulfur mustard (HD)': 'HD  sulfur mustard',
         'lewisite 1 (L)': 'L  lewisite 1', 'hydrogen cyanide (AC)': 'AC  hydrogen cyanide',
         'cyanogen chloride (CK)': 'CK  cyanogen chloride', 'phosgene (CG)': 'CG  phosgene',
         'chloropicrin (PS)': 'PS  chloropicrin', 'VX': 'VX',
         'dichloromethane': 'dichloromethane', 'p-dichlorobenzene': 'p-DCB'}
s = lambda x: SHORT.get(x, x)


def panel_label(ax, L):
    ax.annotate(f'({L})', xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, 6),
                textcoords='offset points', fontsize=FS_LABEL, fontweight='bold',
                va='bottom', ha='left', annotation_clip=False)


def dress(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(color='#DDDDDD', ls='--', lw=0.6); ax.set_axisbelow(True)


T = pd.read_csv('outputs/toxicity_aegl.csv')
H = T[T.AEGL3_10min_mgm3.notna()].sort_values('AEGL3_10min_mgm3', ascending=False)
M = T[T.AEGL3_10min_mgm3.isna()].copy()
order_cls = {'simulant': 0, 'real_agent': 1, 'interferent': 2}
M = M.sort_values(['compound_class', 'compound_name'],
                  key=lambda c: c.map(order_cls) if c.name == 'compound_class' else c)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(24, 10),
                             gridspec_kw=dict(width_ratios=[1.35, 1.0]))
fig.subplots_adjust(left=0.135, right=0.985, top=0.90, bottom=0.14, wspace=0.52)

# ---------------- (a) AEGL-3 사다리 (로그)
y = np.arange(len(H))
ax.barh(y, H.AEGL3_10min_mgm3, color=[CMAP[k] for k in H.compound_class],
        edgecolor='white', lw=1.0, height=0.74, zorder=2)
ax.set_yticks(y); ax.set_yticklabels([s(n) for n in H.compound_name])
ax.set_xscale('log'); ax.set_xlim(8e-3, 6e5)
ax.set_ylim(-1.7, len(H) - 0.35 + 1.25)
ax.spines['left'].set_bounds(0, len(H) - 1)
for yy, r in zip(y, H.itertuples()):
    ax.text(r.AEGL3_10min_mgm3 * 1.35, yy, f'{r.AEGL3_10min_mgm3:.3g}',
            fontsize=FS_FOOT + 1, va='center', ha='left', color='#555555')
# DCM 강조
dy = int(y[(H.compound_name == 'dichloromethane').to_numpy()][0])
dv = float(H.loc[H.compound_name == 'dichloromethane', 'AEGL3_10min_mgm3'].iloc[0])
_vx = float(H.loc[H.compound_name == 'VX', 'AEGL3_10min_mgm3'].iloc[0])
_ac = float(H.loc[H.compound_name == 'hydrogen cyanide (AC)',
                  'AEGL3_10min_mgm3'].iloc[0])
# 주석의 모든 수치와 순위는 데이터에서 계산한다. 하드코딩 금지.
# 2026-10-11 사고: 수치는 재생성됐는데 'least acutely toxic of the thirteen' 이라는
# 문자열만 구판으로 남아, 같은 패널의 막대가 메탄올과 n-헥산을 DCM 위에 그리고
# 본문이 'not the safest' 라고 쓰는 상태에서 그림만 반대를 주장했다.
_n = len(H)
_safer = H[H.AEGL3_10min_mgm3 > dv].compound_name.tolist()
_rank = int((H.AEGL3_10min_mgm3 > dv).sum()) + 1
if _safer:
    _lead = (f'dichloromethane, a designated simulant, ranks {_rank} of {_n}\n'
             f'on this axis: {" and ".join(_safer)} score higher.\n')
else:
    _lead = (f'dichloromethane, a designated simulant, is the least\n'
             f'acutely toxic of the {_n}.\n')
# 상자는 막대가 없는 좌하단에 둔다. 간섭물 구간(>400 mg/m3)은 축의 오른쪽
# 절반에서 시작하므로 아래쪽 행의 왼편이 비어 있다.
_txt_a = fig_annot.record('Fig_toxicity', 'panel_a', _lead +
            f'It is {dv/_vx:,.0f}$\\times$ less acutely toxic than VX\n'
            f'and {dv/_ac:,.0f}$\\times$ less than hydrogen cyanide.')
ax.annotate(_txt_a,
            xy=(dv * 1.15, dy), xycoords='data',
            xytext=(0.015, 0.40), textcoords='axes fraction',
            fontsize=FS_ANNOT, color=C_REF, va='top', ha='left',
            bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.35', alpha=0.96),
            arrowprops=dict(arrowstyle='-', color='none', lw=0))
ax.set_xlabel('EPA AEGL-3, 10 min (mg m$^{-3}$, log scale)\n'
              'lower = lethal at lower airborne concentration')
ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=12, mfc=c, mec='white', label=l)
                   for c, l in [(C_REAL, 'authentic agent'),
                                (C_SIM, 'designated simulant'),
                                (C_INTF, 'interferent / solvent')]],
          loc='upper right', frameon=False, bbox_to_anchor=(1.0, 1.0))
dress(ax); panel_label(ax, 'a')

# ---------------- (b) 결측 구조
y2 = np.arange(len(M))[::-1]
bx.scatter(np.zeros(len(M)), y2, marker='x', s=190,
           c=[CMAP[k] for k in M.compound_class], lw=3.0, zorder=3)
bx.set_yticks(y2[::-1]); bx.set_yticklabels([s(n) for n in M.compound_name][::-1])
REASON = {'[수집 불가: EPA AEGL 미등재]': 'not on the AEGL list at all',
          '[수집 불가: AEGL-3 = NR]': 'listed, but AEGL-3 = "not recommended,\n'
                                      'insufficient data"'}
for yy, r in zip(y2, M.itertuples()):
    bx.text(0.09, yy, REASON.get(r.note, r.note), fontsize=FS_ANNOT - 1,
            va='center', ha='left', color='#555555')
bx.set_xlim(-0.06, 1.0); bx.set_ylim(-3.0, len(M) - 0.35 + 1.25)
bx.spines['left'].set_bounds(0, len(M) - 1)
bx.set_xticks([])
bx.spines['bottom'].set_visible(False)
bx.text(-0.045, -2.05,
        'DMMP, the simulant carrying every valid match in this study, has no AEGL.\n'
        'EPA PPRTV (EPA/690/R-06/014F, 2006) declined to derive an inhalation RfC\n'
        '("data are insufficient") and derived no inhalation unit risk.',
        fontsize=FS_ANNOT, color=C_REF, va='center', ha='left',
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.4', alpha=0.93))
bx.text(-0.045, len(M) - 0.35 + 0.72,
        fig_annot.record('Fig_toxicity', 'panel_b_header',
            f'no regulator-grade acute inhalation value exists for any of these {len(M)}'),
        fontsize=FS_AXIS, color=C_INK, va='center', ha='left')
dress(bx); bx.grid(False); panel_label(bx, 'b')

fig.savefig('outputs/Fig_toxicity.png', dpi=200, bbox_inches='tight', facecolor='white')
print('[fig] outputs/Fig_toxicity.png')
print(f'  (a) {len(H)}종  (b) {len(M)}종')
