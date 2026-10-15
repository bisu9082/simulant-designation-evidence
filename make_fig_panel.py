#!/usr/bin/env python3
"""Figure: 패널 구성이 결론을 정한다 (2026-10-14 신설).

35종을 주 패널로 올리면서, 그 전환이 무엇을 바꾸는지 한 장에 담는다.
(a) 작용제별 최근접 거리. 기존 모사체까지의 거리는 패널을 바꿔도 거의 그대로다.
    바뀌는 것은 '더 가까운 것이 생긴다' 는 사실뿐이다. 이 구분이 논문의 핵심이라
    두 거리를 같은 행에 나란히 놓는다.
(b) 보정 4방식 x 패널 3구성의 최소 조정값. 0.05 를 넘는지 한눈에 보이게 한다.

스타일은 make_fig_main.py 의 정본 설정을 그대로 따른다.
"""
import json

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import fig_annot

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 28, 16, 15, 13, 13, 11
C_REAL, C_SIM, C_INTF, C_REF, C_INK = ('#C94F4A', '#3A6EA8', '#E8943A',
                                       '#B84040', '#333333')
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.facecolor': 'white', 'figure.facecolor': 'white'})
SHORT = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
         'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD',
         'lewisite 1 (L)': 'L', 'hydrogen cyanide (AC)': 'AC',
         'cyanogen chloride (CK)': 'CK', 'phosgene (CG)': 'CG',
         'chloropicrin (PS)': 'PS', 'dichloromethane': 'DCM'}
s = lambda x: SHORT.get(x, x)


def panel_label(ax, L, dy=6):
    ax.annotate(f'({L})', xy=(0.0, 1.0), xycoords='axes fraction',
                xytext=(0, dy), textcoords='offset points', fontsize=FS_LABEL,
                fontweight='bold', va='bottom', ha='left', annotation_clip=False)


def dress(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(color='#DDDDDD', ls='--', lw=0.6)
    ax.set_axisbelow(True)


M = json.load(open('outputs/match_by_panel.json'))
S = json.load(open('outputs/stat_supplement.json'))
B = json.load(open('outputs/stat_supplement2.json'))

inc = {a['agent']: a for a in M['panel35_incumbent3']['agents']}
allc = {a['agent']: a for a in M['panel35_incumbent_plus_lit8']['agents']}
order = sorted(inc, key=lambda a: allc[a]['d_norm'])

fig, (ax, cx, bx) = plt.subplots(1, 3, figsize=(27, 9),
                                 gridspec_kw={'width_ratios': [1.35, 0.85, 1.0]})
fig.subplots_adjust(left=0.058, right=0.988, top=0.90, bottom=0.115, wspace=0.30)

# ---- (a) 작용제별 최근접 거리: 기존 3종 대 전체 후보
y = np.arange(len(order))[::-1]
# 두 마커가 같은 값일 때 겹쳐 하나로 보이므로 행 안에서 위아래로 떼어 놓는다.
DY = 0.17
for yy, a in zip(y, order):
    d_inc, d_all = inc[a]['d_norm'], allc[a]['d_norm']
    ax.plot([d_all, d_inc], [yy - DY, yy + DY], color='#BBBBBB', lw=2.2,
            zorder=1)
    ax.scatter([d_inc], [yy + DY], s=230, c=C_SIM, edgecolors='white', lw=1.8,
               zorder=3)
    ax.scatter([d_all], [yy - DY], s=230, c=C_REAL, edgecolors='white', lw=1.8,
               marker='D', zorder=3)
    ax.text(d_inc + 0.085, yy + DY, s(inc[a]['nearest_simulant']),
            fontsize=FS_ANNOT, va='center', ha='left', color=C_SIM)
    ax.text(d_all - 0.085, yy - DY, s(allc[a]['nearest_simulant']),
            fontsize=FS_ANNOT, va='center', ha='right', color=C_REAL)
ax.axvline(2.0, color=C_REF, ls='--', lw=2.0, zorder=2)
ax.text(2.0, -1.05, ' match threshold', color=C_REF,
        fontsize=FS_ANNOT, va='bottom', ha='left')
ax.set_yticks(y)
ax.set_yticklabels([s(a) for a in order])
ax.set_xlim(0.18, 3.95)
ax.set_ylim(-1.15, len(order) + 0.95)
ax.set_xlabel('distance to nearest candidate, as a multiple of the\n'
              'typical nearest-neighbour separation of the merged space')
ax.legend(handles=[
    Line2D([], [], marker='o', ls='', ms=14, mfc=C_SIM, mec='white', mew=1.8,
           label='nearest of the three incumbent simulants'),
    Line2D([], [], marker='D', ls='', ms=13, mfc=C_REAL, mec='white', mew=1.8,
           label="nearest once the field's five other simulants are included")],
    loc='upper right', frameon=False, bbox_to_anchor=(1.03, 1.012),
    handletextpad=0.5, labelspacing=0.5)
dress(ax)
panel_label(ax, 'a')

# ---- (c) 서술자 부트스트랩: 유효매칭 수 분포
BT = B['distance_bootstrap']
KK = [('incumbent3', 'the three\nincumbent simulants', C_SIM),
      ('incumbent_plus_lit8', "the field's\neight simulants", C_REAL)]
width = 0.38
for i, (key, lab, col) in enumerate(KK):
    o = BT[key]
    cx.barh([0 + (i - 0.5) * 0.9], [o['median']], height=0.0)   # 자리만
positions = np.arange(len(KK))[::-1]
for pos, (key, lab, col) in zip(positions, KK):
    o = BT[key]
    lo, hi = o['ci']
    cx.plot([lo, hi], [pos, pos], color=col, lw=7, solid_capstyle='butt',
            alpha=0.33, zorder=2)
    cx.scatter([o['mean']], [pos], s=320, c=col, edgecolors='white', lw=2.0,
               zorder=4)
    cx.scatter([o['median']], [pos], s=170, marker='|', c=C_INK, lw=3.0,
               zorder=5)
    cx.text(hi + 0.25, pos, f"median {o['median']}\n95% CI {lo} to {hi}",
            fontsize=FS_ANNOT, va='center', ha='left', color=col)
cx.set_yticks(positions)
cx.set_yticklabels([l for _, l, _ in KK])
cx.set_xlim(-0.5, 11.3)
cx.set_ylim(-1.25, len(KK) - 0.35)
cx.set_xlabel('agents of eleven with a valid match,\nover 1000 descriptor'
              ' bootstrap replicates')
d = BT['difference']
cx.text(0.5, -0.95,
        f"difference positive in {100*d['frac_positive']:.0f}% of replicates\n"
        f"(median {d['median']}, 95% CI {d['ci'][0]} to {d['ci'][1]})",
        fontsize=FS_ANNOT, va='center', ha='left', color=C_REF,
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.35', alpha=0.92))
dress(cx)
panel_label(cx, 'b')

# ---- (b) 보정 4방식 x 3구성
CFG = [('panel28_incumbent3', '28 cpd.\nincumbent 3'),
       ('panel35_incumbent3', '35 cpd.\nincumbent 3'),
       ('panel35_incumbent_plus_lit8', '35 cpd.\n+ literature 8')]
SCH = [('min_p', 'smallest nominal $p$', 'o'),
       ('bh_all_min_q', 'BH over 18', 's'),
       ('bh_by_space_min_q', 'BH over 6, per space', '^'),
       ('wy', 'Westfall-Young min$P$', 'v'),
       ('prespecified_bh_min_q', 'statistic fixed in advance', 'D')]
COLS = [C_INK, '#7F7F7F', '#A9A9A9', '#6B8E23', C_REAL]
WY = {k: v['wy_adjusted_p'] for k, v in B['westfall_young'].items()}
xs = np.arange(len(CFG))
for (key, lab, mk), col, off in zip(SCH, COLS,
                                    np.linspace(-0.26, 0.26, len(SCH))):
    v = [(WY.get(c, np.nan) if key == 'wy' else S['family'][c][key])
         for c, _ in CFG]
    bx.scatter(xs + off, v, s=230, marker=mk, c=col, edgecolors='white',
               lw=1.6, zorder=3, label=lab)
bx.axhline(0.05, color=C_REF, ls='--', lw=2.0, zorder=2)
bx.text(2.48, 0.052, r'$\alpha$ = 0.05', color=C_REF, fontsize=FS_ANNOT,
        va='bottom', ha='right')
bx.set_yscale('log')
bx.set_ylim(0.02, 1.5)
bx.set_xticks(xs)
bx.set_xticklabels([l for _, l in CFG])
bx.set_xlim(-0.55, 2.55)
bx.set_ylabel('smallest adjusted value over the 18 tests')
bx.legend(loc='lower right', frameon=False, ncol=1,
          bbox_to_anchor=(1.03, -0.02), handletextpad=0.5, labelspacing=0.45)
dress(bx)
panel_label(bx, 'c')

n0 = M['panel35_incumbent3']['n_match']
n8 = M['panel35_incumbent_plus_lit8']['n_match']
n28 = M['panel28_incumbent3']['n_match']
fig_annot.record('Fig_panel', 'panel_a',
                 f'valid matches: {n28} of 11 on the 28-compound panel, '
                 f'{n0} of 11 on the 35-compound panel with the same three '
                 f'simulants, {n8} of 11 once the five literature simulants '
                 f'are admitted as simulants.')
fig_annot.record('Fig_panel', 'panel_b',
                 f"incumbent three match a median of "
                 f"{BT['incumbent3']['median']} of 11 agents "
                 f"(95% CI {BT['incumbent3']['ci'][0]} to "
                 f"{BT['incumbent3']['ci'][1]}); the field's eight match a "
                 f"median of {BT['incumbent_plus_lit8']['median']} "
                 f"(95% CI {BT['incumbent_plus_lit8']['ci'][0]} to "
                 f"{BT['incumbent_plus_lit8']['ci'][1]}).")
fig_annot.record('Fig_panel', 'panel_c',
                 'only the pre-specified-statistic scheme on the '
                 '28-compound panel falls below 0.05.')

fig.savefig('outputs/Fig_panel.png', dpi=200, bbox_inches='tight',
            facecolor='white')
print('[fig] outputs/Fig_panel.png')
print(f'  유효매칭 28종 {n28}/11, 35종 지정3 {n0}/11, 35종 지정8 {n8}/11')
