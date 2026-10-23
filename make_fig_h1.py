#!/usr/bin/env python3
"""H1 결합해제 검정 + 구조 경보 게이트 검증 그림"""
import json
import numpy as np
import pandas as pd
import matplotlib

import fig_annot
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT = 19, 15, 14, 12, 12
C_REAL, C_SIM, C_INTF, C_REF, C_INK = '#C94F4A', '#3A6EA8', '#E8943A', '#B84040', '#333333'
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.facecolor': 'white', 'figure.facecolor': 'white'})


def panel_label(ax, L):
    ax.annotate(f'({L})', xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, 6),
                textcoords='offset points', fontsize=FS_LABEL, fontweight='bold',
                va='bottom', ha='left', annotation_clip=False)


def dress(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(color='#DDDDDD', ls='--', lw=0.6); ax.set_axisbelow(True)


S = json.load(open('outputs/h1_summary.json'))
D = pd.read_csv('outputs/h1_candidates.csv')
fid, tox = S['fidelity_cut'], S['tox_cut']

# figsize 는 \textwidth(384 pt = 5.33 in) 기준으로 정한다. 24~27 in 로
# 그리면 축소율이 0.2 가 되어 15 pt 틱이 실효 3 pt 로 인쇄된다 (2026-10-20
# 외부 검토 지적, 전 도판 해당). 10 in 안쪽으로 그려 실효 7 pt 이상을 둔다.
# 패널을 세로로 쌓는다. 가로 2열이면 (a) 가 3.1 in 밖에 안 되어 점별 라벨이
# 겹친다. 세로로 쌓으면 (a) 가 전폭(5.33 in)을 받는다.
fig, (ax, bx) = plt.subplots(2, 1, figsize=(10.0, 12.2),
                             gridspec_kw=dict(height_ratios=[1.0, 1.0]))
fig.subplots_adjust(left=0.085, right=0.985, top=0.972, bottom=0.075,
                    hspace=0.52)

# ---------------- (a) H1 산점도
clean = D[~D.alert_flag]
flag = D[D.alert_flag]
from matplotlib.patches import Rectangle
RELAX = tox - 1.0                     # 현행 모사체 DCM 의 1/10 안전 기준
ax.add_patch(Rectangle((-0.12, RELAX), fid + 0.12, 50, fc='#EDF3E9', ec='none',
                       zorder=0))
ax.axvline(fid, color=C_INK, ls=(0, (5, 3)), lw=1.2, zorder=1)
ax.axhline(tox, color=C_INK, ls=(0, (5, 3)), lw=1.2, zorder=1)
ax.axhline(RELAX, color='#4E7A3A', ls=(0, (2, 2)), lw=1.2, zorder=1)
ax.scatter(flag.fidelity_dist, flag.y_tox, s=24, marker='x', c=C_REF, lw=1.5,
           zorder=3, label='structural alert, Tier A/B (excluded)')
ax.scatter(clean.fidelity_dist, clean.y_tox, s=25, c=C_INTF, edgecolors='white',
           lw=0.8, alpha=0.92, zorder=4, label='candidate, passes alert gate')
own = D[D.src == 'study28']
ax.scatter(own.fidelity_dist, own.y_tox, s=59, c=C_SIM, edgecolors='white',
           lw=1.3, zorder=6, label='compound already in this study')
t = [ax.text(r.fidelity_dist, r.y_tox, r.label, fontsize=FS_ANNOT - 1, color=C_INK,
             zorder=7) for r in own.itertuples()]
best = clean[(clean.fidelity_dist <= fid) & (clean.y_tox >= RELAX)] \
    .sort_values('fidelity_dist')
ax.scatter(best.fidelity_dist, best.y_tox, s=56, facecolors='none',
           edgecolors='#4E7A3A', lw=1.7, zorder=5,
           label=f'clears the relaxed bar ({len(best)})')
t += [ax.text(r.fidelity_dist, r.y_tox, r.label.split('(')[0].strip()[:26],
              fontsize=FS_ANNOT - 2, color='#3C5E2E', fontweight='bold', zorder=7)
      for r in best.itertuples()]
XMAX = float(np.percentile(D.fidelity_dist, 97)) * 1.12      # 이상치 클립
n_clip = int((D.fidelity_dist > XMAX).sum())
ax.set_xlim(-0.12, XMAX)
ax.set_ylim(float(D.y_tox.min()) - 0.45, float(D.y_tox.max()) + 2.0)
adjust_text(t, x=D.fidelity_dist.to_numpy(), y=D.y_tox.to_numpy(), ax=ax,
            only_move={'text': 'xy'}, expand=(1.9, 2.2), force_text=(1.0, 1.4),
            force_static=(1.1, 1.5), max_move=90, iter_lim=900,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
YTOP = float(D.y_tox.max()) + 2.0
n_strict = int(((clean.fidelity_dist <= fid) & (clean.y_tox >= tox)).sum())
ax.text(fid * 0.5, YTOP - 0.22, 'FAITHFUL  +  SAFE', fontsize=FS_ANNOT + 2,
        fontweight='bold', color=C_REF, ha='center', va='center')
# 기준선 라벨은 데이터가 비어 있는 오른쪽 끝에 둔다
ax.text(0.995, tox + 0.07, 'dichloromethane, the incumbent simulant',
        transform=ax.get_yaxis_transform(), fontsize=FS_ANNOT - 2, color=C_INK,
        ha='right', va='bottom')
ax.text(0.995, RELAX + 0.07, 'one tenth as safe', fontsize=FS_ANNOT - 2,
        transform=ax.get_yaxis_transform(), color='#3C5E2E', ha='right',
        va='bottom')

ax.set_xlabel('fidelity: standardised electronic-space distance to the\n'
              'nearest authentic agent   (smaller = more faithful)')
ax.set_ylabel(r'safety: measured log$_{10}$ AEGL-3 10 min (mg m$^{-3}$)'
              '\n(larger = safer)')
ax.text(1.0, -0.300, transform=ax.transAxes, s=
        f"Spearman $r$ = {S['spearman_r']:+.3f}   "
        f"95% CI [{S['spearman_ci'][0]:+.3f}, {S['spearman_ci'][1]:+.3f}]   "
        f"$n$ = {S['n_candidates']}\n"
        f"no monotonic association: fidelity explains "
        f"{100*S['r_squared']:.1f}% of the variance in toxicity\n"
        f"faithful and safe at the incumbent bar: {n_strict}.   "
        f"One order of magnitude below it: {len(best)}, all clearing the alert gate\n"
        f"({n_clip} far-field point(s) beyond x = {XMAX:.1f} not shown)",
        fontsize=FS_ANNOT, color=C_INK, ha='right', va='top',
        bbox=dict(fc='white', ec=C_INK, boxstyle='round,pad=0.4', alpha=0.93))
ax.legend(loc='upper right', frameon=False, bbox_to_anchor=(1.0, 1.0))
dress(ax); panel_label(ax, 'a')

# ---------------- (b) 게이트 검증
G = pd.read_csv('descriptor_cache_xtb.csv')[['compound_name', 'compound_class']]
P = D[D.src == 'study28'][['label', 'tier']].rename(columns={'label': 'compound_name'})
ALLT = pd.read_csv('outputs/h1_candidates.csv')
tiers = {}
for _, r in pd.read_csv('outputs/h1_candidates.csv').iterrows():
    tiers[r.label] = r.tier
TIER = json.load(open('outputs/h1_gate_tiers.json'))
G['tier'] = G.compound_name.map(TIER).fillna('-')
ORD = {'real_agent': 0, 'simulant': 1, 'interferent': 2}
G = G.sort_values(['compound_class', 'tier'],
                  key=lambda c: c.map(ORD) if c.name == 'compound_class' else c)
TCOL = {'A': C_REF, 'B': '#E8943A', 'C': '#7FAF6A', '-': '#CCCCCC'}
y = np.arange(len(G))[::-1]
for yy, r in zip(y, G.itertuples()):
    bx.barh(yy, 1.0, color=TCOL[r.tier], edgecolor='white', height=0.78, zorder=2)
    bx.text(1.06, yy, {'A': 'Tier A  excluded', 'B': 'Tier B  excluded',
                       'C': 'Tier C  noted only', '-': 'no alert'}[r.tier],
            fontsize=FS_ANNOT - 4, va='center', ha='left', color='#555555')
bx.set_yticks(y[::-1])
bx.set_yticklabels([f'{n}' for n in G.compound_name][::-1], fontsize=FS_TICK - 5)
for yy, r in zip(y, G.itertuples()):
    bx.text(-0.04, yy, {'real_agent': 'agent', 'simulant': 'simulant',
                        'interferent': ''}[r.compound_class],
            fontsize=FS_ANNOT - 3, va='center', ha='right',
            color={'real_agent': C_REAL, 'simulant': C_SIM,
                   'interferent': '#999999'}[r.compound_class])
bx.set_xlim(-0.5, 2.6); bx.set_ylim(-6.6, len(G) - 0.4 + 1.6)
bx.set_xticks([]); bx.spines['bottom'].set_visible(False)
bx.spines['left'].set_bounds(0, len(G) - 1)
bx.text(-0.5, -1.5, '9 of 11 authentic agents Tier A, phosgene Tier B;\n'
        'chloropicrin not caught (no agent-class motif).\n'
        'Zero false alerts on 14 interferents. DMMP is Tier C: it carries the\n'
        'phosphoryl ester but no leaving group, so it is noted, not excluded.',
        fontsize=FS_ANNOT - 2, color=C_INK, va='top', ha='left',
        bbox=dict(fc='white', ec=C_INK, boxstyle='round,pad=0.4', alpha=0.93))
dress(bx); bx.grid(False); panel_label(bx, 'b')

fig_annot.record('Fig_h1', 'spearman',
                 f"Spearman r = {S['spearman_r']:+.3f}, 95% CI "
                 f"{S['spearman_ci'][0]:+.3f} to {S['spearman_ci'][1]:+.3f}, "
                 f"n = {S['n_candidates']}; no monotonic association.")
fig_annot.record('Fig_h1', 'quadrant',
                 f'faithful and safe at the incumbent bar: {n_strict}; '
                 f'one order of magnitude below it: {len(best)}.')

fig.savefig('outputs/Fig_h1.png', dpi=200, bbox_inches='tight', facecolor='white')
print('[fig] outputs/Fig_h1.png')
