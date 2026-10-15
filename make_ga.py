#!/usr/bin/env python3
"""simulant / Graphical abstract — Elsevier 가로형. 실측 데이터만 [C-3].
A·B = C(17,3)=680 전수 순열 귀무분포 + 실제 지정의 관측값
C   = 충실도 x 안전성 평면 (n=140, 실측 양축)"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

C_REAL, C_SIM, C_INTF, C_REF, C_INK = '#C94F4A', '#3A6EA8', '#E8943A', '#B84040', '#333333'
C_OK, C_NO, C_NULL = '#4E7A3A', '#B84040', '#C8C8C8'
plt.rcParams.update({'font.family': 'sans-serif', 'axes.facecolor': 'white',
                     'figure.facecolor': 'white'})

B = json.load(open('outputs/sanity_checks.json'))['sanity_baseline']
H = json.load(open('outputs/h1_summary.json'))
D = pd.read_csv('outputs/h1_candidates.csv')

fig = plt.figure(figsize=(17.5, 7.0))
gs = fig.add_gridspec(1, 3, left=0.052, right=0.988, top=0.700, bottom=0.115,
                      wspace=0.255, width_ratios=[1.0, 1.0, 1.30])
ax, bx, cx = (fig.add_subplot(g) for g in gs)

FS_T, FS_AX, FS_TK, FS_N, FS_V = 21, 16, 14, 15, 26


def dress(a):
    a.spines[['top', 'right']].set_visible(False)
    a.grid(color='#DDDDDD', ls='--', lw=0.7, axis='y')
    a.set_axisbelow(True)
    a.tick_params(labelsize=FS_TK)


def null_panel(a, cond, title, verdict, ok, tx):
    null = np.load(f'outputs/_null_{cond}.npy')
    obs = B[cond]['observed_diff']
    p = B[cond]['p_exact']
    col = C_OK if ok else C_NO
    a.hist(null, bins=34, color=C_NULL, edgecolor='white', lw=0.8, zorder=2)
    ymax = a.get_ylim()[1]
    a.set_ylim(0, ymax * 1.22)
    a.axvline(obs, color=col, lw=4.0, zorder=5)
    a.annotate('the actual\nthree simulants',
               xy=(obs, ymax * 0.42), xycoords='data',
               xytext=(tx, 0.55), textcoords='axes fraction',
               fontsize=FS_N - 1, color=col, fontweight='bold',
               ha='left', va='center',
               arrowprops=dict(arrowstyle='-|>', color=col, lw=2.4,
                               shrinkA=4, shrinkB=6))
    a.set_title(title, fontsize=FS_T, fontweight='bold', color=C_INK, pad=44)
    a.set_xlabel('mean  d(nearest simulant) $-$ d(nearest interferent)\n'
                 'over all 680 ways of picking three non-agents',
                 fontsize=FS_AX - 3, labelpad=8)
    a.set_ylabel('number of designations', fontsize=FS_AX - 2)
    a.text(0.5, 1.015, verdict, transform=a.transAxes, fontsize=FS_V,
           fontweight='bold', color=col, ha='center', va='bottom',
           clip_on=False)
    a.text(0.62, 0.945, f'exact $p$ = {p:.3f}', transform=a.transAxes,
           fontsize=FS_N + 2, color=C_INK, ha='center', va='center')
    dress(a)


null_panel(ax, 'A_rdkit', 'judged in structural space',
           'beats chance', True, 0.26)
null_panel(bx, 'B_electronic', 'judged in electronic space',
           'no better than chance', False, 0.52)

# ---------------- (c) 두 축의 독립성
fid, tox = H['fidelity_cut'], H['tox_cut']
clean, flag = D[~D.alert_flag], D[D.alert_flag]
XMAX = float(np.percentile(D.fidelity_dist, 97)) * 1.10
YLO, YHI = float(D.y_tox.min()) - 0.45, float(D.y_tox.max()) + 4.40
XMIN = -0.55
cx.add_patch(Rectangle((XMIN, tox), fid - XMIN, YHI - tox, fc='#E9F0E4',
                       ec='none', zorder=0))
cx.axvline(fid, color=C_INK, ls=(0, (5, 3)), lw=2.0, zorder=1)
cx.axhline(tox, color=C_INK, ls=(0, (5, 3)), lw=2.0, zorder=1)
cx.scatter(flag.fidelity_dist, flag.y_tox, s=80, marker='x', c=C_REF, lw=1.9,
           zorder=3, label='excluded by structural-alert gate')
cx.scatter(clean.fidelity_dist, clean.y_tox, s=100, c=C_INTF,
           edgecolors='white', lw=1.1, alpha=0.9, zorder=4,
           label='candidate with both axes measured')
own = D[D.src == 'study28']
cx.scatter(own.fidelity_dist, own.y_tox, s=250, c=C_SIM, edgecolors='white',
           lw=2.0, zorder=6, label='compound already in this study')
cx.set_xlim(XMIN, XMAX)
cx.set_ylim(YLO, YHI)
cx.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
cx.invert_xaxis()
cx.set_title('the two axes are independent', fontsize=FS_T, fontweight='bold',
             color=C_INK, pad=44)
cx.set_xlabel('fidelity: electronic distance to the nearest agent\n'
              r'$\longrightarrow$  more faithful', fontsize=FS_AX - 3, labelpad=8)
cx.set_ylabel('safety: measured log$_{10}$ AEGL-3\n'
              r'$\longrightarrow$  safer', fontsize=FS_AX - 3)
cx.text(fid * 0.99, tox + 0.52 * (YHI - tox), 'FAITHFUL\n+ SAFE',
        fontsize=FS_N - 2, fontweight='bold', color=C_OK, ha='left', va='center')
RELAX = tox - 1.0
_clean = D[~D.alert_flag]
_best = _clean[(_clean.fidelity_dist <= fid) & (_clean.y_tox >= RELAX)]
_strict = int(((_clean.fidelity_dist <= fid) & (_clean.y_tox >= tox)).sum())
cx.axhline(RELAX, color='#4E7A3A', ls=(0, (2, 2)), lw=1.8, zorder=1)
cx.scatter(_best.fidelity_dist, _best.y_tox, s=170, facecolors='none',
           edgecolors='#4E7A3A', lw=2.4, zorder=7,
           label=f'clears the relaxed bar ({len(_best)})')

cx.text(0.030, 0.715,
        f"Spearman $r$ = {H['spearman_r']:+.3f}   "
        f"95% CI [{H['spearman_ci'][0]:+.3f}, {H['spearman_ci'][1]:+.3f}]   "
        f"$n$ = {H['n_candidates']}\n"
        f"faithful + safe: {_strict} at the incumbent simulant's own bar,\n"
        f"{len(_best)} one order of magnitude below it",
        transform=cx.transAxes, fontsize=FS_N - 4, color=C_INK, ha='left',
        va='top', zorder=12,
        bbox=dict(fc='white', ec=C_INK, boxstyle='round,pad=0.40', alpha=1.0))
cx.legend(loc='upper left', bbox_to_anchor=(0.012, 0.985), frameon=False,
          fontsize=FS_N - 5, handletextpad=0.45, labelspacing=0.38)
dress(cx)
cx.grid(color='#DDDDDD', ls='--', lw=0.7)

# ---------------- 머리글
fig.text(0.5, 0.945,
         'Simulant validity is a property of the descriptor space, not of the compound',
         fontsize=28, fontweight='bold', color=C_INK, ha='center', va='center')
fig.text(0.5, 0.865,
         'The same eleven agents and three designated simulants, scored two ways.  '
         'Safety is a third, unlinked axis.',
         fontsize=17, color='#555555', ha='center', va='center')

fig.savefig('outputs/graphical_abstract.png', dpi=190, bbox_inches='tight',
            facecolor='white', pad_inches=0.20)
fig.savefig('outputs/graphical_abstract.pdf', bbox_inches='tight',
            facecolor='white', pad_inches=0.20)
print('[GA] outputs/graphical_abstract.png / .pdf')
