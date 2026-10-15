#!/usr/bin/env python3
"""simulant / STEP 4 — 본 Figure (6패널). 정본 라벨(모사체 3종) 기준."""
import json
import numpy as np
import pandas as pd

import fig_annot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from adjustText import adjust_text

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 28, 16, 15, 13, 13, 11
C_REAL, C_SIM, C_INTF, C_REF, C_INK = '#C94F4A', '#3A6EA8', '#E8943A', '#B84040', '#333333'
CMAP = {'real_agent': C_REAL, 'simulant': C_SIM, 'interferent': C_INTF}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.facecolor': 'white', 'figure.facecolor': 'white'})
SHORT = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
         'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD', 'lewisite 1 (L)': 'L',
         'hydrogen cyanide (AC)': 'AC', 'cyanogen chloride (CK)': 'CK',
         'phosgene (CG)': 'CG', 'chloropicrin (PS)': 'PS',
         'dichloromethane': 'DCM', 'p-dichlorobenzene': 'p-DCB',
         'ethyl acetate': 'EtOAc', 'acetonitrile': 'ACN', 'n-hexane': 'hexane'}
s = lambda x: SHORT.get(x, x)
SPACE_LBL = {'A_rdkit': 'structural / physicochemical',
             'B_electronic': 'GFN2 electronic', 'C_merged': 'merged'}
MARK = {'A_rdkit': 'o', 'B_electronic': 's', 'C_merged': '^'}


def panel_label(ax, L, fs=FS_LABEL, dy=6):
    ax.annotate(f'({L})', xy=(0.0, 1.0), xycoords='axes fraction', xytext=(0, dy),
                textcoords='offset points', fontsize=fs, fontweight='bold',
                va='bottom', ha='left', annotation_clip=False)


def dress(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(color='#DDDDDD', ls='--', lw=0.6); ax.set_axisbelow(True)


df = pd.read_csv('descriptor_cache_xtb.csv')
S = json.load(open('outputs/experiment_summary.json'))
R = pd.read_csv('outputs/distance_table_normalized.csv')
CON = {c['condition']: c for c in S['conditions']}


def scores(cond):
    keep = CON[cond]['kept']
    X = df[keep].to_numpy(float)
    Xz = (X - X.mean(0)) / X.std(0, ddof=1)
    p = PCA(n_components=2, random_state=42).fit(Xz)
    return p.transform(Xz), p.explained_variance_ratio_, keep, p.components_


fig, axes = plt.subplots(2, 3, figsize=(27, 16))
fig.subplots_adjust(left=0.055, right=0.985, top=0.925, bottom=0.06,
                    hspace=0.34, wspace=0.26)
ax, bx, cx = axes[0]
dx, ex, fx = axes[1]


def pca_panel(axis, cond, letter, xlab):
    Y, evr, keep, comp = scores(cond)
    for k in ['interferent', 'simulant', 'real_agent']:
        m = (df.compound_class == k).to_numpy()
        axis.scatter(Y[m, 0], Y[m, 1], s=175, c=CMAP[k], edgecolors='white',
                     lw=1.6, zorder=3,
                     label={'real_agent': 'authentic agent',
                            'simulant': 'incumbent simulant',
                            'interferent': 'interferent / solvent'}[k])
    t = [axis.text(Y[i, 0], Y[i, 1], s(df.compound_name[i]), fontsize=FS_ANNOT - 2,
                   color=C_INK, zorder=4) for i in range(len(df))]
    axis.set_xlim(Y[:, 0].min() - 0.9, Y[:, 0].max() + 0.9)
    axis.set_ylim(Y[:, 1].min() - 1.0, Y[:, 1].max() + 1.2)
    adjust_text(t, x=Y[:, 0], y=Y[:, 1], ax=axis, only_move={'text': 'xy'},
                expand=(2.1, 2.5), force_text=(1.1, 1.5), force_static=(1.3, 1.7),
                force_pull=(0.004, 0.004), max_move=120, iter_lim=1500,
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7))
    axis.set_xlabel(f'PC1 ({evr[0]*100:.1f}%): {xlab}')
    axis.set_ylabel(f'PC2 ({evr[1]*100:.1f}%)')
    dress(axis); panel_label(axis, letter)   # 색 범례는 패널 (f)에 1회만 표시
    return evr


pca_panel(ax, 'A_rdkit', 'a', 'molecular size')
pca_panel(bx, 'B_electronic', 'b', 'frontier orbital energies')

# ---- (c) 정규화 최근접-모사체 거리, 3개 공간
order = (R[R.condition == 'C_merged'].sort_values('d_sim_norm').real_agent.tolist())
yy = {n: len(order) - 1 - i for i, n in enumerate(order)}
for cond in ['A_rdkit', 'B_electronic', 'C_merged']:
    t = R[R.condition == cond]
    for r in t.itertuples():
        ok = (not r.interferent_closer) and r.d_sim_norm <= 2.0
        cx.scatter(r.d_sim_norm, yy[r.real_agent], marker=MARK[cond], s=145,
                   facecolors=(C_SIM if ok else 'white'), edgecolors=C_INK,
                   lw=1.5, zorder=3)
cx.axvline(2.0, color=C_REF, ls=(0, (5, 3)), lw=2.0, zorder=1)
cx.set_yticks(range(len(order)))
cx.set_yticklabels([s(n) for n in order][::-1])
xm = float(R.d_sim_norm.max())
cx.set_xlim(0, xm * 1.16); cx.set_ylim(-4.6, len(order) - 0.45 + 1.15)
cx.spines['left'].set_bounds(0, len(order) - 1)
cx.set_xlabel('nearest incumbent simulant, distance / typical\nnearest-neighbour '
              'distance of that space')
cx.text(2.08, len(order) - 0.6, 'match\nthreshold', fontsize=FS_ANNOT, color=C_REF,
        va='top', ha='left')
cx.legend(handles=[Line2D([], [], marker=MARK[c], ls='', ms=11, mfc='white',
                          mec=C_INK, mew=1.5, label=SPACE_LBL[c])
                   for c in MARK] +
                  [Line2D([], [], marker='o', ls='', ms=11, mfc=C_SIM, mec=C_INK,
                          mew=1.5, label='valid match (also beats every interferent)')],
          loc='lower right', frameon=False, bbox_to_anchor=(1.03, -0.015),
          handletextpad=0.5, labelspacing=0.55)
dress(cx); panel_label(cx, 'c')

# ---- (d) 전수 순열 귀무분포 vs 실제 지정 모사체
SAN = json.load(open('outputs/sanity_checks.json'))['sanity_baseline']
cond_order = ['A_rdkit', 'B_electronic', 'C_merged']
yv = np.arange(len(cond_order))[::-1]
rngj = np.random.default_rng(7)
for y, c in zip(yv, cond_order):
    null = np.load(f'outputs/_null_{c}.npy')
    dx.scatter(null, y + rngj.uniform(-0.17, 0.17, null.size), s=9,
               c='#BBBBBB', alpha=0.55, lw=0, zorder=2)
    dx.scatter([null.mean()], [y], marker='|', s=900, c=C_INK, lw=2.4, zorder=3)
    k = SAN[c]
    dx.scatter([k['observed_diff']], [y], s=300, c=C_REAL, edgecolors='white',
               lw=2.2, zorder=5)
    dx.text(k['observed_diff'], y + 0.40,
            f"exact p = {k['p_exact']:.3f}", fontsize=FS_ANNOT, va='bottom',
            ha='center', color=(C_REF if k['p_exact'] < 0.05 else C_INK),
            fontweight='bold' if k['p_exact'] < 0.05 else 'normal')
# 패널 (d) 에 렌더링되는 p 를 검증 대상에 등재한다. 2026-10-13. 이 세 수치가
# 도판 안에만 있어서 본문이 0.275 로 바뀐 뒤에도 그림은 0.297 을 그리고 있었다.
for c in cond_order:
    fig_annot.record('Fig_main', f'panel_d_{c}',
                     f"exact p = {SAN[c]['p_exact']:.3f}")
dx.set_yticks(yv); dx.set_yticklabels([SPACE_LBL[c] for c in cond_order])
dx.set_xlim(-1.5, 4.3); dx.set_ylim(-1.35, len(cond_order) - 1 + 0.95)
dx.spines['left'].set_bounds(yv.min(), yv.max())
dx.set_xlabel('mean ( d to nearest simulant  \u2212  d to nearest interferent )\n'
              'grey = all 680 ways of designating 3 of the 17 non-agents')
dx.legend(handles=[Line2D([], [], marker='o', ls='', ms=9, mfc='#BBBBBB',
                          mec='none', label='null: 3 compounds picked at random'),
                   Line2D([], [], marker='|', ls='', ms=16, mec=C_INK, mew=2.4,
                          label='null mean'),
                   Line2D([], [], marker='o', ls='', ms=13, mfc=C_REAL,
                          mec='white', mew=2, label='actual incumbent simulants')],
          loc='lower right', frameon=False, bbox_to_anchor=(1.03, -0.015),
          handletextpad=0.6, labelspacing=0.5)
dress(dx); panel_label(dx, 'd')

# ---- (e) 인 중심 국소 반응성
P = df[df.q_P.notna()].reset_index(drop=True)
g3 = P[P.compound_name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])]
gm, gs = g3.fukui_f_plus_P.mean(), g3.fukui_f_plus_P.std(ddof=1)
ex.axhspan(gm - gs, gm + gs, color=C_REAL, alpha=0.13, zorder=0)
ex.axhline(gm, color=C_REF, ls=(0, (5, 3)), lw=1.8, zorder=1)
for r in P.itertuples():
    ex.scatter(r.q_P, r.fukui_f_plus_P, s=290, c=CMAP[r.compound_class],
               edgecolors='white', lw=2.0, zorder=3)
OFF = {'sarin (GB)': (14, 7), 'soman (GD)': (16, -3), 'cyclosarin (GF)': (14, -19),
       'DMMP': (16, 0), 'tabun (GA)': (16, 2), 'VX': (16, 0)}
for r in P.itertuples():
    ex.annotate(s(r.compound_name), xy=(r.q_P, r.fukui_f_plus_P),
                xytext=OFF[r.compound_name], textcoords='offset points',
                fontsize=FS_ANNOT + 1, color=C_INK, va='center', ha='left',
                fontweight='bold' if r.compound_class == 'simulant' else 'normal')
ex.set_xlim(0.28, 0.62); ex.set_ylim(0.10, 0.89)
ex.text(0.292, gm + gs + 0.018,
        f'GB/GD/GF mean $\\pm$ SD  (${gm:.3f} \\pm {gs:.3f}$)',
        fontsize=FS_ANNOT, color=C_REF, va='bottom', ha='left',
        bbox=dict(fc='white', ec='none', alpha=0.9))
ex.annotate('', xy=(0.300, 0.2043), xytext=(0.300, 0.5585),
            arrowprops=dict(arrowstyle='<->', color=C_INK, lw=1.6))
# 수치를 데이터에서 계산한다 (2026-10-12). 종전에는 2.7 이 문자열 리터럴이었다.
_f = dict(zip(P6.compound_name, P6.fukui_f_plus_P)) if 'P6' in dir() else \
     dict(zip(df.compound_name, df.fukui_f_plus_P))
_ratio = _f['DMMP'] / _f['VX']
ex.text(0.392, 0.300,
        fig_annot.record('Fig_main', 'panel_e',
                         'DMMP\'s P centre is more susceptible\n'
                         'to nucleophilic attack than VX\'s\n'
                         f'by ${_ratio:.1f}' + r'\times$'),
        fontsize=FS_ANNOT, color=C_INK, va='center', ha='left',
        bbox=dict(fc='white', ec=C_INK, boxstyle='round,pad=0.3', alpha=0.92))
ex.set_xlabel(r'Mulliken charge on phosphorus,  $q_{\mathrm{P}}$  (e)')
ex.set_ylabel(r'condensed Fukui index at phosphorus,  $f^{+}(\mathrm{P})$'
              '\n(susceptibility to nucleophilic attack)')
dress(ex); panel_label(ex, 'e')

# ---- (f) 전역 omega 순위
E = df.sort_values('electrophilicity_index_eV').reset_index(drop=True)
y2 = np.arange(len(E))
fx.barh(y2, E.electrophilicity_index_eV, color=[CMAP[k] for k in E.compound_class],
        edgecolor='white', height=0.78, lw=1.0, zorder=2)
fx.set_yticks(y2); fx.set_yticklabels([s(n) for n in E.compound_name],
                                      fontsize=FS_TICK - 3)
wm = float(E.electrophilicity_index_eV.max())
fx.set_xlim(0, wm * 1.52); fx.set_ylim(-0.9, len(E) - 0.1)
fx.spines['left'].set_bounds(y2.min(), y2.max())
for y, r in zip(y2, E.itertuples()):
    fx.text(r.electrophilicity_index_eV + wm * 0.012, y,
            f'{r.electrophilicity_index_eV:.1f}', fontsize=FS_FOOT - 1,
            va='center', ha='left', color='#555555')
gm2 = E.compound_name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])
lo, hi = int(y2[gm2].min()), int(y2[gm2].max())
fx.annotate('', xy=(wm * 1.13, lo - 0.45), xytext=(wm * 1.13, hi + 0.45),
            arrowprops=dict(arrowstyle='-', color=C_REF, lw=3.0))
fx.text(wm * 1.16, (lo + hi) / 2,
        'omega is reported here and used in no\ndistance: it is an exact function of\n'
        'the two orbital energies. The three\nG-series agents still rank LOWEST of\n'
        'all 11 authentic agents; GA ranks 8th\nand VX 7th. 9 of 14 interferents\n'
        r'score higher than sarin.',
        fontsize=FS_ANNOT, color=C_REF, va='center', ha='left',
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.3', alpha=0.92))
fx.set_xlabel(r'Parr electrophilicity index,  $\omega = \mu^{2}/2\eta$  (eV)'
              '\n(derived quantity, not a coordinate of any space)')
fx.legend(handles=[Line2D([], [], marker='s', ls='', ms=12, mfc=c, mec='white',
                          label=l) for c, l in
                   [(C_REAL, 'authentic agent'), (C_SIM, 'incumbent simulant'),
                    (C_INTF, 'interferent / solvent')]],
          loc='lower right', frameon=False, bbox_to_anchor=(1.0, 0.02))
dress(fx); panel_label(fx, 'f')

fig.savefig('outputs/Fig_main.png', dpi=200, bbox_inches='tight', facecolor='white')
print('[fig] outputs/Fig_main.png')
