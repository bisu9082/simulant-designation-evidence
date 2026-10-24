#!/usr/bin/env python3
"""simulant / STEP 4 — GFN2-xTB 전자 서술자 Figure (4패널)"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

# ---- Ku 승인 스케일 / 색맹검증 통과 3색 --------------------------------------
FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 28, 16, 15, 13, 13, 11
C_REAL, C_SIM, C_INTF = '#C94F4A', '#3A6EA8', '#E8943A'
C_REF = '#B84040'
CMAP = {'real_agent': C_REAL, 'simulant': C_SIM, 'interferent': C_INTF}
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 15,
                     'axes.labelsize': FS_AXIS, 'xtick.labelsize': FS_TICK,
                     'ytick.labelsize': FS_TICK, 'legend.fontsize': FS_LEG,
                     'axes.facecolor': 'white', 'figure.facecolor': 'white'})


def panel_label(ax, letter, fs=FS_LABEL, dy=6):
    ax.annotate(f'({letter})', xy=(0.0, 1.0), xycoords='axes fraction',
                xytext=(0, dy), textcoords='offset points', fontsize=fs,
                fontweight='bold', va='bottom', ha='left', annotation_clip=False)


def dress(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(color='#DDDDDD', ls='--', lw=0.6)
    ax.set_axisbelow(True)


df = pd.read_csv('outputs/_electronic_with_gap.csv')
Y = np.load('outputs/_pca_scores.npy')
T = pd.read_csv('outputs/electronic_distance_table.csv')
S = json.load(open('outputs/electronic_space_summary.json'))
evr = S['pca']['explained_variance_ratio']

SHORT = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
         'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD',
         'lewisite 1 (L)': 'L', 'hydrogen cyanide (AC)': 'AC',
         'cyanogen chloride (CK)': 'CK', 'phosgene (CG)': 'CG',
         'chloropicrin (PS)': 'PS', 'dichloromethane': 'DCM',
         'p-dichlorobenzene': 'p-DCB', 'ethyl acetate': 'EtOAc',
         'acetonitrile': 'ACN', 'n-hexane': 'hexane'}
short = lambda s: SHORT.get(s, s)

fig, axes = plt.subplots(2, 2, figsize=(22, 15))
fig.subplots_adjust(left=0.07, right=0.98, top=0.93, bottom=0.07,
                    hspace=0.30, wspace=0.22)
(ax, bx), (cx, dx) = axes

# ================================================= (a) 전자공간 PCA
for k in ['interferent', 'simulant', 'real_agent']:
    m = (df.klass == k).to_numpy()
    ax.scatter(Y[m, 0], Y[m, 1], s=190, c=CMAP[k], edgecolors='white',
               linewidths=1.6, zorder=3,
               label={'real_agent': 'authentic agent', 'simulant': 'incumbent simulant',
                      'interferent': 'interferent / solvent'}[k])
txt = [ax.text(Y[i, 0], Y[i, 1], short(df.name[i]), fontsize=FS_ANNOT - 1,
               color='#333333', zorder=4) for i in range(len(df))]
ax.set_xlim(Y[:, 0].min() - 0.9, Y[:, 0].max() + 0.9)
ax.set_ylim(Y[:, 1].min() - 0.5, Y[:, 1].max() + 2.1)      # R-15 헤드룸
adjust_text(txt, ax=ax, only_move={'text': 'xy'}, expand=(1.6, 1.9),
            force_text=(0.6, 0.9), force_static=(0.5, 0.8), max_move=60,
            iter_lim=400,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7))
ax.set_xlabel(f'PC1 ({evr[0]*100:.1f}%):  global electrophilicity  '
              r'($\omega$ $+$, $E_{\mathrm{LUMO}}$ $-$)')
ax.set_ylabel(f'PC2 ({evr[1]*100:.1f}%):  polarity  '
              r'($\mu$ $+$, $E_{\mathrm{HOMO}}$ $-$)')
ax.legend(loc='upper left', frameon=False, ncol=1)
dress(ax); panel_label(ax, 'a')

# ================================================= (b) 최근접 이웃 덤벨
T = T.sort_values('d_simulant', ascending=True).reset_index(drop=True)
yy = np.arange(len(T))[::-1]
for y, r in zip(yy, T.itertuples()):
    bx.plot([r.d_simulant, r.d_interferent], [y, y], color='#BBBBBB',
            lw=2.2, zorder=1, solid_capstyle='round')
bx.scatter(T.d_simulant, yy, s=200, c=C_SIM, edgecolors='white', lw=1.6,
           zorder=3, label='nearest incumbent simulant')
bx.scatter(T.d_interferent, yy, s=200, c=C_INTF, edgecolors='white', lw=1.6,
           zorder=3, label='nearest interferent')
bx.set_yticks(yy); bx.set_yticklabels([short(n) for n in T.real_agent])
xmax = float(max(T.d_simulant.max(), T.d_interferent.max()))
bx.set_xlim(0, xmax * 1.22)
bx.set_ylim(-1.45, len(T) - 0.4 + 1.30)                     # R-15 양측 헤드룸
bx.spines['left'].set_bounds(yy.min(), yy.max())
bx.set_xlabel('standardised Euclidean distance in electronic space')
# ✕ 는 각 행의 가장 오른쪽 점 바로 옆에 — 데이터와 분리되어 뜨지 않도록
for y, r in zip(yy, T.itertuples()):
    if r.interferent_closer:
        bx.text(max(r.d_simulant, r.d_interferent) + xmax * 0.055, y, '✕',
                fontsize=FS_ANNOT + 2, color=C_REF, va='center', ha='left',
                fontweight='bold')
bx.legend(loc='upper right', frameon=False,
          bbox_to_anchor=(1.0, 1.0))
bx.text(xmax * 0.02, -1.08,
        '✕  an interferent lies closer than any incumbent simulant '
        f'({int(T.interferent_closer.sum())} of {len(T)} authentic agents)',
        fontsize=FS_ANNOT, color=C_REF, va='center', ha='left')
dress(bx); panel_label(bx, 'b')

# ================================================= (c) 인 중심 국소 반응성
P = df[df.q_P.notna()].reset_index(drop=True)
gm, gs = S['p_local']['G_series_mean'], S['p_local']['G_series_sd']
cx.axhspan(gm - gs, gm + gs, color=C_REAL, alpha=0.13, zorder=0)
cx.axhline(gm, color=C_REF, ls=(0, (5, 3)), lw=1.8, zorder=1)
for _, r in P.iterrows():
    cx.scatter(r.q_P, r.fukui_f_plus_P, s=300, c=CMAP[r.klass],
               edgecolors='white', lw=2.0, zorder=3)
# 6점뿐이므로 포인트 단위 수동 오프셋으로 확정 배치 (GD/GF 는 y 차이 0.005)
OFF = {'sarin (GB)': (14, 7), 'soman (GD)': (16, -3), 'cyclosarin (GF)': (14, -19),
       'DMMP': (16, 0), 'tabun (GA)': (16, 2), 'VX': (16, 0)}
for _, r in P.iterrows():
    cx.annotate(short(r['name']), xy=(r.q_P, r.fukui_f_plus_P),
                xytext=OFF[r['name']], textcoords='offset points',
                fontsize=FS_ANNOT + 1, color='#333333', va='center', ha='left',
                fontweight='bold' if r.klass == 'simulant' else 'normal')
cx.set_xlim(0.28, 0.62)
cx.set_ylim(0.10, 0.75 + 0.14)                              # R-15 헤드룸
cx.text(0.292, gm + gs + 0.018,
        r'G-series mean $\pm$ SD  ($0.646 \pm 0.019$)',
        fontsize=FS_ANNOT, color=C_REF, va='bottom', ha='left',
        bbox=dict(fc='white', ec='none', alpha=0.9))
cx.annotate('', xy=(0.300, 0.2043), xytext=(0.300, 0.5585),
            arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.6))
cx.text(0.392, 0.300, 'DMMP overstates\nV-agent P-electrophilicity\n'
        r'by $2.7\times$', fontsize=FS_ANNOT, color='#333333', va='center',
        ha='left', bbox=dict(fc='white', ec='#333333',
                             boxstyle='round,pad=0.3', alpha=0.92))
cx.set_xlabel(r'Mulliken charge on phosphorus,  $q_{\mathrm{P}}$  (e)')
cx.set_ylabel(r'local electrophilicity at phosphorus,  $f^{+}(\mathrm{P})$')
dress(cx); panel_label(cx, 'c')

# ================================================= (d) 전역 omega 의 실패
E = df.sort_values('electrophilicity_index_eV').reset_index(drop=True)
yy2 = np.arange(len(E))
cols = [CMAP[k] for k in E.klass]
dx.barh(yy2, E.electrophilicity_index_eV, color=cols, edgecolor='white',
        height=0.78, lw=1.0, zorder=2)
dx.set_yticks(yy2)
dx.set_yticklabels([short(n) for n in E.name], fontsize=FS_TICK - 3)
wmax = float(E.electrophilicity_index_eV.max())
dx.set_xlim(0, wmax * 1.52)
dx.set_ylim(-0.9, len(E) - 0.1)
dx.spines['left'].set_bounds(yy2.min(), yy2.max())
for y, r in zip(yy2, E.itertuples()):
    dx.text(r.electrophilicity_index_eV + wmax * 0.012, y,
            f'{r.electrophilicity_index_eV:.1f}', fontsize=FS_FOOT - 1,
            va='center', ha='left', color='#555555')
gmask = E.name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])
lo, hi = int(yy2[gmask].min()), int(yy2[gmask].max())
dx.annotate('', xy=(wmax * 1.13, lo - 0.45), xytext=(wmax * 1.13, hi + 0.45),
            arrowprops=dict(arrowstyle='-', color=C_REF, lw=3.0))
dx.text(wmax * 1.16, (lo + hi) / 2,
        'the three nerve agents rank\nLOWEST of all 11 authentic agents;\n'
        '10 of 15 interferents score\n'
        r'higher than sarin on global $\omega$',
        fontsize=FS_ANNOT, color=C_REF, va='center',
        ha='left', bbox=dict(fc='white', ec=C_REF,
                             boxstyle='round,pad=0.3', alpha=0.92))
dx.set_xlabel(r'global electrophilicity index,  $\omega = \mu^{2}/2\eta$  (eV)')
dx.legend(handles=[Line2D([], [], marker='s', ls='', ms=12, mfc=c, mec='white',
                          label=l) for c, l in
                   [(C_REAL, 'authentic agent'), (C_SIM, 'incumbent simulant'),
                    (C_INTF, 'interferent / solvent')]],
          loc='lower right', frameon=False, bbox_to_anchor=(1.0, 0.02))
dress(dx); panel_label(dx, 'd')

fig.savefig('outputs/Fig_electronic_space.png', dpi=200,
            bbox_inches='tight', facecolor='white')
print('[fig] outputs/Fig_electronic_space.png')
