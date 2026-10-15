#!/usr/bin/env python3
"""QSAR 독성 모델 진단 그림 — 왜 이 모델로 모사체 안전성 축을 메울 수 없는가."""
import json
import numpy as np
import pandas as pd

import fig_annot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import qsar_tox as Q      # 재사용: desc(), FEAT 구성 로직은 아래서 재계산

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT = 28, 16, 15, 13, 13
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


TR, HO, FEAT = Q.TR, Q.HO, Q.FEAT
X, y = TR[FEAT].to_numpy(float), TR.y.to_numpy(float)
oof = np.zeros_like(y)
for s in Q.SEEDS:
    o = np.empty_like(y)
    for tr, te in KFold(5, shuffle=True, random_state=s).split(X):
        sc = StandardScaler().fit(X[tr])
        m = RandomForestRegressor(n_estimators=600, min_samples_leaf=2,
                                  max_features=0.4, random_state=s,
                                  n_jobs=-1).fit(sc.transform(X[tr]), y[tr])
        o[te] = m.predict(sc.transform(X[te]))
    oof += o / len(Q.SEEDS)

sc_all = StandardScaler().fit(X)
rf = RandomForestRegressor(n_estimators=600, min_samples_leaf=2, max_features=0.4,
                           random_state=11, n_jobs=-1).fit(sc_all.transform(X), y)
Hz = sc_all.transform(HO[FEAT].to_numpy(float))
hp, ht = rf.predict(Hz), HO.y.to_numpy()
hd = Q.knn_dist(Hz)
td = Q.knn_dist(sc_all.transform(X))
AD = Q.AD_T
SH = {'Agent GA (Tabun)': 'GA', 'Agent GD (Soman)': 'GD', 'Agent GB (Sarin)': 'GB',
      'Sulfur Mustard': 'HD', 'Lewisite 1': 'L', 'Agent VX': 'VX',
      'Nitrogen Mustard-1': 'HN-1', 'Nitrogen Mustard-2': 'HN-2',
      'Nitrogen Mustard-3': 'HN-3'}

fig, (ax, bx) = plt.subplots(1, 2, figsize=(22, 9.5))
fig.subplots_adjust(left=0.07, right=0.985, top=0.91, bottom=0.12, wspace=0.26)

# ---- (a) 예측 vs 실측
lim = (-1.6, 4.3)
ax.plot(lim, lim, color=C_INK, lw=1.6, ls=(0, (5, 3)), zorder=1)
ax.scatter(y, oof, s=95, c=C_INTF, edgecolors='white', lw=1.1, alpha=0.9,
           zorder=3, label=f'AEGL training set, out-of-fold (n={len(y)})')
for t, p, n in zip(ht, hp, HO.aegl_name):
    ax.annotate('', xy=(t, p), xytext=(t, t),
                arrowprops=dict(arrowstyle='->', color=C_REF, lw=2.0,
                                shrinkA=0, shrinkB=4), zorder=4)
ax.scatter(ht, hp, s=300, c=C_REAL, edgecolors='white', lw=2.2, zorder=5,
           label=f'chemical warfare agents, held out of training (n={len(ht)})')
from adjustText import adjust_text as _adj
_ta = [ax.text(t, p, SH.get(n, n), fontsize=FS_ANNOT + 1, fontweight='bold',
               color=C_INK, zorder=7) for t, p, n in zip(ht, hp, HO.aegl_name)]
_adj(_ta, x=np.asarray(ht), y=np.asarray(hp), ax=ax, only_move={'text': 'xy'},
     expand=(2.0, 2.4), force_text=(1.0, 1.4), force_static=(1.2, 1.6),
     max_move=110, iter_lim=900,
     arrowprops=dict(arrowstyle='-', color='#999999', lw=0.8))
ax.set_xlim(*lim); ax.set_ylim(lim[0], lim[1] + 0.9)
ax.set_xlabel(r'measured  log$_{10}$ AEGL-3 10 min (mg m$^{-3}$)')
ax.set_ylabel(r'predicted  log$_{10}$ AEGL-3 10 min (mg m$^{-3}$)')
ax.text(-1.45, 4.0, 'above the diagonal = model says SAFER than it is',
        fontsize=FS_ANNOT, color=C_REF, va='center', ha='left',
        bbox=dict(fc='white', ec='none', alpha=0.9))
ax.legend(loc='lower right', frameon=False, bbox_to_anchor=(1.0, 0.0))
dress(ax); panel_label(ax, 'a')

# ---- (b) 오차 vs 적용영역 거리
ax2 = bx
ax2.axhline(0, color=C_INK, lw=1.4, zorder=1)
ax2.axvline(AD, color=C_REF, lw=2.0, ls=(0, (5, 3)), zorder=1)
ax2.scatter(td, oof - y, s=85, c=C_INTF, edgecolors='white', lw=1.0, alpha=0.9,
            zorder=3, label='AEGL training set')
ax2.scatter(hd, hp - ht, s=300, c=C_REAL, edgecolors='white', lw=2.2, zorder=5,
            label='held-out CWAs')
_tb = [ax2.text(d, e, SH.get(n, n), fontsize=FS_ANNOT + 1, fontweight='bold',
                color=C_INK, zorder=7)
       for d, e, n in zip(hd, hp - ht, HO.aegl_name)]
D28 = pd.read_csv('outputs/toxicity_predicted.csv')
dm = D28[D28.compound_name == 'DMMP'].iloc[0]
ax2.axvline(dm.ad_dist, color=C_SIM, lw=2.4, ls=(0, (1, 2)), zorder=2)
ax2.set_xlim(0, 13.2); ax2.set_ylim(-3.2, 4.3)
ax2.text(AD + 0.15, -2.95, 'applicability-domain\nboundary (95th pct)',
         fontsize=FS_ANNOT, color=C_REF, va='bottom', ha='left')
ax2.text(dm.ad_dist - 0.2, 3.15, 'DMMP sits here', fontsize=FS_ANNOT + 1,
         color=C_SIM, fontweight='bold', va='center', ha='right',
         bbox=dict(fc='white', ec='none', alpha=0.9))
ax2.set_xlabel('applicability-domain distance (mean distance to 5 nearest\n'
               'training compounds, standardised descriptor space)')
ax2.set_ylabel(r'prediction error  (predicted $-$ measured, log$_{10}$)')
_err = np.asarray(hp) - np.asarray(ht)
_in = int((np.asarray(hd) <= AD).sum())
# 배수 표기는 본문과 같은 규칙을 쓴다. 작은 쪽은 소수 한 자리, 큰 쪽은
# 유효숫자 두 자리. 종전에는 %d 로 반올림해 2.51 이 '3' 으로 찍혀, 본문의
# '2.5 to 2600' 과 그림의 '3 to 2630' 이 어긋났다.
_lo = 10 ** float(_err.min())
_hi = float('%.2g' % (10 ** float(_err.max())))
_worst_in = float(_err[np.asarray(hd) <= AD].max())
ax2.text(0.35, -1.55,
         fig_annot.record('Fig_qsar_failure', 'panel_b',
             f'{_in} of the {len(hd)} agents pass the domain check,\n'
             f'and all {len(hd)} are wrong by {_err.min():.2f} to {_err.max():.2f} log\n'
             f'(a factor of {_lo:.1f} to {_hi:,.0f}), always in the unsafe direction.\n'
             f'The worst inside the domain is wrong by {_worst_in:.2f} log.'),
         fontsize=FS_ANNOT, color=C_REF, va='top', ha='left',
         bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.35', alpha=0.95))
_adj(_tb, x=np.asarray(hd), y=_err, ax=ax2, only_move={'text': 'xy'},
     expand=(2.0, 2.4), force_text=(1.0, 1.4), force_static=(1.2, 1.6),
     max_move=110, iter_lim=900,
     arrowprops=dict(arrowstyle='-', color='#999999', lw=0.8))
ax2.legend(loc='upper right', frameon=False, bbox_to_anchor=(1.0, 0.92))
dress(ax2); panel_label(ax2, 'b')

fig.savefig('outputs/Fig_qsar_failure.png', dpi=200, bbox_inches='tight',
            facecolor='white')
print('[fig] outputs/Fig_qsar_failure.png')
