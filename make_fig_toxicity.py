#!/usr/bin/env python3
"""simulant — 안전성 축 도판 (본문 Figure 3).

(a) EPA AEGL-3 10분 사다리  (b) 결측 구조  (c) QSAR 예측 대 실측
(d) 오차 대 적용영역 거리. 2026-10-20 에 make_fig_qsar.py 를 흡수했다.
본문 도판을 세 장, 각 네 패널로 맞추기 위한 통합이다.
"""
import json
import textwrap
import numpy as np
import pandas as pd

import fig_annot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text as _adj
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

import qsar_tox as Q

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 19, 15, 14, 12, 12, 10
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

# figsize 는 \textwidth(384 pt = 5.33 in) 기준으로 정한다. 24~27 in 로
# 그리면 축소율이 0.2 가 되어 15 pt 틱이 실효 3 pt 로 인쇄된다 (2026-10-20
# 외부 검토 지적, 전 도판 해당). 10 in 안쪽으로 그려 실효 7 pt 이상을 둔다.
# 2x2 로 둔다. 3행(전폭 2개)으로 그리면 높이가 폭의 1.38배가 되고, 본문
# \textwidth 5.33 in 에 맞추면 그림만 7.3 in 이 되어 \textheight 562 pt
# (7.81 in) 안에 캡션 10행이 들어가지 못한다 (쪽번호를 덮었다, 2026-10-20).
# 2x2 는 높이비 1.1 로 떨어져 그림 5.9 in + 캡션 1.7 in 으로 들어간다.
# 대신 (c)(d) 가 반폭이 되므로 점별 라벨은 본 연구 패널의 작용제로 한정한다.
fig = plt.figure(figsize=(10.0, 11.0))
_gs = fig.add_gridspec(2, 2, height_ratios=[1.06, 1.0],
                       hspace=0.30, wspace=0.42,
                       left=0.135, right=0.986, top=0.952, bottom=0.062)
ax = fig.add_subplot(_gs[0, 0])
bx = fig.add_subplot(_gs[0, 1])
cx = fig.add_subplot(_gs[1, 0])
dx = fig.add_subplot(_gs[1, 1])

# ---------------- (a) AEGL-3 사다리 (로그)
y = np.arange(len(H))
ax.barh(y, H.AEGL3_10min_mgm3, color=[CMAP[k] for k in H.compound_class],
        edgecolor='white', lw=0.7, height=0.74, zorder=2)
SHORT_A = {k: v.split('  ')[0] for k, v in SHORT.items() if '  ' in v}
ax.set_yticks(y)
ax.set_yticklabels([SHORT_A.get(n, n) for n in H.compound_name])
ax.set_xscale('log'); ax.set_xlim(8e-3, 6e5)
ax.set_ylim(-1.7, len(H) - 0.35 + 1.25)
ax.spines['left'].set_bounds(0, len(H) - 1)
for yy, r in zip(y, H.itertuples()):
    ax.text(r.AEGL3_10min_mgm3 * 1.35, yy, f'{r.AEGL3_10min_mgm3:.3g}',
            fontsize=FS_ANNOT, va='center', ha='left', color='#555555')
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
# 패널 안에는 둘 자리가 없다. 막대 19개와 범례가 이미 폭을 다 쓴다.
# 축 위로 올린다 (fig_annot 등록은 유지).
# wrap=True 는 figure 폭 기준이라 옆 패널까지 가로지른다 (2026-10-20).
# 패널 폭(약 5.3 in)에 맞춰 직접 줄을 끊는다.
ax.set_title(textwrap.fill(_txt_a.replace('\n', ' '), 68),
             fontsize=FS_ANNOT - 2, color=C_REF, loc='left', pad=34)
ax.set_xlabel('EPA AEGL-3, 10 min (mg m$^{-3}$, log scale)\n'
              'lower = lethal at lower airborne concentration')
ax.legend(handles=[Line2D([], [], marker='s', ls='', ms=9.1, mfc=c, mec='white', label=l)
                   for c, l in [(C_REAL, 'authentic agent'),
                                (C_SIM, 'designated simulant'),
                                (C_INTF, 'interferent / solvent')]],
          loc='upper right', frameon=False, bbox_to_anchor=(1.0, 1.0))
dress(ax); panel_label(ax, 'a')

# ---------------- (b) 결측 구조
y2 = np.arange(len(M))[::-1]
bx.scatter(np.zeros(len(M)), y2, marker='x', s=43,
           c=[CMAP[k] for k in M.compound_class], lw=2.0, zorder=3)
SHORT_B = {'cyanogen chloride (CK)': 'CK', 'p-dichlorobenzene': 'p-DCB',
           'ethyl acetate': 'ethyl acetate'}
bx.set_yticks(y2[::-1])
bx.set_yticklabels([SHORT_B.get(n, n) for n in M.compound_name][::-1])
REASON = {'[수집 불가: EPA AEGL 미등재]': 'not on the AEGL list at all',
          '[수집 불가: EPA AEGL 보류목록, 값 미개발]':
              'on the holding list; values\nnot being developed',
          '[수집 불가: AEGL-3 = NR]': 'listed, but AEGL-3 = "not recommended,\n'
                                      'insufficient data"'}
for yy, r in zip(y2, M.itertuples()):
    bx.text(0.09, yy, REASON.get(r.note, r.note), fontsize=FS_ANNOT - 1,
            va='center', ha='left', color='#555555')
bx.set_xlim(-0.06, 1.0); bx.set_ylim(-4.8, len(M) - 0.35 + 1.55)
bx.spines['left'].set_bounds(0, len(M) - 1)
bx.set_xticks([])
bx.spines['bottom'].set_visible(False)
bx.text(-0.045, -2.75,
        textwrap.fill('DMMP, the simulant carrying every valid match in this '
                      'study, has no AEGL. EPA PPRTV (EPA/690/R-06/014F, 2006) '
                      'declined to derive an inhalation RfC ("data are '
                      'insufficient") and derived no inhalation unit risk.', 44),
        fontsize=FS_ANNOT - 2, color=C_REF, va='center', ha='left',
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.4', alpha=0.93))
bx.text(-0.045, len(M) - 0.35 + 0.72,
        fig_annot.record('Fig_toxicity', 'panel_b_header',
            f'no regulator-grade acute inhalation\nvalue exists for any of these {len(M)}'),
        fontsize=FS_ANNOT + 1, color=C_INK, va='center', ha='left')
dress(bx); bx.grid(False); panel_label(bx, 'b')

# ===== (c)(d) QSAR 진단. make_fig_qsar.py 흡수 (2026-10-20) =====
TR, HO, FEAT = Q.TR, Q.HO, Q.FEAT
X, y = TR[FEAT].to_numpy(float), TR.y.to_numpy(float)
oof = np.zeros_like(y)
for _sd in Q.SEEDS:
    o = np.empty_like(y)
    for tr, te in KFold(5, shuffle=True, random_state=_sd).split(X):
        sc = StandardScaler().fit(X[tr])
        m = RandomForestRegressor(n_estimators=600, min_samples_leaf=2,
                                  max_features=0.4, random_state=_sd,
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
      'GF Agent': 'GF', 'Sulfur Mustard': 'HD', 'Lewisite 1': 'L', 'Agent VX': 'VX',
      'Hydrogen cyanide': 'AC', 'Phosgene': 'CG', 'Chloropicrin': 'PS',
      'Nitrogen Mustard-1': 'HN-1', 'Nitrogen Mustard-2': 'HN-2',
      'Nitrogen Mustard-3': 'HN-3'}
# 라벨은 SH 에 있는 것만 적는다. 나머지 보류 작용제는 붉은 점으로만 둔다.
LBL = [i for i, n in enumerate(HO.aegl_name) if n in SH]
# 캡션이 '본 연구 패널의 작용제 N종과 질소 머스터드 3종' 이라고 쓰므로 그 수를
# 산출물로 남긴다. 손으로 센 수가 캡션에 남아 썩는 경로를 막는다 (2026-10-20).
PANEL11 = {'GB', 'GD', 'GA', 'GF', 'HD', 'L', 'AC', 'CK', 'CG', 'PS', 'VX'}
_codes = [SH[HO.aegl_name.iloc[i]] for i in LBL]
json.dump({'n_heldout': int(len(HO)), 'n_labelled': len(LBL),
           'n_unlabelled': int(len(HO)) - len(LBL),
           'n_from_panel': sum(c in PANEL11 for c in _codes),
           'n_other_labelled': sum(c not in PANEL11 for c in _codes),
           'labelled_codes': _codes},
          open('outputs/fig_label_counts.json', 'w'), indent=1)


# ---- (c) 예측 vs 실측
lim = (-1.6, 4.3)
cx.plot(lim, lim, color=C_INK, lw=1.2, ls=(0, (5, 3)), zorder=1)
cx.scatter(y, oof, s=25, c=C_INTF, edgecolors='white', lw=0.8, alpha=0.9,
           zorder=3, label=f'training set, out-of-fold (n={len(y)})')
for t, p, n in zip(ht, hp, HO.aegl_name):
    cx.annotate('', xy=(t, p), xytext=(t, t),
                arrowprops=dict(arrowstyle='->', color=C_REF, lw=1.5,
                                shrinkA=0, shrinkB=4), zorder=4)
cx.scatter(ht, hp, s=79, c=C_REAL, edgecolors='white', lw=1.6, zorder=5,
           label=f'agents held out of training (n={len(ht)})')
_ta = [cx.text(ht[i], hp[i], SH[HO.aegl_name.iloc[i]], fontsize=FS_ANNOT + 1,
               fontweight='bold', color=C_INK, zorder=7) for i in LBL]
_adj(_ta, x=np.asarray(ht), y=np.asarray(hp), ax=cx, only_move={'text': 'xy'},
     expand=(2.0, 2.4), force_text=(1.0, 1.4), force_static=(1.2, 1.6),
     max_move=110, iter_lim=900,
     arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
cx.set_xlim(lim[0] - 0.45, lim[1] + 0.35); cx.set_ylim(lim[0] - 0.55, lim[1] + 1.1)
cx.set_xlabel('measured  ' r'log$_{10}$ AEGL-3 10 min' '\n' r'(mg m$^{-3}$)')
cx.set_ylabel('predicted\n' r'(log$_{10}$ mg m$^{-3}$)')
dress(cx); panel_label(cx, 'c')

# ---- (d) 오차 vs 적용영역 거리
dx.axhline(0, color=C_INK, lw=1.0, zorder=1)
dx.axvline(AD, color=C_REF, lw=1.5, ls=(0, (5, 3)), zorder=1)
dx.scatter(td, oof - y, s=24, c=C_INTF, edgecolors='white', lw=0.7, alpha=0.9,
            zorder=3, label='AEGL training set')
dx.scatter(hd, hp - ht, s=79, c=C_REAL, edgecolors='white', lw=1.6, zorder=5,
            label='held-out CWAs')
_tb = [dx.text(hd[i], (hp - ht)[i], SH[HO.aegl_name.iloc[i]],
               fontsize=FS_ANNOT + 1, fontweight='bold', color=C_INK, zorder=7)
       for i in LBL]
D28 = pd.read_csv('outputs/toxicity_predicted.csv')
dm = D28[D28.compound_name == 'DMMP'].iloc[0]
dx.axvline(dm.ad_dist, color=C_SIM, lw=1.8, ls=(0, (1, 2)), zorder=2)
dx.set_xlim(-0.35, 13.5); dx.set_ylim(-3.5, 4.7)
dx.text(AD + 0.15, -2.95, 'applicability-domain\nboundary (95th pct)',
         fontsize=FS_ANNOT, color=C_REF, va='bottom', ha='left')
# 가로로 쓰면 왼쪽 축 라벨 위로 넘어간다 (반폭 패널, 2026-10-20). 점선을
# 따라 세로로 세운다.
dx.text(dm.ad_dist - 0.22, -3.35, 'DMMP sits here', fontsize=FS_ANNOT,
        color=C_SIM, fontweight='bold', rotation=90, va='bottom', ha='right')
dx.set_xlabel('applicability-domain distance\n'
              '(mean distance to 5 nearest training compounds)')
dx.set_ylabel('prediction error\n' r'(predicted $-$ measured, log$_{10}$)')
_err = np.asarray(hp) - np.asarray(ht)
_in = int((np.asarray(hd) <= AD).sum())
# 배수 표기는 본문과 같은 규칙을 쓴다. 작은 쪽은 소수 한 자리, 큰 쪽은
# 유효숫자 두 자리. 종전에는 %d 로 반올림해 2.51 이 '3' 으로 찍혀, 본문의
# '2.5 to 2600' 과 그림의 '3 to 2630' 이 어긋났다.
_lo = 10 ** float(_err.min())
_hi = float('%.2g' % (10 ** float(_err.max())))
_worst_in = float(_err[np.asarray(hd) <= AD].max())
# 종전에는 이 문장을 패널 안에 그려서 화합물 라벨을 덮었다. 축 밖으로
# 내보내되 fig_annot 등록은 유지해야 검증기가 계속 수치를 대조한다.
_msg = fig_annot.record('Fig_toxicity', 'panel_d',
    f'{_in} of the {len(hd)} agents pass the domain check, '
    f'and all {len(hd)} are wrong by {_err.min():.2f} to {_err.max():.2f} log '
    f'(a factor of {_lo:.1f} to {_hi:,.0f}), always in the unsafe direction. '
    f'The worst inside the domain is wrong by {_worst_in:.2f} log.')
# 축 밖 여백은 bbox_inches='tight' 가 잘라내므로 패널 제목으로 올린다.
dx.set_title(textwrap.fill(_msg, 66), fontsize=FS_ANNOT - 2, color=C_REF,
             loc='left', pad=34)
_adj(_tb, x=np.asarray(hd), y=_err, ax=dx, only_move={'text': 'xy'},
     expand=(2.0, 2.4), force_text=(1.0, 1.4), force_static=(1.2, 1.6),
     max_move=110, iter_lim=900,
     arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
dress(dx); panel_label(dx, 'd')


fig.savefig('outputs/Fig_toxicity.png', dpi=200, bbox_inches='tight', facecolor='white')
print('[fig] outputs/Fig_toxicity.png')
print(f'  (a) {len(H)}종  (b) {len(M)}종  (c)(d) 보류 {len(hd)}종')
