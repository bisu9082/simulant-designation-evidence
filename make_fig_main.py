#!/usr/bin/env python3
"""simulant / STEP 4 — 본문 도판 Figure 1, Figure 2 와 보충 Fig (인 국소반응성).

2026-10-20 에 make_fig_panel.py 를 흡수했다. 본문 도판을 세 장, 각 네 패널로
맞추면서 한 스크립트에서 세 figure 를 모두 그린다. 정본 라벨(모사체 3종) 기준.
"""
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

FS_LABEL, FS_AXIS, FS_TICK, FS_LEG, FS_ANNOT, FS_FOOT = 19, 15, 14, 12, 12, 10
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


# figsize 는 \textwidth(384 pt = 5.33 in) 기준으로 정한다. 24~27 in 로
# 그리면 축소율이 0.2 가 되어 15 pt 틱이 실효 3 pt 로 인쇄된다 (2026-10-20
# 외부 검토 지적, 전 도판 해당). 10 in 안쪽으로 그려 실효 7 pt 이상을 둔다.
# 2026-10-20 재구성. \textwidth 가 5.33 in 이므로 가로 3패널은 패널당 1.7 in
# 이 되어 라벨이 들어가지 않는다. 동시에 패널 2개짜리 도판은 이 저널의 관행과
# 어긋난다 (J. Hazard. Mater. 501 (2026) 140703 은 8쪽에 도판 5장, 패널 1/5/4/2/9).
# 2x2 를 기본 격자로 삼아 본문 도판 3장을 각 4패널로 구성한다.
#   Figure 1  비교 구성이 결론을 정한다      (a) 작용제별 거리 (b) 부트스트랩
#                                             (c) 공간별 거리 (d) 순열 귀무분포
#   Figure 2  공간과 명목값의 값어치         (a)(b) 주성분 (c) 부분공간 분포
#                                             (d) 보정 체계
#   Fig S     인 국소반응성과 omega 순위     (a)(b)
# Figure 1 은 (a) 작용제별 거리와 (b) 부트스트랩이 가로폭을 요구하므로 (둘 다
# 마커 오른쪽에 글자가 붙는다) 그 둘을 전폭 행에 두고, (c)(d) 를 아래에서 나눈다.
fig = plt.figure(figsize=(10.0, 12.0))
_g1 = fig.add_gridspec(3, 2, height_ratios=[1.22, 0.66, 1.05],
                       hspace=0.46, wspace=0.40,
                       left=0.098, right=0.980, top=0.962, bottom=0.062)
gx = fig.add_subplot(_g1[0, :])
hx = fig.add_subplot(_g1[1, :])
cx = fig.add_subplot(_g1[2, 0])
dx = fig.add_subplot(_g1[2, 1])

fig2, _f2 = plt.subplots(2, 2, figsize=(10.0, 8.8))
fig2.subplots_adjust(left=0.090, right=0.982, top=0.950, bottom=0.078,
                     hspace=0.44, wspace=0.32)
ax, bx = _f2[0]
sx, ix = _f2[1]

# (b) 는 막대 19개에 화합물 이름과 값이 모두 붙으므로 세로 여유가 필요하다.
# 4.6 in 으로 그렸더니 행 간격이 13 pt 가 되어 라벨이 서로 겹쳤다 (2026-10-20).
figS = plt.figure(figsize=(10.0, 6.9))
_gsS = figS.add_gridspec(1, 2, width_ratios=[1.0, 1.18], wspace=0.30,
                         left=0.088, right=0.984, top=0.940, bottom=0.145)
ex = figS.add_subplot(_gsS[0, 0])
fx = figS.add_subplot(_gsS[0, 1])


def pca_panel(axis, cond, letter, xlab):
    Y, evr, keep, comp = scores(cond)
    for k in ['interferent', 'simulant', 'real_agent']:
        m = (df.compound_class == k).to_numpy()
        axis.scatter(Y[m, 0], Y[m, 1], s=38, c=CMAP[k], edgecolors='white',
                     lw=1.1, zorder=3,
                     label={'real_agent': 'authentic agent',
                            'simulant': 'incumbent simulant',
                            'interferent': 'interferent / solvent'}[k])
    # 2026-10-20. 반폭 패널에서 28개 라벨은 10 pt 로 그려도 인쇄 실효가
    # 5.3 pt 라 읽히지 않는다. 작용제와 모사체 14개만 적고 간섭물은 색으로만
    # 구분한다. 색 대응은 캡션에 적는다.
    lab_i = [i for i in range(len(df))
             if df.compound_class[i] in ('real_agent', 'simulant')]
    t = [axis.text(Y[i, 0], Y[i, 1], s(df.compound_name[i]), fontsize=FS_ANNOT,
                   color=C_INK, zorder=4) for i in lab_i]
    axis.set_xlim(Y[:, 0].min() - 1.0, Y[:, 0].max() + 1.0)
    axis.set_ylim(Y[:, 1].min() - 1.1, Y[:, 1].max() + 1.3)
    adjust_text(t, x=Y[:, 0], y=Y[:, 1], ax=axis, only_move={'text': 'xy'},
                expand=(1.9, 2.3), force_text=(1.1, 1.5), force_static=(1.3, 1.7),
                force_pull=(0.004, 0.004), max_move=140, iter_lim=2000,
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.6))
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
        cx.scatter(r.d_sim_norm, yy[r.real_agent], marker=MARK[cond], s=32,
                   facecolors=(C_SIM if ok else 'white'), edgecolors=C_INK,
                   lw=1.0, zorder=3)
cx.axvline(2.0, color=C_REF, ls=(0, (5, 3)), lw=1.3, zorder=1)
cx.set_yticks(range(len(order)))
cx.set_yticklabels([s(n) for n in order][::-1])
xm = float(R.d_sim_norm.max())
cx.set_xlim(0, xm * 1.16); cx.set_ylim(-4.6, len(order) - 0.45 + 1.15)
cx.spines['left'].set_bounds(0, len(order) - 1)
cx.set_xlabel('distance to nearest incumbent simulant\n'
              '(typical nearest-neighbour units)')
cx.text(2.08, len(order) - 0.6, 'match\nthreshold', fontsize=FS_ANNOT, color=C_REF,
        va='top', ha='left')
cx.legend(handles=[Line2D([], [], marker=MARK[c], ls='', ms=8.2, mfc='white',
                          mec=C_INK, mew=1.0, label=SPACE_LBL[c])
                   for c in MARK] +
                  [Line2D([], [], marker='o', ls='', ms=8.2, mfc=C_SIM, mec=C_INK,
                          mew=1.0, label='valid match (also beats every interferent)')],
          loc='lower right', frameon=False, fontsize=FS_LEG - 3, bbox_to_anchor=(1.03, -0.015),
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
    dx.scatter([null.mean()], [y], marker='|', s=196, c=C_INK, lw=1.6, zorder=3)
    k = SAN[c]
    dx.scatter([k['observed_diff']], [y], s=65, c=C_REAL, edgecolors='white',
               lw=1.4, zorder=5)
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
dx.set_xlabel(r'mean ( $d_\mathrm{sim} - d_\mathrm{intf}$ )'
              '\ngrey: all 680 designations of 3 of 17')
dx.legend(handles=[Line2D([], [], marker='o', ls='', ms=6.7, mfc='#BBBBBB',
                          mec='none', label='null: 3 compounds picked at random'),
                   Line2D([], [], marker='|', ls='', ms=11.9, mec=C_INK, mew=1.6,
                          label='null mean'),
                   Line2D([], [], marker='o', ls='', ms=9.6, mfc=C_REAL,
                          mec='white', mew=1.3, label='actual incumbent simulants')],
          loc='lower right', frameon=False, fontsize=FS_LEG - 3, bbox_to_anchor=(1.03, -0.015),
          handletextpad=0.6, labelspacing=0.5)
dress(dx); panel_label(dx, 'd')

# ---- (e) 인 중심 국소 반응성
P = df[df.q_P.notna()].reset_index(drop=True)
g3 = P[P.compound_name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])]
gm, gs = g3.fukui_f_plus_P.mean(), g3.fukui_f_plus_P.std(ddof=1)
ex.axhspan(gm - gs, gm + gs, color=C_REAL, alpha=0.13, zorder=0)
ex.axhline(gm, color=C_REF, ls=(0, (5, 3)), lw=1.2, zorder=1)
for r in P.itertuples():
    ex.scatter(r.q_P, r.fukui_f_plus_P, s=63, c=CMAP[r.compound_class],
               edgecolors='white', lw=1.3, zorder=3)
OFF = {'sarin (GB)': (14, 9), 'soman (GD)': (16, -3), 'cyclosarin (GF)': (14, -19),
       'DMMP': (16, 0), 'tabun (GA)': (16, 2), 'VX': (16, 0)}
for r in P.itertuples():
    ex.annotate(s(r.compound_name), xy=(r.q_P, r.fukui_f_plus_P),
                xytext=OFF[r.compound_name], textcoords='offset points',
                fontsize=FS_ANNOT + 1, color=C_INK, va='center', ha='left',
                fontweight='bold' if r.compound_class == 'simulant' else 'normal')
ex.set_xlim(0.28, 0.62); ex.set_ylim(0.10, 0.89)
# 한 줄로 쓰면 오른쪽 끝이 GB 라벨에 닿아 흰 bbox 가 그것을 덮었다 (2026-10-20).
ex.text(0.292, gm + gs + 0.018,
        f'GB/GD/GF mean $\\pm$ SD\n(${gm:.3f} \\pm {gs:.3f}$)',
        fontsize=FS_ANNOT, color=C_REF, va='bottom', ha='left', zorder=2,
        bbox=dict(fc='white', ec='none', alpha=0.9))
ex.annotate('', xy=(0.300, 0.2043), xytext=(0.300, 0.5585),
            arrowprops=dict(arrowstyle='<->', color=C_INK, lw=1.1))
# 수치를 데이터에서 계산한다 (2026-10-12). 종전에는 2.7 이 문자열 리터럴이었다.
_f = dict(zip(P6.compound_name, P6.fukui_f_plus_P)) if 'P6' in dir() else \
     dict(zip(df.compound_name, df.fukui_f_plus_P))
_ratio = _f['DMMP'] / _f['VX']
ex.text(0.352, 0.300,
        fig_annot.record('Fig_plocal', 'panel_a',
                         'DMMP\'s P centre is more susceptible\n'
                         'to nucleophilic attack than VX\'s\n'
                         f'by ${_ratio:.1f}' + r'\times$'),
        fontsize=FS_ANNOT, color=C_INK, va='center', ha='left',
        bbox=dict(fc='white', ec=C_INK, boxstyle='round,pad=0.3', alpha=0.92))
ex.set_xlabel(r'Mulliken charge on phosphorus,  $q_{\mathrm{P}}$  (e)')
ex.set_ylabel(r'condensed Fukui index, $f^{+}(\mathrm{P})$'
              '\n(nucleophilic attack at P)')
dress(ex); panel_label(ex, 'a')

# ---- (f) 전역 omega 순위
E = df.sort_values('electrophilicity_index_eV').reset_index(drop=True)
y2 = np.arange(len(E))
fx.barh(y2, E.electrophilicity_index_eV, color=[CMAP[k] for k in E.compound_class],
        edgecolor='white', height=0.78, lw=0.7, zorder=2)
fx.set_yticks(y2); fx.set_yticklabels([s(n) for n in E.compound_name],
                                      fontsize=FS_TICK - 3)
wm = float(E.electrophilicity_index_eV.max())
fx.set_xlim(0, wm * 1.66); fx.set_ylim(-0.9, len(E) - 0.1)
fx.spines['left'].set_bounds(y2.min(), y2.max())
for y, r in zip(y2, E.itertuples()):
    fx.text(r.electrophilicity_index_eV + wm * 0.012, y,
            f'{r.electrophilicity_index_eV:.1f}', fontsize=FS_FOOT - 1,
            va='center', ha='left', color='#555555')
gm2 = E.compound_name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])
lo, hi = int(y2[gm2].min()), int(y2[gm2].max())
# 2026-10-20. 종전에는 1.13*wm (약 45 eV) 에 그려, 막대가 6 eV 에서 끝나는
# 세 행 옆에 아무 데도 붙지 않은 붉은 선처럼 보였다. 값 라벨 바로 뒤로 옮긴다.
fx.annotate('', xy=(wm * 0.215, lo - 0.45), xytext=(wm * 0.215, hi + 0.45),
            arrowprops=dict(arrowstyle='-', color=C_REF, lw=2.0))
fx.text(wm * 0.235, (lo + hi) / 2, 'the three\nG agents', fontsize=FS_ANNOT - 2,
        color=C_REF, va='center', ha='left')
fx.text(wm * 1.16, len(E) - 1.0,
        r'$\omega$ is reported here and enters no'
        '\ndistance: it is an exact function of the\n'
        'two orbital energies. The three G-series\n'
        'agents still rank lowest of all 11 agents;\n'
        'GA ranks 8th and VX 7th. 9 of 14\n'
        'interferents score higher than sarin.',
        fontsize=FS_ANNOT - 2, color=C_REF, va='top', ha='left',
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.3', alpha=0.95))
fx.set_xlabel(r'Parr electrophilicity index, $\omega = \mu^{2}/2\eta$ (eV)'
              '\n(derived, not a coordinate of any space)')
fx.legend(handles=[Line2D([], [], marker='s', ls='', ms=8.9, mfc=c, mec='white',
                          label=l) for c, l in
                   [(C_REAL, 'authentic agent'), (C_SIM, 'incumbent simulant'),
                    (C_INTF, 'interferent / solvent')]],
          loc='lower right', frameon=True, fontsize=FS_LEG - 3,
          bbox_to_anchor=(1.01, 0.01), framealpha=0.95, edgecolor='#CCCCCC',
          borderpad=0.35, handletextpad=0.4, labelspacing=0.3)
dress(fx); panel_label(fx, 'b')

# =====================================================================
# Figure 1 (a)(b) 와 Figure 2 (c)(d). 2026-10-20 에 make_fig_panel.py 를
# 흡수했다. 본문 도판을 세 장 x 네 패널로 맞추면서 스크립트를 하나로 합쳤다.
# =====================================================================
M = json.load(open('outputs/match_by_panel.json'))
B = json.load(open('outputs/stat_supplement2.json'))
PF = json.load(open('outputs/pesticide_family.json'))
SUB = json.load(open('outputs/subspace_pvalues.json'))
SS = json.load(open('outputs/stat_supplement.json'))

_inc = {a['agent']: a for a in M['panel35_incumbent3']['agents']}
_allc = {a['agent']: a for a in M['panel35_incumbent_plus_lit8']['agents']}
_order = sorted(_inc, key=lambda a: _allc[a]['d_norm'])

# ---- Fig 1 (a) 작용제별 최근접 거리: 기존 3종 대 전체 후보 (35종 패널)
_y = np.arange(len(_order))[::-1]
DY = 0.17   # 두 마커가 같은 값이면 겹쳐 하나로 보이므로 행 안에서 떼어 놓는다
for yy, a in zip(_y, _order):
    d_inc, d_all = _inc[a]['d_norm'], _allc[a]['d_norm']
    gx.plot([d_all, d_inc], [yy - DY, yy + DY], color='#BBBBBB', lw=1.3, zorder=1)
    gx.scatter([d_inc], [yy + DY], s=78, c=C_SIM, edgecolors='white', lw=1.0,
               zorder=3)
    gx.scatter([d_all], [yy - DY], s=78, c=C_REAL, edgecolors='white', lw=1.0,
               marker='D', zorder=3)
    gx.text(d_inc + 0.085, yy + DY, s(_inc[a]['nearest_simulant']),
            fontsize=FS_ANNOT, va='center', ha='left', color=C_SIM)
    gx.text(d_all - 0.085, yy - DY, s(_allc[a]['nearest_simulant']),
            fontsize=FS_ANNOT, va='center', ha='right', color=C_REAL)
gx.axvline(2.0, color=C_REF, ls='--', lw=1.4, zorder=2)
gx.text(2.04, -1.28, 'match threshold', color=C_REF, fontsize=FS_ANNOT - 1,
        va='bottom', ha='left')
gx.set_yticks(_y); gx.set_yticklabels([s(a) for a in _order])
gx.set_xlim(0.18, 3.95); gx.set_ylim(-1.35, len(_order) + 1.75)
gx.set_xlabel('distance to nearest candidate, as a multiple of the typical '
              'nearest-neighbour separation of the merged space')
gx.legend(handles=[
    Line2D([], [], marker='o', ls='', ms=8, mfc=C_SIM, mec='white', mew=1.0,
           label='nearest of the three incumbent simulants'),
    Line2D([], [], marker='D', ls='', ms=7.5, mfc=C_REAL, mec='white', mew=1.0,
           label="nearest once the field's other five simulants are in")],
    loc='upper right', frameon=False, bbox_to_anchor=(1.0, 1.0),
    fontsize=FS_LEG - 1, handletextpad=0.4, labelspacing=0.3, borderaxespad=0.3)
dress(gx); panel_label(gx, 'a')

# ---- Fig 1 (b) 서술자 부트스트랩: 유효매칭 수 분포
BT = B['distance_bootstrap']
KK = [('incumbent3', 'the three incumbent\nsimulants', C_SIM),
      ('incumbent_plus_lit8', "the field's eight\nsimulants", C_REAL)]
_pos = np.arange(len(KK))[::-1]
for pos, (key, lab, col) in zip(_pos, KK):
    o = BT[key]; lo, hi = o['ci']
    hx.plot([lo, hi], [pos, pos], color=col, lw=5, solid_capstyle='butt',
            alpha=0.33, zorder=2)
    hx.scatter([o['mean']], [pos], s=110, c=col, edgecolors='white', lw=1.2,
               zorder=4)
    hx.scatter([o['median']], [pos], s=90, marker='|', c=C_INK, lw=2.0, zorder=5)
    hx.text(hi + 0.22, pos, f"median {o['median']}, 95% CI {lo} to {hi}",
            fontsize=FS_ANNOT, va='center', ha='left', color=col)
hx.set_yticks(_pos); hx.set_yticklabels([l for _, l, _ in KK])
hx.set_xlim(-0.5, 13.6); hx.set_ylim(-1.05, len(KK) - 0.40)
hx.set_xlabel('agents of eleven with a valid match, over 1000 descriptor '
              'bootstrap replicates')
_d = BT['difference']
hx.text(5.4, -0.76,
        f"difference positive in {100*_d['frac_positive']:.0f}% of replicates "
        f"(median {_d['median']}, 95% CI {_d['ci'][0]} to {_d['ci'][1]})",
        fontsize=FS_ANNOT, va='center', ha='left', color=C_REF,
        bbox=dict(fc='white', ec=C_REF, boxstyle='round,pad=0.32', alpha=0.95))
dress(hx); panel_label(hx, 'b')

# ---- Fig 2 (c) 무작위 3-서술자 부분공간 2000회에서의 p 분포
_BINS = np.logspace(np.log10(1.3e-3), 0.0, 34)
for _k, _c, _lab in [('raw', '#9AA7B4', 'as measured'),
                     ('whitened', '#6B8E23', 'after ZCA whitening')]:
    hv = np.asarray(SUB[_k], float)
    sx.hist(hv, bins=_BINS, color=_c, alpha=0.52, label=_lab,
            edgecolor='white', lw=0.5, zorder=2)
_obs = SUB['observed']['struct_full']; _pc = SUB['percentiles']['struct_full']
sx.axvline(_obs, color=C_REAL, lw=1.8, zorder=4)
sx.axvline(0.05, color=C_REF, ls=(0, (5, 3)), lw=1.3, zorder=3)
sx.set_xscale('log')
sx.set_xlim(1.3e-3, 1.0); sx.set_ylim(0, 318)
# 주석은 패널 안에 머물러야 한다. 오른쪽으로 두면 옆 패널의 축 라벨을 덮고,
# 위로 두면 범례와 겹친다. 범례 아래 빈 띠에 왼쪽 정렬로 세 줄을 넣는다.
sx.annotate(f'the structural space\nused here gives $p$ = {_obs:.3f},\n'
            f'the {_pc:.0f}th percentile of these draws',
            xy=(_obs, 112), xytext=(1.45e-3, 248), fontsize=FS_ANNOT - 1,
            color=C_REAL, ha='left', va='top',
            arrowprops=dict(arrowstyle='-', color=C_REAL, lw=0.9))
sx.text(0.0565, 3, r'$\alpha$ = 0.05', color=C_REF, fontsize=FS_ANNOT - 1,
        rotation=90, va='bottom', ha='left')
sx.set_xlabel('smallest $p$ of the six statistics,\n'
              'random three-descriptor subspace')
sx.set_ylabel('draws')
sx.legend(loc='upper left', frameon=False, fontsize=FS_LEG - 2,
          handletextpad=0.5, labelspacing=0.28, borderaxespad=0.25)
dress(sx); panel_label(sx, 'c')

# 분류 색 범례는 캡션에 적는다. 28개 점이 패널 폭을 다 쓴다.

# ---- Fig 2 (d) 보정 5체계 x 패널 4구성
CFG = [('panel28_incumbent3', '28 cpd.\ninc. 3'),
       ('panel35_incumbent3', '35 cpd.\ninc. 3'),
       ('panel35_incumbent_plus_lit8', '35 cpd.\n+ lit. 8'),
       ('PESTICIDE', '35 cpd.\nrelab. 10')]
SCH = [('min_p', 'smallest nominal $p$', 'o'),
       ('bh_all_min_q', 'BH over 18', 's'),
       ('bh_by_space_min_q', 'BH over 6, per space', '^'),
       ('wy', 'Westfall-Young min$P$', 'v'),
       ('prespecified_bh_min_q', 'statistic fixed in advance', 'D')]
COLS = [C_INK, '#7F7F7F', '#A9A9A9', '#6B8E23', C_REAL]
PFMAP = {'min_p': 'min_p', 'bh_all_min_q': 'bh_all_min_q',
         'bh_by_space_min_q': 'bh_by_space_min_q',
         'wy': 'westfall_young_adjusted_p',
         'prespecified_bh_min_q': 'prespecified_bh_min_q'}
WY = {k: v['wy_adjusted_p'] for k, v in B['westfall_young'].items()}
_xs = np.arange(len(CFG))
for (key, lab, mk), col, off in zip(SCH, COLS,
                                    np.linspace(-0.26, 0.26, len(SCH))):
    v = [PF[PFMAP[key]] if c == 'PESTICIDE' else
         (WY[c] if key == 'wy' else SS['family'][c][key]) for c, _ in CFG]
    ix.scatter(_xs + off, v, s=78, marker=mk, c=col, edgecolors='white',
               lw=1.0, zorder=3, label=lab)
ix.axhline(0.05, color=C_REF, ls='--', lw=1.4, zorder=2)
ix.text(3.50, 0.062, r'$\alpha$ = 0.05', color=C_REF, fontsize=FS_ANNOT - 1,
        va='bottom', ha='right')
ix.set_yscale('log'); ix.set_ylim(3e-5, 2.6)
ix.set_xticks(_xs); ix.set_xticklabels([l for _, l in CFG], fontsize=FS_TICK - 2)
ix.set_xlim(-0.55, 3.55)
ix.set_ylabel('smallest adjusted value\nover the 18 tests')
ix.legend(loc='lower left', frameon=True, ncol=1, fontsize=FS_LEG - 3,
          bbox_to_anchor=(-0.012, -0.012), handletextpad=0.4, labelspacing=0.28,
          framealpha=0.95, edgecolor='#CCCCCC', borderpad=0.35)
dress(ix); panel_label(ix, 'd')

_n0 = M['panel35_incumbent3']['n_match']
_n8 = M['panel35_incumbent_plus_lit8']['n_match']
_n28 = M['panel28_incumbent3']['n_match']
fig_annot.record('Fig_main', 'panel_a',
                 f'valid matches: {_n28} of 11 on the 28-compound panel, '
                 f'{_n0} of 11 on the 35-compound panel with the same three '
                 f'simulants, {_n8} of 11 once the five literature simulants '
                 f'are admitted as simulants.')
fig_annot.record('Fig_main', 'panel_b',
                 f"incumbent three match a median of "
                 f"{BT['incumbent3']['median']} of 11 agents "
                 f"(95% CI {BT['incumbent3']['ci'][0]} to "
                 f"{BT['incumbent3']['ci'][1]}); the field's eight match a "
                 f"median of {BT['incumbent_plus_lit8']['median']} "
                 f"(95% CI {BT['incumbent_plus_lit8']['ci'][0]} to "
                 f"{BT['incumbent_plus_lit8']['ci'][1]}).")
fig_annot.record('Fig_spaces', 'panel_c',
                 f'the structural space used here gives p = {_obs:.3f}, the '
                 f'{_pc:.1f}th percentile of 2000 random three-descriptor '
                 f'subspaces; {100*(np.asarray(SUB["raw"], float) < 0.05).mean():.0f}% '
                 f'of those draws fall below 0.05 as measured and '
                 f'{100*(np.asarray(SUB["whitened"], float) < 0.05).mean():.0f}% '
                 f'after whitening.')
fig_annot.record('Fig_spaces', 'panel_d',
                 'on the primary panel no scheme falls below 0.05; on the '
                 '28-compound panel only the pre-specified-statistic scheme '
                 'does; the fourth configuration, which relabels the two '
                 'pesticides as simulants, falls below 0.05 under every '
                 'scheme.')


fig.savefig('outputs/Fig_main.png', dpi=200, bbox_inches='tight',
            facecolor='white')
fig2.savefig('outputs/Fig_spaces.png', dpi=200, bbox_inches='tight',
             facecolor='white')
figS.savefig('outputs/Fig_plocal.png', dpi=200, bbox_inches='tight',
             facecolor='white')
print('[fig] outputs/Fig_spaces.png')
print('[fig] outputs/Fig_plocal.png')
print('[fig] outputs/Fig_main.png')
