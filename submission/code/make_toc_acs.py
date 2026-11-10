#!/usr/bin/env python3
"""ACS Sensors TOC 도판 (2026-11-02 신설).

규격은 ACS Sensors 저자 안내에서 실측했다: 8.25 cm × 4.45 cm
(3.25 in × 1.75 in), 컬러 300 dpi 이상, 최종 크기에서 글자 4.5 pt 이상.
따라서 figsize 를 실척으로 잡고 폰트를 포인트 그대로 쓴다. 확대·축소하면
글자 크기 보증이 깨지므로 bbox_inches='tight' 를 쓰지 않는다.

구판(make_ga.py)은 17.5 in 폭에 JHM 프레임(구조공간 대 전자공간, 독성축)
이었다. ACS 판본의 주장은 "서술자 유사도가 어레이 응답을 예측하지 못한다"
이므로 내용을 바꾼다.

수치는 전부 outputs/deltaE_premise.json 에서 읽는다. 손으로 적지 않는다.
"""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker

# 워밍 머티드 어스톤 팔레트 (2026-04-13 Ku 지시)
C_TEAL = '#4AACB0'
C_DARK_TEAL = '#2E7D7A'
C_RED = '#C94F4A'
C_SAGE = '#7FAF6A'
C_BROWN = '#8B6347'
C_TAN = '#C9A882'
C_GREY = '#8A8279'

# ------------------------------------------------------------------ 데이터
# 2026-11-10 (3층 재설계). 종전 TOC 는 2층(배열 응답 전제)을 그렸다. 그
# 분석이 전부 SI 로 내려갔으므로 독자가 가장 먼저 보는 그림이 본문 논지와
# 맞지 않게 됐다. 본문은 1층(문서 기록)과 3층(분석자 자유도)이다. 그래서
# (a) 패널 구성에 따른 유효 매칭 수, (b) 구성에 따른 보정 후 q 로 바꾼다.
# 수치는 전부 산출물에서 읽는다. 손으로 적지 않는다.
MB = json.load(open('outputs/match_by_panel.json'))
FAM = json.load(open('outputs/stat_supplement.json'))['family']
PF = json.load(open('outputs/pesticide_family.json'))

# (a) 유효 매칭 수 / 11. 세 패널 구성.
A_KEYS = [('panel35_incumbent3', 'full\npanel'),
          ('panel28_incumbent3', 'omit\nseven'),
          ('panel35_incumbent_plus_lit8', 'eight\ndesig.')]
MATCH = [MB[k]['n_match'] for k, _ in A_KEYS]
A_LAB = [lab for _, lab in A_KEYS]
N_AGENT = len(MB['panel35_incumbent3']['agents'])

# (b) 18개 검정에 대한 Benjamini-Hochberg 최소 q. 네 구성.
B_KEYS = [('panel35_incumbent3', 'incumbent 3, full panel'),
          ('panel35_incumbent_plus_lit8', 'eight designated'),
          ('panel28_incumbent3', 'incumbent 3, seven omitted')]
QS = [FAM[k]['bh_all_min_q'] for k, _ in B_KEYS] + [PF['bh_all_min_q']]
B_LAB = [lab for _, lab in B_KEYS] + ['two labels moved']

FS = 5.6          # 본문 글자. 최종 크기에서 4.5 pt 하한을 넉넉히 넘긴다.
FS_SM = 5.0
FS_TITLE = 6.4

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': FS,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 1.8,
    'ytick.major.size': 1.8,
    'xtick.labelsize': FS_SM,
    'ytick.labelsize': FS_SM,
})

fig = plt.figure(figsize=(3.25, 1.75), dpi=600)

# 머리글 한 줄. 높이 예산 0.17 in.
# 2026-11-06. GATE 8 9차. 종전 머리글은 'Descriptor similarity does not
# predict array response' 였다. 제목은 검정력 한계 때문에 'does not predict'
# 에서 'no resolvable agreement' 로 이미 내렸는데 TOC 만 구판이었다. 독자가
# 가장 먼저 보는 자리에 철회한 주장이 남아 있었으므로 제목과 동기화한다.
# 2026-11-07 (방향 A). 머리글을 제목과 함께 측정 소견으로 옮긴다. 종전
# 머리글은 귀무였고, 그것이 세 회차 연속 편집장 지적의 핵심이었다.
# 폭 실측: FS_TITLE 6.4 에서 2.72 in (그림 3.25 in, 상한 3.10 in).
# 2026-11-07 (12차 Step 5). 방향 A 머리글 철회. 폭 실측은 아래 출력에서
# 확인한다 (상한 3.10 in).
# 2026-11-10 (3층 재설계). 종전 머리글 '13 of 34 cases are nearest a
# different compound' 은 2층 수치였고 이미 철회된 읽기였다(관측 13 이 귀무
# 29.8 보다 낮으므로 그 수만으로는 소견을 과장한다). 본문 결과 절의 제목과
# 같은 문장으로 둔다. 수치를 넣지 않으므로 구판이 될 수 없다.
# 폭 실측은 아래 출력에서 확인한다 (상한 3.10 in).
HEADLINE = 'The comparison decides the answer'
fig.text(0.5, 0.955, HEADLINE,
         ha='center', va='top', fontsize=FS_TITLE, fontweight='bold',
         color='#2B2B2B')

# 패널 배치. (a) 세로 막대는 좁게, (b) 는 가로 막대로 둬서 범주 이름을
# 축 왼쪽에 한 줄로 적는다. 2026-11-10: 처음에 (b) 도 세로로 그렸더니 네
# 범주 이름이 서로 겹쳤다. 렌더해서 눈으로 보고 가로로 바꿨다.
axA = fig.add_axes([0.098, 0.225, 0.250, 0.545])
axB = fig.add_axes([0.605, 0.225, 0.370, 0.545])

# ---------------------------------------------------------------- (a)
_cols = [C_RED, C_TAN, C_DARK_TEAL]
axA.bar(range(len(MATCH)), MATCH, width=0.62, color=_cols,
        edgecolor='none', zorder=3)
# 값 라벨. 0 은 막대가 없으므로 축 바로 위에, 나머지는 막대 안쪽 상단에.
for i, v in enumerate(MATCH):
    if v == 0:
        axA.text(i, 0.20, '0', ha='center', va='bottom',
                 fontsize=5.4, color=C_RED, fontweight='bold')
    else:
        axA.text(i, v - 0.22, f'{v}', ha='center', va='top',
                 fontsize=5.4, color='white', fontweight='bold')
axA.set_xticks(range(len(MATCH)))
axA.set_xticklabels(A_LAB, fontsize=4.3, linespacing=1.12)
axA.tick_params(axis='x', pad=1.4, length=0)
axA.set_xlim(-0.60, len(MATCH) - 0.40)
axA.set_ylim(0, 6.3)
axA.set_yticks([0, 2, 4, 6])
axA.set_ylabel(f'agents matched of {N_AGENT}', fontsize=FS_SM, labelpad=1.5)
axA.grid(axis='y', color='#DDD9D2', lw=0.35, zorder=0)
axA.set_axisbelow(True)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)
axA.set_title('valid matches', fontsize=FS_SM, pad=2.2, color='#2B2B2B')

# ---------------------------------------------------------------- (b)
_colsB = [C_RED, C_TAN, C_BROWN, C_DARK_TEAL]
_y = np.arange(len(QS))[::-1]          # 위에서 아래로 입력 순서대로
axB.barh(_y, QS, height=0.60, color=_colsB, edgecolor='none', zorder=3)
axB.set_xscale('log')
axB.axvline(0.05, color='#555555', lw=0.6, ls=(0, (2.2, 1.6)), zorder=4)
# 값 라벨은 막대 오른쪽 끝 바깥에 둔다 (우측이 1순위 자리, 겹침 없음).
for yy, v in zip(_y, QS):
    _s = f'{v:.3f}' if v >= 0.001 else f'{v:.4f}'
    axB.text(v * 1.45, yy, _s, ha='left', va='center',
             fontsize=4.5, color='#2B2B2B')
axB.set_yticks(_y)
axB.set_yticklabels([l.replace('\n', ' ') for l in B_LAB], fontsize=4.3)
axB.tick_params(axis='y', pad=1.2, length=0)
axB.set_ylim(-0.75, len(QS) - 0.25)
axB.set_xlim(3.0e-4, 6.0)
axB.set_xticks([1e-3, 1e-2, 1e-1, 1e0])
axB.set_xticklabels(['0.001', '0.01', '0.1', '1'], fontsize=4.5)
axB.xaxis.set_minor_locator(mticker.NullLocator())
axB.set_xlabel('smallest $q$ over 18 tests', fontsize=FS_SM, labelpad=1.5)
# 역치 설명은 파선 바로 위, 맨 위 막대 위의 빈 띠에 둔다.
axB.text(0.056, len(QS) - 0.52, '$q$ = 0.05', ha='left', va='center',
         fontsize=4.4, color='#2B2B2B', style='italic')
axB.grid(axis='x', color='#DDD9D2', lw=0.35, zorder=0)
axB.set_axisbelow(True)
for s in ('top', 'right', 'left'):
    axB.spines[s].set_visible(False)
axB.set_title('after correction', fontsize=FS_SM, pad=2.2, color='#2B2B2B')

# 2026-11-07 (11차). TOC 주석을 등재부에 기록한다.
#
# 왜. 9차에서 TOC 머리글이 본문이 이미 철회한 주장('does not predict')을
# 들고 있었고, 독자가 가장 먼저 보는 자리인데 어느 검사도 보지 않았다.
# 그때는 눈으로 발견했다. TOC 자신의 글자를 등재해 앵커를 이 도판에 둔다.
import fig_annot

fig_annot.record('toc_acs', 'headline', HEADLINE)
fig_annot.record('toc_acs', 'panel_a',
                 'valid matches of %d agents: %s'
                 % (N_AGENT, ', '.join(
                     '%d (%s)' % (m, l.replace('\n', ' '))
                     for m, l in zip(MATCH, A_LAB))))
fig_annot.record('toc_acs', 'panel_b',
                 'smallest Benjamini-Hochberg q over 18 tests: %s'
                 % ', '.join('%.4f (%s)' % (q, l.replace('\n', ' '))
                             for q, l in zip(QS, B_LAB)))

fig.savefig('outputs/toc_acs.png', dpi=600, facecolor='white')
fig.savefig('outputs/toc_acs.pdf', facecolor='white')
print('[TOC] outputs/toc_acs.png / .pdf 생성')
print(f'[TOC] 실척 3.25 x 1.75 in, 600 dpi -> {int(3.25*600)} x {int(1.75*600)} px')
print(f'[TOC] (a) 유효 매칭 {MATCH} / {N_AGENT}')
print(f'[TOC] (b) 최소 q ' + ', '.join(f'{q:.4f}' for q in QS))

# 머리글 실측 폭. 상한 3.10 in 을 넘으면 즉시 알린다.
_r = fig.canvas.get_renderer() if hasattr(fig.canvas, 'get_renderer') else None
if _r is None:
    fig.canvas.draw()
    _r = fig.canvas.get_renderer()
_t = fig.texts[0]
_w = _t.get_window_extent(renderer=_r).width / fig.dpi
print(f'[TOC] 머리글 실측 폭 {_w:.2f} in (상한 3.10)')
if _w > 3.10:
    raise SystemExit(f'[TOC] 머리글이 그림 폭을 넘는다: {_w:.2f} in')
