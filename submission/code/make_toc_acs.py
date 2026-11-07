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

# 워밍 머티드 어스톤 팔레트 (2026-04-13 Ku 지시)
C_TEAL = '#4AACB0'
C_DARK_TEAL = '#2E7D7A'
C_RED = '#C94F4A'
C_SAGE = '#7FAF6A'
C_BROWN = '#8B6347'
C_TAN = '#C9A882'
C_GREY = '#8A8279'

D = json.load(open('outputs/deltaE_premise.json'))
# 2026-11-07 (12차 Step 5). 패널 (a) 를 교체한다. 종전 패널은 DMMP 대
# 작용제 코사인 다섯 개였고, 중심화가 음의 부호를 만들므로 기저선 없이
# 읽으면 표현의 성질을 화합물의 성질로 보여 준다. 역전 계수로 바꾼다.
R = json.load(open('outputs/reversal_audit.json'))

# --- (a) 자기보다 남이 가까운 사례, 34건 중 ---------------------------------
RT = R['totals']
REV = [RT['reversals_uncentred'], RT['reversals_centred']]
REV_NULL = [RT['null_mean_sum_uncentred'], RT['null_mean_sum_centred']]
N_INST = RT['n_instances']
REV_LAB = ['before\ncentring', 'after\ncentring']

# --- (b) 서술자거리 대 응답거리, 12셀 ---------------------------------------
st = D['space_test']
rho = [r['spearman'] for r in st]
lo = [r['ci_lo'] for r in st]
hi = [r['ci_hi'] for r in st]

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
HEADLINE = '13 of 34 cases are nearest a different compound'
fig.text(0.5, 0.955, HEADLINE,
         ha='center', va='top', fontsize=FS_TITLE, fontweight='bold',
         color='#2B2B2B')

# 패널 배치. 좌 52%, 우 48%. 라벨 겹침을 피하려고 여백을 명시한다.
axA = fig.add_axes([0.105, 0.195, 0.355, 0.60])
axB = fig.add_axes([0.620, 0.195, 0.350, 0.60])

# ---------------------------------------------------------------- (a)
axA.bar(range(2), REV, width=0.52, color=[C_RED, C_DARK_TEAL],
        edgecolor='none', zorder=3)
# 라벨순열 기대값을 열린 표식으로 겹쳐 둔다. 관측이 우연보다 낮다는 것을
# 막대 하나당 점 하나로 보인다.
axA.scatter(range(2), REV_NULL, s=7.0, facecolors='none',
            edgecolors='#2B2B2B', linewidths=0.6, zorder=5)
# 표식 설명은 한 줄로 모은다. 막대마다 붙이면 위쪽 라벨이 패널 제목에
# 닿고 같은 말이 두 번 나온다. 자리는 라벨 우선순위대로 우하단이다.
# 우하단이되 막대 오른쪽 빈 띠 안이다. 우측 정렬로 두면 글자가 막대
# 위로 넘어온다 (막대는 x 1.26 까지다).
axA.text(1.30, 1.2, u'\u25cb = chance', ha='left', va='bottom',
         fontsize=4.4, color='#2B2B2B', style='italic')
axA.set_xticks(range(2))
axA.set_xticklabels(REV_LAB, fontsize=FS_SM)
axA.set_xlim(-0.62, 1.78)
axA.set_ylim(0, N_INST + 2.8)
axA.set_yticks([0, 10, 20, 30])
axA.set_ylabel(f'cases of {N_INST}', fontsize=FS_SM, labelpad=1.5)
axA.grid(axis='y', color='#DDD9D2', lw=0.35, zorder=0)
axA.set_axisbelow(True)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)
for i, v in enumerate(REV):
    axA.text(i, v - 1.1, f'{v}', ha='center', va='top',
             fontsize=5.4, color='white', fontweight='bold')
axA.set_title('nearest neighbour is a\ndifferent compound',
              fontsize=FS_SM, pad=2.0, color='#2B2B2B')

# ---------------------------------------------------------------- (b)
y = np.arange(len(rho))
axB.axvline(0, color='#555555', lw=0.6, zorder=2)
for i in range(len(rho)):
    axB.plot([lo[i], hi[i]], [y[i], y[i]], color=C_TAN, lw=0.9,
             solid_capstyle='butt', zorder=3)
axB.scatter(rho, y, s=5.0, color=C_BROWN, zorder=4, linewidths=0)
axB.set_ylim(-1.9, len(rho) + 0.1)
axB.set_yticks([])
axB.set_xlim(-0.85, 1.0)
axB.set_xticks([-0.5, 0.0, 0.5, 1.0])
axB.set_xlabel('rank correlation, 95% CI', fontsize=FS_SM, labelpad=1.5)
axB.grid(axis='x', color='#DDD9D2', lw=0.35, zorder=0)
axB.set_axisbelow(True)
for s in ('top', 'right', 'left'):
    axB.spines[s].set_visible(False)
axB.set_title('descriptor distance against\nresponse distance, 12 cells',
              fontsize=FS_SM, pad=2.0, color='#2B2B2B')
# 모든 구간이 0 을 지난다는 사실을 한 줄로.
# 2026-11-02. 자리를 두 번 옮겼다. 우하단은 맨 아래 CI 막대 위에 얹혔고,
# 좌상단은 맨 위 CI 막대 시작점(-0.3)에 닿았다. y 하한을 -1.9 로 내려
# 맨 아래 행 밑에 완전히 빈 띠를 만들고 거기 한 줄로 둔다.
axB.text(0.02, 0.015, 'none below $q$ = 0.05', transform=axB.transAxes,
         ha='left', va='bottom', fontsize=4.8, color=C_GREY, style='italic')

# 2026-11-07 (11차). TOC 주석을 등재부에 기록한다.
#
# 왜. 9차에서 TOC 머리글이 본문이 이미 철회한 주장('does not predict')을
# 들고 있었고, 독자가 가장 먼저 보는 자리인데 어느 검사도 보지 않았다.
# 그때는 눈으로 발견했다. 종전에는 구판 make_ga.py 의 graphical_abstract
# 주석만 등재돼 있었고, 그 스크립트는 ACS 빌드가 부르지 않으므로 빈
# 저장소를 클론하면 그 앵커가 통째로 꺼진다(11차 냉클론 시험에서 확인).
# TOC 자신의 글자를 등재해 앵커를 이 도판으로 옮긴다.
import fig_annot

fig_annot.record('toc_acs', 'headline', HEADLINE)
fig_annot.record('toc_acs', 'panel_a',
                 f'nearest neighbour is a different compound in '
                 f'{REV[0]} of {N_INST} cases before centring and '
                 f'{REV[1]} of {N_INST} after; label-permutation chance '
                 f'{REV_NULL[0]:.1f} and {REV_NULL[1]:.1f}')
fig_annot.record('toc_acs', 'panel_b',
                 f'rank correlation, 12 cells: {min(rho):+.3f} to '
                 f'{max(rho):+.3f}; none below q = 0.05')

fig.savefig('outputs/toc_acs.png', dpi=600, facecolor='white')
fig.savefig('outputs/toc_acs.pdf', facecolor='white')
print('[TOC] outputs/toc_acs.png / .pdf 생성')
print(f'[TOC] 실척 3.25 x 1.75 in, 600 dpi -> {int(3.25*600)} x {int(1.75*600)} px')
print(f'[TOC] (a) 역전 {REV[0]}/{N_INST} (중심화 전), {REV[1]}/{N_INST} (후); 귀무 {REV_NULL[0]:.2f} / {REV_NULL[1]:.2f}')
print(f'[TOC] (b) rho 범위 {min(rho):+.3f} ~ {max(rho):+.3f}, n={len(rho)}')
