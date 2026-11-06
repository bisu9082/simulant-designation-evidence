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

# --- (a) DMMP 대 작용제 코사인, UV 채널 -------------------------------------
uv = next(r for r in D['dmmp'] if r['mode'] == 'uv')
AG = ['GA', 'GB', 'GD', 'GF', 'VX']
cos = [uv[f'cos_{a}'] for a in AG]

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
# 폭 실측: FS_TITLE 6.4 에서 2.98 in (그림 3.25 in, 여백 포함 상한 3.10 in).
# 제목 전문은 3.73 in 로 넘쳐서 'distance' 와 'array' 를 뺐다. 두 낱말은
# 좌우 패널의 축 설명이 이미 적고 있다.
fig.text(0.5, 0.955, 'No resolvable agreement between descriptor and response',
         ha='center', va='top', fontsize=FS_TITLE, fontweight='bold',
         color='#2B2B2B')

# 패널 배치. 좌 52%, 우 48%. 라벨 겹침을 피하려고 여백을 명시한다.
axA = fig.add_axes([0.105, 0.195, 0.355, 0.60])
axB = fig.add_axes([0.620, 0.195, 0.350, 0.60])

# ---------------------------------------------------------------- (a)
colA = [C_DARK_TEAL if c > 0 else C_RED for c in cos]
axA.bar(range(5), cos, width=0.68, color=colA, edgecolor='none', zorder=3)
axA.axhline(0, color='#555555', lw=0.6, zorder=4)
axA.set_xticks(range(5))
axA.set_xticklabels(AG, fontsize=FS_SM)
axA.set_ylim(-1.05, 1.00)
axA.set_yticks([-0.8, -0.4, 0.0, 0.4, 0.8])
axA.set_ylabel('cosine to DMMP', fontsize=FS_SM, labelpad=1.5)
axA.grid(axis='y', color='#DDD9D2', lw=0.35, zorder=0)
axA.set_axisbelow(True)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)
# 값 라벨은 막대 바깥쪽에 둔다 (막대 위에 얹으면 겹친다).
for i, c in enumerate(cos):
    off = 0.055 if c > 0 else -0.055
    axA.text(i, c + off, f'{c:+.2f}', ha='center',
             va='bottom' if c > 0 else 'top',
             fontsize=4.8, color=colA[i], fontweight='bold')
axA.set_title('the standard G-agent simulant,\nmeasured beside the agents',
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

fig.savefig('outputs/toc_acs.png', dpi=600, facecolor='white')
fig.savefig('outputs/toc_acs.pdf', facecolor='white')
print('[TOC] outputs/toc_acs.png / .pdf 생성')
print(f'[TOC] 실척 3.25 x 1.75 in, 600 dpi -> {int(3.25*600)} x {int(1.75*600)} px')
print(f'[TOC] (a) UV 코사인 {dict(zip(AG, [round(c,3) for c in cos]))}')
print(f'[TOC] (b) rho 범위 {min(rho):+.3f} ~ {max(rho):+.3f}, n={len(rho)}')
