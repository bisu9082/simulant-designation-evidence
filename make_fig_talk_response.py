#!/usr/bin/env python3
"""발표용 응답공간 도판 (2026-10-28 신설).

왜. ACS Sensors 프레임의 헤드라인은 "어레이는 화합물을 구별하는데 서술자
거리가 응답 거리로 전이되지 않는다" 이다. 그런데 원고의 도판 다섯 장은 전부
서술자공간과 독성축이고 응답공간 도판이 없다. 발표에서 이 축을 앞에 세우려면
그림이 필요하다.

세 패널.
 (a) 판별력. 물질내 대 물질간 코사인, 여섯 셀. 어레이가 작동한다는 전제.
 (b) 전이. Mantel 전수치환 12셀의 BH q. 0.05 아래가 없다.
 (c) DMMP 대 작용제 코사인을 농도별로. 표준 G 모사체가 GD, GF 와 반상관.

출력: outputs/Fig_talk_response.png
"""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# feedback_figure_color_palette 의 워밍 머티드 어스톤 팔레트
C_TEAL = '#4AACB0'
C_DARK_TEAL = '#2E7D7A'
C_RED = '#C94F4A'
C_AMBER = '#E8943A'
C_SAGE = '#7FAF6A'
C_BROWN = '#8B6347'
C_TAN = '#C9A882'
C_DUSTY = '#5B8DB8'
C_GOLD = '#F0C040'

plt.rcParams.update({
    'font.size': 11, 'axes.labelsize': 11.5, 'axes.titlesize': 12,
    'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 9.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.9, 'figure.dpi': 200,
})


def main():
    pre = json.load(open('outputs/deltaE_premise.json'))
    roll = json.load(open('outputs/rollback_c.json'))['deltaE']
    ext = json.load(open('outputs/deltaE_extended.json'))

    fig = plt.figure(figsize=(14.5, 4.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25],
                          wspace=0.34, left=0.055, right=0.985,
                          bottom=0.17, top=0.86)

    # ---------------------------------------------------------- (a) 판별력
    ax = fig.add_subplot(gs[0, 0])
    disc = pre['discrimination']
    lab = [f"{d['dataset']}\n{d['mode']}" for d in disc]
    y = np.arange(len(disc))
    w = 0.36
    ax.barh(y + w / 2, [d['within'] for d in disc], w, color=C_TEAL,
            edgecolor='white', lw=0.6)
    ax.barh(y - w / 2, [d['between'] for d in disc], w, color=C_TAN,
            edgecolor='white', lw=0.6)
    ax.axvline(0, color='0.35', lw=0.9)
    # p 값을 한 열에 정렬한다. 막대 끝에 붙이면 열이 들쭉날쭉해 읽기 어렵다.
    for i, d in enumerate(disc):
        ax.text(0.60, i, f"p = {d['p']:.1e}", va='center', ha='left',
                fontsize=8.6, color=C_DARK_TEAL)
    # 색 키를 왼쪽 빈 영역에 둔다. 윗행에 직접 붙였더니 축 눈금과 겹쳤고,
    # 범례 상자로 두었더니 p 열을 가렸다. 여기는 어느 막대도 닿지 않는다.
    ax.text(-0.70, 1.05, 'same compound,\nreplicates', va='center',
            ha='left', fontsize=8.6, color=C_DARK_TEAL, fontweight='bold',
            linespacing=1.25)
    ax.text(-0.70, 0.05, 'different\ncompounds', va='center', ha='left',
            fontsize=8.6, color=C_BROWN, fontweight='bold', linespacing=1.25)
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=9)
    ax.set_xlim(-0.72, 1.02)
    ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax.set_xlabel('mean cosine between response patterns')
    ax.set_title('(a)  the array separates compounds',
                 loc='left', fontweight='bold')

    # ---------------------------------------------------------- (b) 전이
    ax = fig.add_subplot(gs[0, 1])
    cells = roll['mantel']
    ps = sorted(c['p_mantel'] for c in cells)
    m = len(ps)
    q = [min(1.0, p * m / (i + 1)) for i, p in enumerate(ps)]
    for i in range(m - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    xs = np.arange(1, m + 1)
    ax.scatter(xs, q, s=52, color=C_BROWN, zorder=3,
               edgecolor='white', lw=0.7, label='12 Mantel cells (BH $q$)')
    fis = ext['combined']['cells']
    fq = sorted(r['bh_q'] for r in fis)
    ax.scatter(np.linspace(2.5, 9.5, len(fq)), fq, s=52, marker='D',
               color=C_AMBER, zorder=3, edgecolor='white', lw=0.7,
               label='6 Fisher combinations (BH $q$)')
    ax.axhline(0.05, color=C_RED, ls='--', lw=1.3)
    ax.text(12.6, 0.068, r'$q=0.05$', color=C_RED, fontsize=9.5,
            ha='right', va='bottom')
    ax.text(6.5, 0.015, 'no test reaches it', color=C_RED, fontsize=9.5,
            ha='center', va='bottom', fontweight='bold')
    ax.set_ylim(0, 1.03)
    ax.set_xlim(0.2, 12.8)
    ax.set_xticks([])
    ax.set_ylabel('adjusted $p$')
    ax.set_xlabel('tests of descriptor distance against response distance')
    ax.set_title('(b)  but descriptor distance does not transfer',
                 loc='left', fontweight='bold')
    ax.legend(loc='upper left', frameon=False, handlelength=1.1)

    # ---------------------------------------------------------- (c) DMMP
    ax = fig.add_subplot(gs[0, 2])
    conc = ext['concentration']
    agents = ['GA', 'GB', 'GD', 'GF', 'VX']
    col = {'GA': C_AMBER, 'GB': C_RED, 'GD': C_DARK_TEAL,
           'GF': C_DUSTY, 'VX': C_BROWN}
    # 반사 채널은 뺀다. 열 줄이 겹쳐 15초 안에 읽히지 않는다. 부록으로 돌린다.
    for ag in agents:
        r = next((x for x in conc
                  if x['mode'] == 'uv' and x['agent'] == ag), None)
        if r is None:
            continue
        ax.plot(r['conc'], r['cosines'], '-', marker='o', ms=5.0, lw=1.9,
                color=col[ag])
    ax.axhline(0, color='0.35', lw=0.9)
    ax.axhspan(-0.85, 0, color=C_RED, alpha=0.055, lw=0)
    ax.text(9.5, -0.80, 'anti-correlated with the agent', fontsize=9,
            color=C_RED, va='bottom', ha='left')
    ax.set_xscale('log')
    ax.set_xticks([10, 50, 100, 500])
    ax.set_xticklabels(['10', '50', '100', '500'])
    ax.set_xlim(8, 1150)
    ax.set_ylim(-0.85, 0.95)
    ax.set_xlabel('concentration (ppm)')
    ax.set_ylabel('cosine, DMMP pattern vs agent pattern')
    ax.set_title('(c)  the standard G-agent simulant,\n'
                 '      against the G agents (UV channel)',
                 loc='left', fontweight='bold')
    # GD 와 GF 의 500 ppm 값이 0.028 밖에 안 떨어져 라벨이 겹친다. 수직으로 민다.
    dy = {'GA': 0, 'GB': 0, 'GD': -9, 'GF': +9, 'VX': 0}
    for ag in agents:
        r = next((x for x in conc
                  if x['mode'] == 'uv' and x['agent'] == ag), None)
        if r:
            ax.annotate(ag, (r['conc'][-1], r['cosines'][-1]),
                        xytext=(7, dy[ag]), textcoords='offset points',
                        color=col[ag], fontsize=10, fontweight='bold',
                        va='center')

    fig.savefig('outputs/Fig_talk_response.png', dpi=200,
                bbox_inches='tight', facecolor='white')
    print('[done] outputs/Fig_talk_response.png')
    return 0


if __name__ == '__main__':
    sys.exit(main())
