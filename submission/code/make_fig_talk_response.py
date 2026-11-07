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
    # 2026-11-07 (10차 R1-1). 이 패널이 본문이 철회한 두 통계를 그리고
    # 있었다 — 부호뒤집기 p 와 -1/(n-1) 바닥선. 본문이 지금 의지하는
    # 라벨순열 귀무로 바꾼다.
    pn = {(c['dataset'], c['mode']): c
          for c in json.load(open('outputs/premise_null.json'))['cells']}

    fig = plt.figure(figsize=(14.5, 4.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25],
                          wspace=0.34, left=0.055, right=0.985,
                          bottom=0.17, top=0.86)

    # ---------------------------------------------------------- (a) 판별력
    ax = fig.add_subplot(gs[0, 0])
    disc = pre['discrimination']
    # 2026-11-02 (GATE 8 7차). 본문은 first/second/third run 과
    # reflectance 로 쓰는데 도판만 d1/d2/d3·normal 이었다. 용어를 맞춘다.
    RUN = {'d1': 'first run', 'd2': 'second run', 'd3': 'third run'}
    CH = {'uv': 'ultraviolet', 'normal': 'reflectance'}
    lab = [f"{RUN[d['dataset']]}\n{CH[d['mode']]}" for d in disc]
    y = np.arange(len(disc))
    w = 0.36
    ax.barh(y + w / 2, [d['within'] for d in disc], w, color=C_TEAL,
            edgecolor='white', lw=0.6)
    ax.barh(y - w / 2, [d['between'] for d in disc], w, color=C_TAN,
            edgecolor='white', lw=0.6)
    ax.axvline(0, color='0.35', lw=0.9)
    # 2026-11-07. 종전에는 -1/(n-1) 바닥선을 그렸는데 본문이 그 논거를
    # 철회했다 (재정규화가 영합을 깨므로 바닥이 아니다). 올바른 기준은
    # 라벨순열 귀무가 같은 between 에서 기대하는 within 이다. 그 값을
    # mean_between + null_mean 으로 그린다.
    for i, d in enumerate(disc):
        c = pn[(d['dataset'], d['mode'])]
        exp_within = c['mean_between'] + c['null_mean']
        ax.plot([exp_within, exp_within], [i, i + w], color=C_RED, lw=1.5,
                solid_capstyle='butt', zorder=5)
    ax.plot([], [], color=C_RED, lw=1.5,
            label='within expected under\nlabel permutation')
    # 범례를 축 영역 안 왼쪽아래에 두었더니 x 눈금과 겹쳤다 (렌더 확인).
    # y 하한을 내려 맨 아래 막대 밑에 빈 띠를 만들고 거기 둔다.
    # p 값을 한 열에 정렬한다. 막대 끝에 붙이면 열이 들쭉날쭉해 읽기 어렵다.
    for i, d in enumerate(disc):
        c = pn[(d['dataset'], d['mode'])]
        # 표집 셀은 바닥에 검열돼 있으므로 점추정으로 적지 않는다.
        txt = (f"p = {c['p_one_sided']:.4f}" if c['enumerated']
               else f"p < {1.0 / c['n_null']:.0e}")
        ax.text(0.60, i, txt, va='center', ha='left',
                fontsize=8.6, color=C_DARK_TEAL)
    # 색 키를 왼쪽 빈 영역에 둔다. 윗행에 직접 붙였더니 축 눈금과 겹쳤고,
    # 범례 상자로 두었더니 p 열을 가렸다. 여기는 어느 막대도 닿지 않는다.
    # 2026-11-02 (GATE 8 7차). 'replicates' 는 수행하지 않은 측정이다.
    # 원자료는 (조명, 염료, 작용제, 농도) 조합당 한 행뿐이고 본문도 그렇게
    # 적는다. 같은 화합물의 서로 다른 '농도' 사이 코사인이다.
    ax.text(-0.70, 1.05, 'same compound,\nacross doses', va='center',
            ha='left', fontsize=8.6, color=C_DARK_TEAL, fontweight='bold',
            linespacing=1.25)
    ax.text(-0.70, 0.05, 'different\ncompounds', va='center', ha='left',
            fontsize=8.6, color=C_BROWN, fontweight='bold', linespacing=1.25)
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=9)
    ax.set_xlim(-0.72, 1.02)
    ax.set_ylim(-1.15, len(disc) - 0.45)
    ax.legend(loc='lower left', fontsize=8.0, frameon=False,
              handlelength=1.4, borderpad=0.0, labelspacing=0.2)
    ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax.set_xlabel('mean cosine between response patterns')
    # 제목이 본문의 유보와 충돌했다. 본문: "Most of the apparent separation
    # is therefore a property of the centring". 측정한 것만 적는다.
    # 2026-11-07. 'more than with others' 는 짝별 전칭으로 읽힌다. 본문은
    # 34개 중 13개에서 그것이 거짓임을 적으므로 평균 진술로 맞춘다.
    ax.set_title('(a)  a compound agrees with itself across doses\n'
                 '      more than with the mean of the others',
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
    ax.set_title('(b)  descriptor distance does not resolve it',
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
    # 단위는 µM 이다. 세 원논문 Methods 실측 (2026-11-02). ppm 은 철회.
    ax.set_xlabel(r'concentration ($\mu$M)')
    ax.set_ylabel('cosine, DMMP pattern vs agent pattern')
    ax.set_title('(c)  the standard G-agent simulant,\n'
                 '      against the G agents (UV channel)',
                 loc='left', fontweight='bold')
    # GD 와 GF 의 500 µM 값이 0.028 밖에 안 떨어져 라벨이 겹친다. 수직으로 민다.
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
