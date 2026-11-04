#!/usr/bin/env python3
"""발표 전용 보조 도판 두 장 (2026-10-28 신설).

왜. 원고 도판은 4패널이라 슬라이드 크기에서 글자가 안 읽힌다. 발표에서 쓸
두 장만 따로 그린다. 수치는 전부 정본 산출물에서 읽는다.

A. Fig_talk_dmmp.png   DMMP 대 G 작용제 코사인을 농도별로 (단독 확대)
B. Fig_talk_tox.png    AEGL 실측 공백 + 모델의 일방 실패 2패널

출력: outputs/Fig_talk_dmmp.png, outputs/Fig_talk_tox.png
"""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C_TEAL = '#4AACB0'
C_DARK_TEAL = '#2E7D7A'
C_RED = '#C94F4A'
C_AMBER = '#E8943A'
C_SAGE = '#7FAF6A'
C_BROWN = '#8B6347'
C_TAN = '#C9A882'
C_DUSTY = '#5B8DB8'

plt.rcParams.update({
    'font.size': 13, 'axes.labelsize': 13.5, 'axes.titlesize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 11.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 1.0, 'figure.dpi': 200,
})


def fig_dmmp():
    ext = json.load(open('outputs/deltaE_extended.json'))
    conc = ext['concentration']
    agents = ['GA', 'GB', 'GD', 'GF', 'VX']
    col = {'GA': C_AMBER, 'GB': C_RED, 'GD': C_DARK_TEAL,
           'GF': C_DUSTY, 'VX': C_BROWN}
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.axhspan(-0.9, 0, color=C_RED, alpha=0.06, lw=0)
    for ag in agents:
        r = next((x for x in conc
                  if x['mode'] == 'uv' and x['agent'] == ag), None)
        if r is None:
            continue
        ax.plot(r['conc'], r['cosines'], '-', marker='o', ms=7.5, lw=2.6,
                color=col[ag])
    ax.axhline(0, color='0.3', lw=1.1)
    dy = {'GA': 0, 'GB': 0, 'GD': -11, 'GF': +11, 'VX': 0}
    for ag in agents:
        r = next((x for x in conc
                  if x['mode'] == 'uv' and x['agent'] == ag), None)
        if r:
            ax.annotate(ag, (r['conc'][-1], r['cosines'][-1]),
                        xytext=(10, dy[ag]), textcoords='offset points',
                        color=col[ag], fontsize=14, fontweight='bold',
                        va='center')
    ax.text(9.3, 0.045, 'resembles the agent', color=C_SAGE, fontsize=12,
            va='bottom', ha='left', fontweight='bold')
    ax.text(9.3, -0.055, 'anti-correlated with the agent', color=C_RED,
            fontsize=12, va='top', ha='left', fontweight='bold')
    ax.set_xscale('log')
    ax.set_xticks([10, 50, 100, 500])
    ax.set_xticklabels(['10', '50', '100', '500'])
    ax.set_xlim(8.6, 1250)
    ax.set_ylim(-0.9, 1.0)
    # 단위는 µM 이다 (2026-11-02 원논문 실측). ppm 은 철회.
    ax.set_xlabel(r'concentration ($\mu$M)')
    ax.set_ylabel('cosine, DMMP response pattern vs agent')
    fig.tight_layout()
    fig.savefig('outputs/Fig_talk_dmmp.png', dpi=200,
                bbox_inches='tight', facecolor='white')
    print('[done] outputs/Fig_talk_dmmp.png')


def fig_tox():
    pan = pd.read_csv('outputs/toxicity_predicted.csv')
    ho = pd.read_csv('outputs/qsar_heldout.csv')
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.2, 5.4))

    # (a) AEGL-3 실측: 작용제 대 모사체 대 간섭물
    d = pan[pan.AEGL3_10min_mgm3.notna()].copy()
    d = d.sort_values('AEGL3_10min_mgm3')
    cmap = {'real_agent': C_RED, 'simulant': C_DUSTY,
            'interferent': C_AMBER}
    y = np.arange(len(d))
    a1.barh(y, d.AEGL3_10min_mgm3, color=[cmap[c] for c in d.compound_class],
            edgecolor='white', lw=0.6)
    a1.set_yticks(y)
    a1.set_yticklabels([n.split(' (')[0] for n in d.compound_name],
                       fontsize=10.5)
    a1.set_xscale('log')
    a1.set_xlabel('EPA AEGL-3, 10 min (mg m$^{-3}$)\n'
                  'lower = lethal at lower airborne concentration')
    missing = ['DMMP', 'DPGME']
    a1.text(0.985, 0.012,
            'no regulator-derived value exists for\n'
            'DMMP or DPGME, two of the three incumbents',
            transform=a1.transAxes, ha='right', va='bottom', fontsize=11.5,
            color=C_RED, fontweight='bold')
    h = [plt.Rectangle((0, 0), 1, 1, color=C_RED, label='authentic agent'),
         plt.Rectangle((0, 0), 1, 1, color=C_DUSTY,
                       label='incumbent simulant'),
         plt.Rectangle((0, 0), 1, 1, color=C_AMBER,
                       label='interferent / solvent')]
    a1.legend(handles=h, loc='lower right', frameon=False, fontsize=11,
              bbox_to_anchor=(1.0, 0.145), handlelength=1.1)
    a1.set_title('(a)  the measured axis', loc='left', fontweight='bold')

    # (b) 모델의 일방 실패
    ok = ho[ho.in_domain]
    no = ho[~ho.in_domain]
    a2.axhspan(0, 4.0, color=C_RED, alpha=0.06, lw=0)
    a2.scatter(ok.ad_distance, ok.error_log10, s=95, color=C_RED,
               edgecolor='white', lw=0.9, zorder=3,
               label='inside the applicability domain (14)')
    a2.scatter(no.ad_distance, no.error_log10, s=95, marker='^',
               facecolor='white', edgecolor=C_BROWN, lw=1.6, zorder=3,
               label='outside the domain (6)')
    a2.axhline(0, color='0.3', lw=1.3)
    a2.set_ylim(-0.4, 4.0)
    a2.set_xlabel('applicability-domain distance')
    a2.set_ylabel('predicted $-$ measured  ($\\log_{10}$)')
    a2.text(0.03, 0.95, 'all 20 agents predicted SAFER than measured',
            transform=a2.transAxes, ha='left', va='top', fontsize=12.5,
            color=C_RED, fontweight='bold')
    # 라벨은 군집 반대쪽으로 민다. 겹치면 어느 점의 이름인지 안 보인다.
    tags = {'Agent GB (Sarin)': ('GB', (-6, 15)),
            'Agent VX': ('VX', (13, -3)),
            'Sulfur Mustard': ('HD', (-8, -17))}
    for nm, (short, off) in tags.items():
        r = ho[ho.aegl_name == nm]
        if len(r):
            a2.annotate(short, (float(r.ad_distance.iloc[0]),
                                float(r.error_log10.iloc[0])),
                        xytext=off, textcoords='offset points',
                        fontsize=12, color=C_BROWN, fontweight='bold',
                        arrowprops=dict(arrowstyle='-', lw=0.8,
                                        color=C_BROWN, alpha=0.6))
    a2.legend(loc='lower right', frameon=False, fontsize=11)
    a2.set_title('(b)  and a model cannot fill the gap',
                 loc='left', fontweight='bold')

    fig.tight_layout()
    fig.savefig('outputs/Fig_talk_tox.png', dpi=200,
                bbox_inches='tight', facecolor='white')
    print('[done] outputs/Fig_talk_tox.png')


def main():
    fig_dmmp()
    fig_tox()
    return 0


if __name__ == '__main__':
    sys.exit(main())
