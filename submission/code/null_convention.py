#!/usr/bin/env python3
"""귀무의 비교군 구성 규약 — 다섯 번째 분석자 자유도 (2026-11-13 신설).

왜 필요한가. 본문은 사전 고정하지 않은 선택으로 패널·공간·통계량·가중
네 가지를 열거한다. 그런데 순열검정의 귀무를 어떻게 짜느냐가 다섯 번째
선택이고, 그것이 이 논문의 유일한 명목 유의값을 0.05 양쪽으로 넘긴다.
GATE 8 11차에서 편집장 페르소나가 적발했고, 그때까지 이 선택은
perm_statistics.py 의 주석에만 있었다.

세 규약:
  complement  — 한 라벨링에서 지정된 k 종이 모사체이고 비작용제 풀의
                나머지 전부가 비교군이다. 원고가 쓴 규약.
  fixed       — 비교군을 실제 간섭물 집합에 고정한다. 라벨링이 바뀌어도
                비교군이 바뀌지 않는다.
  minus       — 실제 간섭물 집합에서 그 라벨링이 지정한 것을 뺀다.

[C-3] 전수열거이고 예치 캐시만 읽는다. 새 계산 없음.
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, ELECTRONIC_DIST, META_NUM, OMEGA_DERIVED
from spaces import PROTECTED, pick

PANEL = 'descriptor_cache_xtb.csv'


def run(df, cols, conv):
    kl = df.compound_class.tolist()
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k == 'simulant']
    i_int = [i for i, k in enumerate(kl) if k == 'interferent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true)))

    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    Ds = D[np.ix_(i_real, i_pool)]
    int_mask = np.zeros(len(i_pool), bool)
    for g in i_int:
        int_mask[pos[g]] = True

    stat = np.empty(len(ALL))
    for a, combo in enumerate(ALL):
        sim = np.zeros(len(i_pool), bool)
        sim[list(combo)] = True
        if conv == 'complement':
            cmp_ = ~sim
        elif conv == 'fixed':
            cmp_ = int_mask
        elif conv == 'minus':
            cmp_ = int_mask & ~sim
        ds = np.where(sim, Ds, np.inf).min(1)
        di = np.where(cmp_, Ds, np.inf).min(1)
        stat[a] = (ds - di).mean()
    t = ALL.index(true)
    return float((stat <= stat[t]).sum() / len(ALL)), len(ALL)


def spaces_of(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return {'A_rdkit': [c for c in num if c not in ELECTRONIC_28],
            'B_electronic': [c for c in num if c in ELECTRONIC_DIST],
            'C_merged': [c for c in num if c != OMEGA_DERIVED]}


def main():
    df = pd.read_csv(PANEL)
    SP = spaces_of(df)
    out = {'generated_by': 'null_convention.py',
           'panel': 'panel28_incumbent3',
           'statistic': 'mean_diff (the prespecified one)',
           'note': ('귀무의 비교군 구성 규약 세 가지. 원고는 complement 를 '
                    '쓴다. 이 선택은 사전 고정하지 않은 다섯 번째 자유도다.'),
           'conventions': {}}
    for conv in ('complement', 'fixed', 'minus'):
        out['conventions'][conv] = {}
        for sp, cols in SP.items():
            p, n = run(df, cols, conv)
            out['conventions'][conv][sp] = {'p_value': round(p, 6),
                                            'n_perm': n}
            print(f'  {conv:<11} {sp:<12} p = {p:.4f}  (n={n})')
    with open('outputs/null_convention.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print('  -> outputs/null_convention.json')


if __name__ == '__main__':
    main()
