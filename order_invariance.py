#!/usr/bin/env python3
"""상관필터의 열 순서 의존성 측정 (2026-11-02 신설, T2-1).

왜 필요한가. GATE 8 7차에서 리뷰어 R4 가 pick() 이 입력 열 순서의 함수임을
지적했다. pick() 은 |r|>0.95 쌍을 (i, j) 이중 루프로 훑으며 뒤쪽 열을
버리므로, 같은 상관군에서 어느 열이 살아남는지가 CSV 의 열 배치로 결정된다.
헤드라인인 '유효매칭 4건'이 그 배치의 함수라면 그것은 결과가 아니다.

이 스크립트는 주장하지 않고 센다. 입력 열 순서를 N 회 무작위 치환하고,
치환마다 (차원 수, 유지 열 집합, 작용제별 매칭, 유효매칭 수) 를 기록해
분포를 낸다. 결정적 필터로 교체한 뒤 같은 스크립트를 다시 돌려 분포가
한 점으로 모이는지 확인한다.

사용:
  python3 order_invariance.py --n 200
  python3 order_invariance.py --n 200 --out outputs/order_invariance.json
"""
import argparse
import json
import random
from collections import Counter

import numpy as np
import pandas as pd

from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, _pick_legacy, pick

RULES = {'deterministic': pick, 'legacy': _pick_legacy}

CASES = [
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
]
THRESHOLD = 2.0


def merged_space(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return [c for c in num if c != OMEGA_DERIVED]


def one_order(df, cols, sim_cls, rule=pick):
    """match_by_panel.one_panel 과 같은 거리·매칭 규칙. 열 순서만 바꿔 호출한다."""
    keep, _ = rule(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))

    kl = df.compound_class.tolist()
    name = df.compound_name.tolist()
    i_agent = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    i_other = [i for i, k in enumerate(kl)
               if k != 'real_agent' and k not in sim_cls]

    matched, dn = [], {}
    for i in i_agent:
        js = min(i_sim, key=lambda j: D[i, j])
        d_norm = D[i, js] / nn_med
        ok = bool(d_norm <= THRESHOLD and all(D[i, js] < D[i, j]
                                              for j in i_other))
        dn[name[i]] = round(float(d_norm), 4)
        if ok:
            matched.append(name[i])
    return sorted(keep), sorted(matched), dn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--seed', type=int, default=20261102)
    ap.add_argument('--out', default='outputs/order_invariance.json')
    ap.add_argument('--rule', choices=sorted(RULES) + ['both'],
                    default='both')
    a = ap.parse_args()

    rng = random.Random(a.seed)
    names = sorted(RULES) if a.rule == 'both' else [a.rule]
    report = {'n_permutations': a.n, 'seed': a.seed, 'threshold': THRESHOLD,
              'rules': names, 'cases': {}}

    for rname in names:
      rule = RULES[rname]
      for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        base = merged_space(df)
        dims, matchsets, nmatch = Counter(), Counter(), Counter()
        dnorm_runs = {}
        keepsets = Counter()

        orders = [list(base)]                      # 1회는 원래 순서
        for _ in range(a.n - 1):
            o = list(base)
            rng.shuffle(o)
            orders.append(o)

        for k, o in enumerate(orders):
            keep, matched, dn = one_order(df, o, sim_cls, rule)
            dims[len(keep)] += 1
            keepsets[tuple(keep)] += 1
            matchsets[tuple(matched)] += 1
            nmatch[len(matched)] += 1
            for ag, v in dn.items():
                dnorm_runs.setdefault(ag, []).append(v)

        agent_spread = {ag: {'min': min(v), 'max': max(v),
                             'n_distinct': len(set(v))}
                        for ag, v in dnorm_runs.items()}
        report['cases'][f'{rname}/{label}'] = {
            'n_input_columns': len(base),
            'dims': dict(sorted(dims.items())),
            'n_distinct_kept_sets': len(keepsets),
            'n_match': dict(sorted(nmatch.items())),
            'match_sets': {('|'.join(k) if k else '(none)'): v
                           for k, v in matchsets.most_common()},
            'invariant': len(keepsets) == 1 and len(nmatch) == 1,
            'agent_dnorm_spread': agent_spread,
        }
        c = report['cases'][f'{rname}/{label}']
        print(f"[{rname}] [{label}] 입력열 {len(base)}  "
              f"유지집합 {c['n_distinct_kept_sets']}종  "
              f"차원 {dict(sorted(dims.items()))}  "
              f"매칭수 {dict(sorted(nmatch.items()))}  "
              f"불변 {c['invariant']}")

    json.dump(report, open(a.out, 'w'), indent=1, ensure_ascii=False)
    print(f'-> {a.out}')


if __name__ == '__main__':
    main()
