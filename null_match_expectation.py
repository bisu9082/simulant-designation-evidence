#!/usr/bin/env python3
"""지정 집합 크기별 유효매칭 수의 귀무기댓값 (2026-10-23 신설).

왜 필요한가. 본문과 초록은 '기존 3종 0매칭 대 문헌포함 8종 5매칭'을 나란히
놓는다. 그런데 유효매칭 기준의 둘째 조건(비지정 전체보다 가깝다)은 지정
집합이 클수록 통과하기 쉬우므로, 원시 차이 5 안에는 군크기에서 오는 몫이
섞여 있다. 2.4절이 짝지은 검정을 버리고 전수 순열로 간 이유와 같은 교란이
이 비교에는 남아 있었다 (2026-10-23 외부 검토 지적).

그래서 각 k 에 대해 비작용제 중 k 종을 지정하는 모든(또는 표집한) 경우의
유효매칭 수 분포를 내고, 관측값이 그 기댓값에서 얼마나 벗어나는지 보고한다.
결론을 바꾸지는 않지만, 비교가 군크기로 설명되지 않는다는 것을 수치로 건다.

출력: outputs/null_match_expectation.json
"""
import itertools
import json

import numpy as np
import pandas as pd

from match_by_panel import merged_space, zdist

PANEL = 'descriptor_cache_xtb_35.csv'
THRESHOLD = 2.0
N_SAMPLE = 20000          # k 가 커서 전수가 불가능할 때의 표집 수
SEED = 0xC0DE


def match_count(D, med, agents, designated, non_designated, thr=THRESHOLD):
    """유효매칭: 정규화거리 <= thr AND 비지정 전체보다 가깝다."""
    n = 0
    for a in agents:
        d_sim = min(D[a, s] for s in designated)
        if d_sim / med <= thr and d_sim < min(D[a, x] for x in non_designated):
            n += 1
    return n


def main():
    df = pd.read_csv(PANEL)
    D, dim = zdist(df, merged_space(df))
    np.fill_diagonal(D, np.inf)
    med = float(np.median(D.min(1)))
    kl = df.compound_class.tolist()
    agents = [i for i, c in enumerate(kl) if c == 'real_agent']
    non = [i for i in range(len(df)) if i not in agents]

    SETS = {'incumbent3': [i for i, c in enumerate(kl) if c == 'simulant'],
            'incumbent_plus_lit8': [i for i, c in enumerate(kl)
                                    if c in ('simulant', 'literature_simulant')]}
    rng = np.random.default_rng(SEED)
    out = {'panel': PANEL, 'dim': dim, 'nn_median': med, 'threshold': THRESHOLD,
           'n_agents': len(agents), 'n_non_agents': len(non), 'sets': {}}

    for lab, desig in SETS.items():
        k = len(desig)
        obs = match_count(D, med, agents, desig,
                          [i for i in non if i not in desig])
        total = 1
        for j in range(k):
            total = total * (len(non) - j) // (j + 1)
        exhaustive = total <= 50000
        combos = (list(itertools.combinations(non, k)) if exhaustive else
                  [tuple(rng.choice(non, k, replace=False))
                   for _ in range(N_SAMPLE)])
        counts = np.array([match_count(D, med, agents, list(c),
                                       [i for i in non if i not in c])
                           for c in combos])
        p = float((counts >= obs).mean())
        out['sets'][lab] = {
            'k': k, 'observed': int(obs),
            'null_mean': round(float(counts.mean()), 3),
            'null_sd': round(float(counts.std(ddof=1)), 3),
            'deviation': round(float(obs - counts.mean()), 3),
            'p_at_least_observed': round(p, 4),
            'enumerated': bool(exhaustive), 'n_designations': len(combos),
            'ceiling_11k_over_n': round(len(agents) * k / len(non), 3)}
        print(f'  {lab:22s} k={k} 관측 {obs} | 귀무평균 '
              f'{counts.mean():.3f} | 편차 {obs - counts.mean():+.3f} | p={p:.3f}')

    a, b = out['sets']['incumbent3'], out['sets']['incumbent_plus_lit8']
    out['raw_gap'] = b['observed'] - a['observed']
    out['null_gap'] = round(b['null_mean'] - a['null_mean'], 3)
    out['deviation_gap'] = round(b['deviation'] - a['deviation'], 3)
    print(f"  원시차 {out['raw_gap']} = 귀무기댓값차 {out['null_gap']} "
          f"+ 편차차 {out['deviation_gap']}")

    with open('outputs/null_match_expectation.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('[null-match] outputs/null_match_expectation.json')


if __name__ == '__main__':
    main()
