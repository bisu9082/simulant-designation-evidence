#!/usr/bin/env python3
"""outputs/threshold_sweep.json 생성기 (2026-11-08 신설).

왜 필요한가. 본문 4.3 의 임계값 스윕 문단이 공간별·임계값별 매칭 집합을
서술하는데, 그 수치를 담은 산출물이 없었다. 에디터 게이트 1회차에서
'어느 예치물에도 없다' 로 걸린 자리다. 문장을 고치려면 먼저 산출해
예치해야 한다.

정의는 새로 만들지 않는다. 거리는 match_by_panel.zdist (분산0 제거 ->
|r|>0.95 상관필터 -> 표준화 -> 유클리드), 정규화 분모는 그 공간의 패널
전체 최근접거리 중앙값, 둘째 조건은 넓은 판독(작용제도 지정모사체도
아닌 전부)이다. 공간 정의는 rollback_c 와 같다.
"""
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick
from match_by_panel import zdist

PANEL = 'descriptor_cache_xtb.csv'
SIM_CLS = ['simulant']
THRESHOLDS = [1.5, 1.75, 2.0, 2.5, 3.0]
CODE = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
        'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD',
        'lewisite 1 (L)': 'L', 'hydrogen cyanide (AC)': 'AC',
        'cyanogen chloride (CK)': 'CK', 'phosgene (CG)': 'CG',
        'chloropicrin (PS)': 'PS', 'VX': 'VX'}


def main():
    df = pd.read_csv(PANEL)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    elec3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]
    SPACES = {
        'structural': [c for c in num if c not in ELECTRONIC_28],
        'electronic3': [c for c in num if c in elec3],
        'merged': [c for c in num if c != OMEGA_DERIVED],
    }
    kl = df.compound_class.tolist()
    nm = df.compound_name.tolist()
    i_agent = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in SIM_CLS]
    i_other = [i for i, k in enumerate(kl)
               if k != 'real_agent' and k not in SIM_CLS]

    out = {'generated_by': 'threshold_sweep.py', 'panel': PANEL,
           'reading': 'wide', 'thresholds': THRESHOLDS, 'spaces': {}}
    for sp, cols in SPACES.items():
        D, nk = zdist(df, cols)
        Dd = D.copy()
        np.fill_diagonal(Dd, np.inf)
        nn_med = float(np.median(Dd.min(1)))
        rows = []
        for ia in i_agent:
            ds = [(D[ia, j], nm[j]) for j in i_sim]
            d_sim, who = min(ds)
            d_other = min(D[ia, j] for j in i_other)
            rows.append({
                'agent': nm[ia], 'code': CODE[nm[ia]],
                'nearest_simulant': who,
                'd_norm': round(d_sim / nn_med, 4),
                'beats_comparison_set': bool(d_sim < d_other),
                'blocker': (None if d_sim < d_other else
                            min(((D[ia, j], nm[j]) for j in i_other))[1]),
            })
        sets = {}
        for t in THRESHOLDS:
            sel = sorted(r['code'] for r in rows
                         if r['beats_comparison_set'] and r['d_norm'] <= t)
            sets[f'{t}'] = sel
        sets['unbounded'] = sorted(r['code'] for r in rows
                                   if r['beats_comparison_set'])
        out['spaces'][sp] = {'n_descriptors': nk,
                             'nn_median': round(nn_med, 4),
                             'agents': rows, 'matched_sets': sets,
                             'counts': {k: len(v) for k, v in sets.items()}}
    # 단조성 점검. 임계값이 오르면 매칭 집합은 줄 수 없다.
    for sp, b in out['spaces'].items():
        prev = 0
        for t in [f'{x}' for x in THRESHOLDS] + ['unbounded']:
            c = b['counts'][t]
            if c < prev:
                raise SystemExit(f'단조성 위반 {sp} {t}: {c} < {prev}')
            prev = c
    json.dump(out, open('outputs/threshold_sweep.json', 'w'),
              indent=1, ensure_ascii=False)
    for sp, b in out['spaces'].items():
        print(f'{sp:12} dims={b["n_descriptors"]:3} nn_med={b["nn_median"]}')
        for t, v in b['matched_sets'].items():
            print(f'   {t:>9}: {len(v)}  {{{", ".join(v)}}}')


if __name__ == '__main__':
    main()
