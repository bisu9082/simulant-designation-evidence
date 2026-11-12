#!/usr/bin/env python3
"""형상 서술자 축약이 유효매칭 수를 바꾸는가 (2026-11-12 신설).

왜 필요한가. outputs/shape_reduction.json 이 순열 p 는 담지만 유효매칭
수는 담지 않는다. 원고의 헤드라인 대비("28종에서 3, 35종에서 0")가 축약
후에도 서는지는 매칭 수로만 답할 수 있다.

축약 규칙은 shape_reduction.json 의 선언을 그대로 쓴다 -- 유도 형상 7개를
빼고 원시 관성량 PMI1/2/3 만 남긴다. 사유는 그 파일의 rationale 과 같다:
일곱이 관성텐서 3-벡터의 비(ratio) 함수이므로 선형 상관필터가 못 잡는다.

매칭 정의는 match_by_panel.py 와 같다: 정규화 거리 <= 2 이고, 지정집합
밖 모든 화합물보다 가까울 것.
"""
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick

SR = json.load(open('outputs/shape_reduction.json'))
DROP = SR['shape_derived_dropped']

CASES = [
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
]


def spaces_of(df, drop_shape):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    cols = [c for c in num if c != OMEGA_DERIVED]
    if drop_shape:
        cols = [c for c in cols if c not in DROP]
    keep, _ = pick(df, cols, PROTECTED)
    elec3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]
    elec = [c for c in keep if c in elec3]
    struct = [c for c in keep if c not in elec3]
    return {'structural': struct, 'electronic3': elec, 'merged': keep}


def dmat(df, cols):
    X = df[cols].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D


def matches(D, kl, names, sim_cls):
    nn = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    i_cmp = [i for i, k in enumerate(kl)
             if k != 'real_agent' and i not in i_sim]
    out = []
    for a in i_real:
        s = min(i_sim, key=lambda j: D[a, j])
        d = D[a, s] / nn
        beats = all(D[a, s] <= D[a, o] for o in i_cmp)
        out.append(dict(agent=names[a], nearest=names[s], d_norm=round(d, 4),
                        within=bool(d <= 2.0), beats=bool(beats),
                        match=bool(d <= 2.0 and beats)))
    return nn, out


def main():
    res = {'generated_by': 'shape_match.py',
           'dropped': DROP,
           'rationale': SR['rationale']}
    for lab, path, cls in CASES:
        df = pd.read_csv(path)
        kl = df.compound_class.tolist()
        names = df.compound_name.tolist()
        res[lab] = {}
        print(f'\n### {lab}')
        for var, drop in (('full', False), ('reduced', True)):
            sp = spaces_of(df, drop)
            blk = {}
            for s in ('structural', 'merged'):
                nn, rows = matches(dmat(df, sp[s]), kl, names, cls)
                n = sum(r['match'] for r in rows)
                blk[s] = dict(n_dims=len(sp[s]), nn_median=round(nn, 4),
                              n_match=n, agents=rows)
                print(f'  {var:<8} {s:<11} {len(sp[s]):>3}차원  '
                      f'nn_med {nn:6.3f}  유효매칭 {n}/11')
            res[lab][var] = blk
    with open('outputs/shape_match.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('\n[shape_match] outputs/shape_match.json 기록')


if __name__ == '__main__':
    main()
