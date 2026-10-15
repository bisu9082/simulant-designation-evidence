#!/usr/bin/env python3
"""κ 병리값(음수 형태지수) 제외 조건에서 C(17,3)=680 전수 치환검정 재실행.

원고가 "all analyses were repeated with them excluded and the conclusions are
unchanged" 라고 쓰면서 헤드라인 통계인 치환 p 를 제시하지 않았다는 지적에 대응.
"""
import json
from itertools import combinations

import numpy as np
import pandas as pd

df = pd.read_csv('descriptor_cache_xtb.csv')
S = json.load(open('outputs/experiment_summary_nopatho.json'))
CON = {c['condition']: c for c in S['conditions']}

klass = df.compound_class.to_numpy()
i_real = np.flatnonzero(klass == 'real_agent')
i_pool = np.flatnonzero(klass != 'real_agent')
i_sim = set(np.flatnonzero(klass == 'simulant').tolist())
pool_sim = tuple(sorted(list(i_pool).index(j) for j in i_sim))


def zmat(cols):
    X = df[cols].to_numpy(float)
    return (X - X.mean(0)) / X.std(0, ddof=1)


def dmat(Z):
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D


ALL = list(combinations(range(len(i_pool)), 3))
assert len(ALL) == 680
out = {}
print(f'{"condition":14s} {"kept":>5s} {"obs":>9s} {"1-sided":>8s} '
      f'{"2-sided":>8s} {"intf closer":>12s}')
for cond in ['A_rdkit', 'B_electronic', 'C_merged']:
    D = dmat(zmat(CON[cond]['kept']))
    Dsub = D[np.ix_(i_real, i_pool)]
    diff = np.empty(680)
    cnt = np.empty(680)
    ti = None
    for t, combo in enumerate(ALL):
        m = np.zeros(len(i_pool), bool)
        m[list(combo)] = True
        ds = Dsub[:, m].min(1)
        di = Dsub[:, ~m].min(1)
        diff[t] = (ds - di).mean()
        cnt[t] = (di < ds).sum()
        if combo == pool_sim:
            ti = t
    obs = diff[ti]
    p1 = float((diff <= obs).mean())
    mu = diff.mean()
    p2 = float((np.abs(diff - mu) >= abs(obs - mu)).mean())
    out[cond] = dict(n_kept=len(CON[cond]['kept']), observed_diff=float(obs),
                     p_one_sided=p1, p_two_sided_abs=p2,
                     observed_count=float(cnt[ti]), null_mean_count=float(cnt.mean()))
    print(f'{cond:14s} {len(CON[cond]["kept"]):5d} {obs:+9.4f} {p1:8.4f} '
          f'{p2:8.4f} {cnt[ti]:6.0f} / 11')
json.dump(out, open('outputs/permutation_nopatho.json', 'w'), indent=1)
print('\n[done] outputs/permutation_nopatho.json')
