#!/usr/bin/env python3
"""블록 균등 가중 거리 민감도 (2026-10-17 신설).

제기된 쟁점. 병합공간은 구조 서술자 38개와 전자 서술자 3개를 한 유클리드
거리에 넣는다. 변수별 표준화는 끝났지만 블록당 기여도는 38 대 3 이고,
따라서 병합공간의 거리는 사실상 구조공간의 거리다. 본문도 이 사실을
적어 두었으나 수치로 확인한 적은 없다.

검정 방식. 블록당 기대 기여도를 같게 만든 거리를 쓴다.

    d_bal^2 = (1/n_s) * sum_i (dz_s,i)^2 + (1/n_e) * sum_j (dz_e,j)^2

그리고 본문과 같은 전수 순열검정을 다시 돌린다. 결론이 유지되면 병합공간
결과가 가중에 의존하지 않는다는 뜻이고, 뒤집히면 병합공간에 기초한 진술을
제한해야 한다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick
from perm_statistics import STATS

CASES = [
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
]
ELEC3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]


def blocks(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    merged = [c for c in num if c != OMEGA_DERIVED]
    keep, _ = pick(df, merged, PROTECTED)
    elec = [c for c in keep if c in ELEC3]
    struct = [c for c in keep if c not in ELEC3]
    return struct, elec


def dist(df, struct, elec, balanced):
    def Z(cols):
        X = df[cols].to_numpy(float)
        return (X - X.mean(0)) / X.std(0, ddof=1)
    Zs, Ze = Z(struct), Z(elec)
    ds2 = ((Zs[:, None, :] - Zs[None, :, :]) ** 2).sum(-1)
    de2 = ((Ze[:, None, :] - Ze[None, :, :]) ** 2).sum(-1)
    if balanced:
        return np.sqrt(ds2 / len(struct) + de2 / len(elec))
    return np.sqrt(ds2 + de2)


def exact_perm(D, kl, sim_classes):
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    null = {s: np.empty(len(ALL)) for s in STATS}
    npool = len(i_pool)
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
        null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
        null['mean_logratio'][a0:a0 + b] = np.log(ds / di).mean(0)
        null['mean_dsim'][a0:a0 + b] = ds.mean(0)
        null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
        null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    ti = ALL.index(true_combo)
    return ({s: float((null[s] <= null[s][ti]).sum() / len(ALL)) for s in STATS},
            len(ALL), nn_med,
            int(-null['n_valid'][ti]))


def main():
    out = {}
    for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        struct, elec = blocks(df)
        row = {'n_struct': len(struct), 'n_elec': len(elec)}
        print(f'\n### {label}  (구조 {len(struct)} / 전자 {len(elec)})')
        for tag, bal in [('equal_per_descriptor', False), ('block_balanced', True)]:
            D = dist(df, struct, elec, bal)
            ps, nperm, nnm, nvalid = exact_perm(D, df.compound_class.tolist(),
                                                sim_cls)
            row[tag] = dict(p=ps, n_perm=nperm, nn_median=nnm, n_valid=nvalid)
            s = ', '.join(f'{k}={ps[k]:.3f}' for k in STATS)
            print(f'  {tag:<22} 최소p {min(ps.values()):.4f}  유효매칭 {nvalid}/11')
            print(f'    {s}')
        out[label] = row
    with open('outputs/blockbalance.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print('\n[blockbalance] outputs/blockbalance.json 기록')


if __name__ == '__main__':
    main()
