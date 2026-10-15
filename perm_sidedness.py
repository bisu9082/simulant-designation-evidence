#!/usr/bin/env python3
"""단측·양측 p 를 생성한다 (2026-10-13 신설).

outputs/permutation_sidedness.json 은 verify_consistency.py 가 본문 대조에 쓰는
사실 파일인데, 지금까지 이 파일을 만드는 스크립트가 없었다. 손으로 한 번 써 놓고
아무도 갱신하지 않았으므로, omega 를 거리공간에서 빼는 변경이 본문에는 반영되고
검증기의 사실값에는 반영되지 않았다. 검증기가 통과하면서 낡은 값을 들고 있는
상태였다. 생성기를 둬서 그 경로를 막는다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM,
                           OMEGA_DERIVED, PROTECTED, pick)

PATH = 'descriptor_cache_xtb.csv'


def main():
    df = pd.read_csv(PATH)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    SP = {'A_rdkit': [c for c in num if c not in ELECTRONIC_28],
          'B_electronic': [c for c in num if c in ELECTRONIC_DIST],
          'C_merged': [c for c in num if c != OMEGA_DERIVED]}
    kl = df.compound_class.tolist()
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k == 'simulant']
    pos = {g: j for j, g in enumerate(i_pool)}
    true = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true)))

    out = {}
    for sp, cols in SP.items():
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        D = np.sqrt(((Z[:, None] - Z[None]) ** 2).sum(-1))
        np.fill_diagonal(D, np.inf)
        Ds = D[np.ix_(i_real, i_pool)]
        C = np.asarray(ALL, dtype=np.int32)
        npool = len(i_pool)
        null = np.empty(len(ALL))
        for a0 in range(0, len(C), 8000):
            cb = C[a0:a0 + 8000]
            b = len(cb)
            m = np.zeros((b, npool), bool)
            m[np.arange(b)[:, None], cb] = True
            A3 = Ds[:, None, :]
            ds = np.where(m[None], A3, np.inf).min(2)
            di = np.where(~m[None], A3, np.inf).min(2)
            null[a0:a0 + b] = (ds - di).mean(0)
        obs = null[ALL.index(true)]
        lo = float((null <= obs).mean())
        hi = float((null >= obs).mean())
        out[sp] = dict(
            n_kept=len(keep),
            n_perm=len(ALL),
            one_sided=lo,
            two_sided_2min=float(min(1.0, 2 * min(lo, hi))),
            two_sided_abs=float((np.abs(null - null.mean())
                                 >= abs(obs - null.mean())).mean()))
        print(f'  {sp:<13} {len(keep):>2}차원  단측 {lo:.4f}  '
              f'양측(2x) {out[sp]["two_sided_2min"]:.4f}  '
              f'양측(편차) {out[sp]["two_sided_abs"]:.4f}')
    with open('outputs/permutation_sidedness.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print('[perm_sidedness] outputs/permutation_sidedness.json 기록')


if __name__ == '__main__':
    main()
