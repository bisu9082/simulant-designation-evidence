#!/usr/bin/env python3
"""'라벨 둘' 주장의 실제 구성 (2026-11-12 신설).

왜 필요한가. main.tex 가 "Classing them as simulants instead, which moves
two labels and changes nothing else, gives ... a Benjamini-Hochberg value of
0.0008" 라고 적고, 기준을 "nothing at q=0.630" 으로 명시한다. q=0.630 은
35종 패널에 기존 3종만 지정한 구성이다.

그런데 0.0008 을 내는 구성(pesticide_panel.py)의 지정집합은 10종이다 --
기존 3 + 문헌 5 + 농약 2. 즉 기준에서 **일곱 라벨**이 움직인다. 농약 둘만
옮긴 구성은 아무도 돌린 적이 없다.

이 스크립트가 그것을 돌린다. 35종 패널, 지정 = 기존 3 + 농약 2 = 5종,
C(24,5)=42,504 전수열거. 정의는 pesticide_panel.py 를 그대로 따른다 --
다른 정의를 쓰면 열끼리 비교가 깨진다.
"""
import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick
from perm_statistics import STATS
from pesticide_panel import spaces_of, dmat, bh_min, PATH, PESTICIDES, CHUNK


def exact_perm(D, kl, sim_idx, label):
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in sim_idx))
    k = len(true_combo)
    ALL = list(itertools.combinations(range(len(i_pool)), k))
    C = np.asarray(ALL, dtype=np.int32)
    Ds = D[np.ix_(i_real, i_pool)]
    npool = len(i_pool)
    null = {s: np.empty(len(ALL)) for s in STATS}
    t0 = time.time()
    for a0 in range(0, len(C), CHUNK):
        cb = C[a0:a0 + CHUNK]
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
    el = time.time() - t0
    print(f'  [{label}] {len(ALL):,} 지정 전수, {el:.1f}s')
    ti = ALL.index(true_combo)
    p = {s: float((null[s] <= null[s][ti]).sum() / len(ALL)) for s in STATS}
    return p, len(ALL), nn_med, int(-null['n_valid'][ti])


def main():
    df = pd.read_csv(PATH)
    kl = df.compound_class.tolist()
    names = df.compound_name.tolist()
    sp = spaces_of(df)

    pest = [i for i, n in enumerate(names) if n in PESTICIDES]
    if len(pest) != len(PESTICIDES):
        sys.exit(f'농약을 찾지 못했다: {[names[i] for i in pest]}')
    # 기존 3종 + 농약 2종. 문헌 5종은 넣지 않는다 -- 그것이 요점이다.
    sim_idx = [i for i, k in enumerate(kl) if k == 'simulant'] + pest
    print(f'패널 {len(df)}종, 작용제 {kl.count("real_agent")}종, '
          f'지정 {len(sim_idx)}종 '
          f'({", ".join(names[i] for i in sorted(sim_idx))})')

    out = {}
    for lab, cols in sp.items():
        p, nperm, nnm, nvalid = exact_perm(dmat(df, cols), kl, sim_idx, lab)
        out[lab] = dict(p=p, n_perm=nperm, nn_median=nnm, n_valid=nvalid,
                        n_dims=len(cols))
        print(f'    {lab:<12} 유효매칭 {nvalid}/11  최소p {min(p.values()):.4f}')

    allp = [out[s]['p'][k] for s in sp for k in STATS]
    res = dict(panel='panel35_incumbent3_plus_pesticides2',
               note='main.tex 가 "moves two labels" 라 부르는 구성의 실제 값. '
                    '0.0008 을 내는 구성은 문헌 5종이 함께 지정된 10종이다.',
               n_simulants=len(sim_idx),
               simulant_names=sorted(names[i] for i in sim_idx),
               spaces=out,
               min_p=float(min(allp)),
               bh_min_q=float(bh_min(allp)))
    print(f'\n  18검정 최소 p {res["min_p"]:.6f}  BH 최소 q {res["bh_min_q"]:.4f}')
    with open('outputs/twolabel_panel.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('[twolabel_panel] outputs/twolabel_panel.json 기록')


if __name__ == '__main__':
    main()
