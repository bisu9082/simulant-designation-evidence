#!/usr/bin/env python3
"""네 번째 패널 — 농약 2종을 모사체로 분류했을 때 (2026-10-20 신설).

왜 필요한가. 본문 §2.1 은 dichlorvos 와 malathion 을 간섭물로 분류하면서
그 결정을 '쟁점의 여지가 있다'고 스스로 적었고, Limitations 은 '그 둘을
모사체로 분류한 독자는 우리가 보고하지 않는 수치를 얻을 것'이라고 적었다.
외부 검토 두 건이 같은 자리를 짚었다. 보고하지 않은 수치를 실제로 계산해
그 문장을 손짓이 아니라 데이터로 바꾼다.

구성. 35종 패널 그대로, class 라벨만 바꾼다.
  모사체 10종 = 현행 3 + 문헌 5 + 농약 2
  비작용제 풀 24종, 순열 C(24,10) = 1,961,256 — 전수 열거 가능.

통계량·거리·매칭 규칙은 본문과 동일하다. 새 규칙을 쓰면 비교가 깨진다.
"""
import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick
from perm_statistics import STATS

PATH = 'descriptor_cache_xtb_35.csv'
PESTICIDES = ['dichlorvos', 'malathion']
CHUNK = 4000


def spaces_of(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    merged = [c for c in num if c != OMEGA_DERIVED]
    keep, _ = pick(df, merged, PROTECTED)
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
        if (a0 // CHUNK) % 25 == 0 or a0 + CHUNK >= len(C):
            done = min(a0 + CHUNK, len(C))
            el = time.time() - t0
            eta = el / done * (len(C) - done)
            sys.stdout.write(f'\r  [{label}] {done:,}/{len(ALL):,} '
                             f'({100*done/len(ALL):5.1f}%)  경과 {el:5.1f}s '
                             f' 남은 {eta:5.1f}s')
            sys.stdout.flush()
    print()
    ti = ALL.index(true_combo)
    p = {s: float((null[s] <= null[s][ti]).sum() / len(ALL)) for s in STATS}
    return p, len(ALL), nn_med, int(-null['n_valid'][ti])


def bh_min(ps):
    o = sorted(ps)
    m = len(o)
    q = [o[i] * m / (i + 1) for i in range(m)]
    for i in range(m - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    return q[0]


def main():
    df = pd.read_csv(PATH)
    kl = df.compound_class.tolist()
    names = df.compound_name.tolist()
    sp = spaces_of(df)

    # 농약 2종을 모사체로 승격. 본문 라벨은 건드리지 않고 여기서만 바꾼다.
    pest = [i for i, n in enumerate(names) if n in PESTICIDES]
    if len(pest) != len(PESTICIDES):
        found = [names[i] for i in pest]
        sys.exit(f'농약을 찾지 못했다: 기대 {PESTICIDES}, 발견 {found}')
    sim_idx = [i for i, k in enumerate(kl)
               if k in ('simulant', 'literature_simulant')] + pest
    print(f'패널 {len(df)}종, 작용제 {kl.count("real_agent")}종, '
          f'모사체 {len(sim_idx)}종 '
          f'({", ".join(names[i] for i in sorted(sim_idx))})')
    print(f'차원: 구조 {len(sp["structural"])}, 전자 {len(sp["electronic3"])}, '
          f'병합 {len(sp["merged"])}')

    out = {}
    for lab, cols in sp.items():
        p, nperm, nnm, nvalid = exact_perm(dmat(df, cols), kl, sim_idx, lab)
        out[lab] = dict(p=p, n_perm=nperm, nn_median=nnm, n_valid=nvalid,
                        n_dims=len(cols))
        print(f'    {lab:<12} 유효매칭 {nvalid}/11  최소p {min(p.values()):.4f}  '
              + '  '.join(f'{k}={p[k]:.3f}' for k in STATS))

    allp = [out[s]['p'][k] for s in sp for k in STATS]
    res = dict(panel='panel35_incumbent3_plus_lit5_plus_pesticides2',
               n_simulants=len(sim_idx),
               simulant_names=sorted(names[i] for i in sim_idx),
               spaces=out,
               min_p=float(min(allp)),
               bh_min_q=float(bh_min(allp)),
               n_sig_bh=int(sum(1 for x in [bh_min(allp)] if x < 0.05)))
    print(f'\n  18검정 최소 p {res["min_p"]:.4f}  BH 최소 q {res["bh_min_q"]:.3f}')
    with open('outputs/pesticide_panel.json', 'w') as fh:
        json.dump(res, fh, indent=1)
    print('[pesticide_panel] outputs/pesticide_panel.json 기록')


if __name__ == '__main__':
    main()
