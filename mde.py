#!/usr/bin/env python3
"""최소검출효과 (minimum detectable effect) — 음성 결과의 검정력 (2026-10-17 신설).

왜 필요한가. 본문 헤드라인은 음성이다. 주 패널에서 작용제 11종 전부가 유효
매칭 0 이고 어떤 검정도 p=0.15 밑으로 내려가지 않는다. 그런데 '없다'와
'볼 수 없었다'는 다르다. 이 스크립트는 후자를 배제한다.

방법. 기존 모사체 3종이 실제보다 delta 만큼(그 공간의 전형 최근접거리 단위)
작용제에 가까이 앉아 있었다면 검정이 무엇을 반환했을지 계산한다. 패널의
나머지는 건드리지 않는다. 귀무분포는 수정하지 않은 거리행렬에서 만든 그대로
쓴다. 해석은 조건부다. '이 패널에서, 지정된 3종만 delta 만큼 더 충실했다면
검정은 p<0.05 를 반환했을 것이다.'

통계량은 mean_diff 하나다. 본문이 '사전에 지정했어야 할 단일 통계량'으로
식별한 것이고, 검정력 진술은 하나의 통계량에 대해서만 뜻이 있다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, pick

CASES = [
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
]
ALPHA = 0.05
GRID = np.round(np.arange(0.0, 1.501, 0.01), 3)


def merged_D(path):
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    cols = [c for c in num if c != OMEGA_DERIVED]
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return df, D, len(keep)


def null_mean_diff(D, kl, sim_cls):
    """수정 전 거리행렬에서 만든 mean_diff 귀무분포와 실제 지정의 위치."""
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    out = np.empty(len(ALL))
    npool = len(i_pool)
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        out[a0:a0 + b] = (ds - di).mean(0)
    return out, ALL.index(true_combo), i_real, i_pool, i_sim


def observed_at(D, i_real, i_sim, i_other, nn_med, delta):
    """지정 3종의 작용제까지 거리만 delta*nn_med 줄였을 때의 mean_diff."""
    ds = np.array([min(max(D[a, s] - delta * nn_med, 1e-9) for s in i_sim)
                   for a in i_real])
    di = np.array([min(D[a, o] for o in i_other) for a in i_real])
    return float(np.mean(ds - di)), ds, di


def run(label, path, sim_cls):
    df, D, nk = merged_D(path)
    kl = df.compound_class.tolist()
    name = df.compound_name.tolist()
    nn_med = float(np.median(D.min(1)))
    null, ti, i_real, i_pool, i_sim = null_mean_diff(D, kl, sim_cls)
    i_other = [i for i in i_pool if i not in i_sim]

    p0 = float((null <= null[ti]).sum() / len(null))
    print(f'\n### {label}  ({nk}차원, 귀무 {len(null):,}, 전형최근접 {nn_med:.3f})')
    print(f'  delta=0 에서 observed {null[ti]:+.4f}, p = {p0:.4f}  (본문값과 대조)')

    rows, mde = [], None
    for d in GRID:
        obs, _, _ = observed_at(D, i_real, i_sim, i_other, nn_med, d)
        p = float((null <= obs).sum() / len(null))
        rows.append(dict(delta=float(d), observed=obs, p=p))
        if mde is None and p < ALPHA:
            mde = float(d)
    if mde is None:
        print(f'  delta {GRID[-1]} 까지 p<{ALPHA} 미달 — 이 패널에서는 '
              f'전형거리 {GRID[-1]}배 접근으로도 검출 불가')
    else:
        print(f'  최소검출효과: delta = {mde:.2f} 전형최근접거리')
        print(f'    = 지정 3종이 작용제에 {mde:.2f} 더 가까이 앉아 있었다면 '
              f'p<{ALPHA} 였다')

    # 실제로 존재하는 변위 폭과 비교한다. 이게 해석의 핵심.
    obs0, ds0, di0 = observed_at(D, i_real, i_sim, i_other, nn_med, 0.0)
    gap = (ds0 - di0) / nn_med
    print(f'    실제 격차(정규화): 중위 {np.median(gap):+.3f}, '
          f'최소 {gap.min():+.3f} ({name[i_real[int(np.argmin(gap))]]})')
    return dict(label=label, n_dims=nk, n_perm=len(null), nn_median=nn_med,
                p_at_zero=p0, mde_delta=mde, alpha=ALPHA,
                observed_at_zero=float(null[ti]),
                gap_median=float(np.median(gap)), gap_min=float(gap.min()),
                curve=rows)


def main():
    res = [run(*c) for c in CASES]
    with open('outputs/mde.json', 'w') as fh:
        json.dump({r['label']: r for r in res}, fh, indent=1)
    print('\n[mde] outputs/mde.json 기록')
    return res


if __name__ == '__main__':
    main()
