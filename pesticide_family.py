#!/usr/bin/env python3
"""네 번째 구성의 다중검정 체계 5종 (2026-10-20 신설).

왜 필요한가. 표 2 는 다섯 체계를 세 구성에만 적용했는데, 보정을 통과하는
유일한 구성은 네 번째다. 유의한 쪽만 체계 비교에서 빠져 있으면 '왜 그것만
안 돌렸나' 가 첫 질문이 된다. 같은 정의로 네 번째 열을 채운다.

정의는 stat_supplement.py(BH, Li-Ji meff)와 stat_supplement2.py(Westfall-Young
단일단계 minP)에서 그대로 가져온다. 새 정의를 쓰면 열끼리 비교가 깨진다.

메모리. C(24,10)=1,961,256 순열 x 18 검정이다. float32 로 담아 약 141 MB.
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


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def meff_li_ji(R):
    lam = np.clip(np.linalg.eigvalsh(R), 0, None)
    return float(sum((l >= 1) + (l - np.floor(l)) for l in lam))


def spaces_of(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    merged = [c for c in num if c != OMEGA_DERIVED]
    keep, _ = pick(df, merged, PROTECTED)
    elec3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]
    return {'structural': [c for c in keep if c not in elec3],
            'electronic': [c for c in keep if c in elec3],
            'merged': keep}


def null_matrix(D, kl, sim_idx, label):
    """(n_perm, 6) 귀무행렬과 실제 지정의 행 번호."""
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in sim_idx))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    C = np.asarray(ALL, dtype=np.int32)
    Ds = D[np.ix_(i_real, i_pool)]
    npool = len(i_pool)
    N = np.empty((len(ALL), 6), dtype=np.float32)
    t0 = time.time()
    for a0 in range(0, len(C), CHUNK):
        cb = C[a0:a0 + CHUNK]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        N[a0:a0 + b, 0] = (ds - di).mean(0)
        N[a0:a0 + b, 1] = np.median(ds - di, axis=0)
        N[a0:a0 + b, 2] = np.log(ds / di).mean(0)
        N[a0:a0 + b, 3] = ds.mean(0)
        N[a0:a0 + b, 4] = -(ds < di).sum(0)
        N[a0:a0 + b, 5] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
        if (a0 // CHUNK) % 40 == 0 or a0 + CHUNK >= len(C):
            done = min(a0 + CHUNK, len(C))
            el = time.time() - t0
            sys.stdout.write(f'\r  [{label}] {done:,}/{len(ALL):,} '
                             f'({100*done/len(ALL):5.1f}%)  {el:5.1f}s')
            sys.stdout.flush()
    print()
    return N, ALL.index(true_combo)


def main():
    df = pd.read_csv(PATH)
    kl = df.compound_class.tolist()
    names_c = df.compound_name.tolist()
    pest = [i for i, n in enumerate(names_c) if n in PESTICIDES]
    assert len(pest) == 2, pest
    sim_idx = [i for i, k in enumerate(kl)
               if k in ('simulant', 'literature_simulant')] + pest

    cols_null, obs_p, names, ti = [], [], [], None
    for sp, cols in spaces_of(df).items():
        X = df[cols].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        N, ti = null_matrix(D, kl, sim_idx, sp)
        for j, st in enumerate(STATS):
            col = N[:, j]
            good = np.isfinite(col)
            obs_p.append(float((col[good] <= col[ti]).sum() / good.sum()))
            cols_null.append(col)
            names.append(f'{sp}:{st}')
        del N
    P = np.array(obs_p)
    M = np.vstack(cols_null).T          # (n_perm, 18) float32
    del cols_null
    nperm = M.shape[0]
    print(f'\n귀무행렬 {M.shape}, {M.nbytes/1e6:.0f} MB')

    # (a) 18검정 BH / (b) 공간별 6검정 BH / (c) 통계량 사전지정 3검정
    q_all = bh(P)
    q_by_space = np.concatenate([bh(P[i:i + 6]) for i in (0, 6, 12)])
    idx_md = [names.index(f'{s}:mean_diff')
              for s in ('structural', 'electronic', 'merged')]
    q_pre = bh(P[idx_md])

    # (d) Westfall-Young 단일단계 minP — 18검정의 의존을 순열로 정확히 처리
    PM = np.empty_like(M)
    for j in range(M.shape[1]):
        col = M[:, j]
        good = np.isfinite(col)
        r = np.empty(nperm, dtype=np.float32)
        r[~good] = 1.0
        v = col[good]
        order = np.argsort(v, kind='mergesort')
        sv = v[order]
        ranks = np.empty(len(v))
        ranks[order] = np.searchsorted(sv, sv, side='right')
        r[good] = ranks / good.sum()
        PM[:, j] = r
        sys.stdout.write(f'\r  [WY] 열 {j+1}/18')
        sys.stdout.flush()
    print()
    minp_perm = PM.min(1)
    obs_minp = float(PM[ti].min())
    wy = float((minp_perm <= obs_minp).sum() / nperm)
    del PM

    # (e) Li-Ji 실질 검정 수
    R = np.nan_to_num(np.corrcoef(M, rowvar=False), nan=0.0)
    np.fill_diagonal(R, 1.0)
    me = meff_li_ji(R)

    res = dict(
        label='panel35_incumbent3_plus_lit5_plus_pesticides2',
        n_tests=int(len(P)), n_perm=int(nperm),
        min_p=float(P.min()),
        bh_all_min_q=float(q_all.min()),
        bh_all_n_sig=int((q_all < 0.05).sum()),
        bh_by_space_min_q=float(q_by_space.min()),
        bh_by_space_n_sig=int((q_by_space < 0.05).sum()),
        westfall_young_adjusted_p=wy,
        observed_min_p_across_tests=obs_minp,
        meff=me,
        meff_adjusted_min_p=float(min(1.0, P.min() * me)),
        prespecified_min_p=float(P[idx_md].min()),
        prespecified_bh_min_q=float(q_pre.min()),
        per_test_p=dict(zip(names, map(float, P))))

    print(f'\n### 네 번째 구성 (모사체 10종)')
    print(f'  최소 명목 p (18검정)        {res["min_p"]:.5f}')
    print(f'  (a) BH over 18              {res["bh_all_min_q"]:.5f}  '
          f'유의 {res["bh_all_n_sig"]}건')
    print(f'  (b) BH over 6, per space    {res["bh_by_space_min_q"]:.5f}  '
          f'유의 {res["bh_by_space_n_sig"]}건')
    print(f'  (c) Westfall-Young min-P    {res["westfall_young_adjusted_p"]:.5f}')
    print(f'  (d) m_eff (Li-Ji)           {res["meff"]:.2f} / 18  '
          f'-> Bonferroni {res["meff_adjusted_min_p"]:.5f}')
    print(f'  (e) 통계량 사전지정 (3검정) {res["prespecified_bh_min_q"]:.5f}')
    with open('outputs/pesticide_family.json', 'w') as fh:
        json.dump(res, fh, indent=1)
    print('[pesticide_family] outputs/pesticide_family.json 기록')


if __name__ == '__main__':
    main()
