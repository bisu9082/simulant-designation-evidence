#!/usr/bin/env python3
"""통계 보완 2차 (2026-10-14). GATE 8 5차 리뷰어 1 R1-4 / 리뷰어 2 R2-4 대응.

R1-4  Figure 2 의 결론은 작은 거리차에 걸려 있다. 주 패널에서 GB 까지의 거리가
      DMMP 1.19, DFP 1.07 로 0.12 차이다. 서술자 선택·표준화가 이 정도는 움직일
      수 있으므로 거리에 불확실성을 붙여야 한다. 병합공간 서술자를 복원추출로
      재표집해 (a) 작용제별 두 거리의 차이의 구간, (b) 유효매칭 0/11 이 유지되는
      복제 비율을 낸다. 후자가 헤드라인의 안정성을 직접 말한다.

R2-4  meff (Li & Ji) 는 원래 검정통계량의 상관에 쓰는 방법인데 우리는 귀무분포
      벡터의 상관에 썼다. 유사하지만 같지 않다. 18검정이 같은 순열 배열을
      공유하므로, 의존을 근사하지 말고 결합 귀무분포에서 직접 다룰 수 있다.
      Westfall-Young 단일단계 minP: 각 순열마다 18검정의 p 를 구해 최솟값을 취하고,
      그 분포에서 관측 최소 p 의 위치를 조정 p 로 읽는다. 의존이 정확히 반영된다.

출력: outputs/stat_supplement2.json
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM, OMEGA_DERIVED,
                    PROTECTED, pick)
from stat_supplement import STATS, null_matrix, spaces_for, zdist

SEED = 0xC0DE
N_BOOT = 1000


# ------------------------------------------------- R2-4 Westfall-Young minP
def westfall_young(path, sim_classes, label):
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    kl = df.compound_class.tolist()
    cols_null, obs_idx, names = [], None, []
    for sp, cols in spaces_for(df, num).items():
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        N, ti = null_matrix(zdist(Z), kl, sim_classes)
        obs_idx = ti
        for j, st in enumerate(STATS):
            cols_null.append(N[:, j])
            names.append(f'{sp}:{st}')
    M = np.vstack(cols_null).T                       # (n_perm, 18)
    nperm = M.shape[0]

    # 각 열에서 '그 값 이하인 순열의 비율' 이 그 순열의 단측 p 다.
    # 작을수록 모사체가 좋다는 방향으로 통일돼 있으므로 오름차순 순위를 쓴다.
    P = np.empty_like(M)
    for j in range(M.shape[1]):
        col = M[:, j]
        good = np.isfinite(col)
        r = np.empty(nperm)
        r[~good] = 1.0
        v = col[good]
        order = np.argsort(v, kind='mergesort')
        ranks = np.empty(len(v))
        sv = v[order]
        # 동점은 같은 p 를 받아야 한다 (<= 비율)
        counts = np.searchsorted(sv, sv, side='right')
        ranks[order] = counts
        r[good] = ranks / good.sum()
        P[:, j] = r
    minp_null = P.min(axis=1)
    obs_minp = float(P[obs_idx].min())
    adj = float((minp_null <= obs_minp).sum() / nperm)
    which = names[int(np.argmin(P[obs_idx]))]
    print(f'  {label:<30} 관측 최소 p {obs_minp:.4f} ({which})  '
          f'-> WY 조정 p {adj:.4f}   순열 {nperm:,}')
    return dict(n_perm=int(nperm), n_tests=int(M.shape[1]),
                observed_min_p=obs_minp, which_test=which,
                wy_adjusted_p=adj)


# ------------------------------------------- R1-4 서술자 부트스트랩 거리 구간
def distance_bootstrap(path, label):
    """병합공간 서술자를 복원추출로 재표집해 유효매칭 수에 구간을 붙인다.

    점추정 0/11 은 서술자 선택에 얼마나 기대고 있는가. 같은 부트스트랩으로
    기존 3종 구성과 문헌 포함 8종 구성을 함께 재서, 두 구성의 차까지 낸다.
    차가 모든 복제에서 양수라면 결론은 서술자 선택에 기대지 않는다.
    """
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    keep, _ = pick(df, [c for c in num if c != OMEGA_DERIVED], PROTECTED)
    X = df[keep].to_numpy(float)
    kl = df.compound_class.tolist()
    nm = df.compound_name.tolist()
    ia = [i for i, k in enumerate(kl) if k == 'real_agent']
    SETS = {'incumbent3': ['simulant'],
            'incumbent_plus_lit8': ['simulant', 'literature_simulant']}
    IDX = {lab: ([i for i, k in enumerate(kl) if k in sc],
                 [i for i, k in enumerate(kl)
                  if k not in sc and k != 'real_agent'])
           for lab, sc in SETS.items()}

    rng = np.random.default_rng(SEED)
    counts = {lab: np.empty(N_BOOT, int) for lab in SETS}
    unbounded = {lab: np.empty(N_BOOT, int) for lab in SETS}
    retain = {lab: {nm[a]: 0 for a in ia} for lab in SETS}
    gaps = np.empty((N_BOOT, len(ia)))
    for b in range(N_BOOT):
        c = rng.integers(0, X.shape[1], X.shape[1])
        Xb = X[:, c]
        sd = Xb.std(0, ddof=1)
        ok = sd > 1e-12
        Zb = (Xb[:, ok] - Xb[:, ok].mean(0)) / sd[ok]
        D = zdist(Zb)
        np.fill_diagonal(D, np.inf)
        nn = float(np.median(D.min(1)))
        for lab in SETS:
            sim, out = IDX[lab]
            k = u = 0
            for a in ia:
                ds = D[a, sim].min()
                nearest = ds <= D[a, out].min()
                u += nearest
                if nearest and ds / nn <= 2.0:
                    k += 1
                    retain[lab][nm[a]] += 1
            counts[lab][b] = k
            unbounded[lab][b] = u
        s3, o3 = IDX['incumbent3']
        s8, _ = IDX['incumbent_plus_lit8']
        for t, a in enumerate(ia):
            gaps[b, t] = (D[a, s3].min() - D[a, s8].min()) / nn

    diff = counts['incumbent_plus_lit8'] - counts['incumbent3']
    out = dict(n_boot=N_BOOT)
    for lab in SETS:
        v = counts[lab]
        out[lab] = dict(
            median=int(np.median(v)), mean=float(v.mean()),
            ci=[int(np.percentile(v, 2.5)), int(np.percentile(v, 97.5))],
            frac_zero=float((v == 0).mean()),
            frac_le1=float((v <= 1).mean()),
            unbounded_median=int(np.median(unbounded[lab])),
            retention={k: v2 / N_BOOT for k, v2 in retain[lab].items() if v2})
    out['difference'] = dict(
        median=int(np.median(diff)), mean=float(diff.mean()),
        ci=[int(np.percentile(diff, 2.5)), int(np.percentile(diff, 97.5))],
        frac_positive=float((diff > 0).mean()))
    rows = []
    for t, a in enumerate(ia):
        lo, hi = np.percentile(gaps[:, t], [2.5, 97.5])
        rows.append(dict(agent=nm[a], gap_median=float(np.median(gaps[:, t])),
                         gap_ci=[float(lo), float(hi)],
                         gap_excludes_zero=bool(lo > 0)))
    out['per_agent_gap'] = rows
    out['n_gap_excluding_zero'] = int(sum(r['gap_excludes_zero'] for r in rows))

    print(f'\n  {label}: 서술자 부트스트랩 {N_BOOT}회 (임계값 2.0 포함)')
    for lab in SETS:
        o = out[lab]
        print(f"    {lab:<22} 중위 {o['median']}  평균 {o['mean']:.2f}  "
              f"95% 구간 [{o['ci'][0]}, {o['ci'][1]}]  0건 {100*o['frac_zero']:.1f}%")
    d = out['difference']
    print(f"    차 (8종 - 3종)        중위 {d['median']}  "
          f"95% 구간 [{d['ci'][0]}, {d['ci'][1]}]  "
          f"양수 비율 {100*d['frac_positive']:.1f}%")
    print('    작용제별 매칭 유지율 (8종 구성):')
    for k, v in sorted(out['incumbent_plus_lit8']['retention'].items(),
                       key=lambda x: -x[1])[:5]:
        print(f'      {k:<24} {100*v:5.1f}%')
    return out


if __name__ == '__main__':
    res = {'seed': SEED, 'n_boot': N_BOOT}
    print('=' * 78)
    print('[R2-4] Westfall-Young 단일단계 minP — 18검정의 의존을 순열로 정확히 처리')
    print('=' * 78)
    res['westfall_young'] = {}
    for lab, path, sc in [
            ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
            ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
             ['simulant', 'literature_simulant']),
            ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant'])]:
        res['westfall_young'][lab] = westfall_young(path, sc, lab)

    print()
    print('=' * 78)
    print('[R1-4] 서술자 부트스트랩 — Figure 2 의 거리차에 불확실성을 붙인다')
    print('=' * 78)
    res['distance_bootstrap'] = distance_bootstrap(
        'descriptor_cache_xtb_35.csv', 'panel35')

    with open('outputs/stat_supplement2.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('\n[stat_supplement2] outputs/stat_supplement2.json 기록')
