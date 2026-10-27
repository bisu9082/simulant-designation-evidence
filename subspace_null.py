#!/usr/bin/env python3
"""부분공간 대조의 귀무 보정 (2026-10-27 신설).

왜. 본문은 "무작위 3서술자 부분공간 2000개 중 36.0% 가 p<0.05 에 이르고
관측값은 그 분포의 25.9 백분위이므로 둘 다 평범한 draw 다" 라고 쓴다.
그런데 36.0% 가 귀무 하에서 얼마여야 하는지는 계산하지 않았다.

고정된 부분공간에서 p 는 C(17,3)=680 개 지정에 대한 순위이므로, 지정이
임의라면 p 는 균일분포이고 P(p<0.05) 의 기대값은 0.05 다. 따라서 36% 는
평범한 값이 아니라 귀무율의 약 7배다.

여기서는 그것을 시늄레이션이 아니라 전수로 확인한다. 한 부분공간을 한 번
열거하면 680개 지정 전부의 통계량이 나오므로, 같은 비용으로 모든 지정의
p 를 읽을 수 있다. 실제 지정의 36% 가 680개 지정의 분포에서 어디에
앉는지 그대로 낸다.

출력: outputs/subspace_null.json
"""
import json
import sys
from itertools import combinations

import numpy as np
import pandas as pd

from spaces import pick, PROTECTED

PANEL = 'descriptor_cache_xtb.csv'
N_SUB = 500
SEED = 49374
ELEC = {'dipole_moment_D', 'E_HOMO_eV', 'E_LUMO_eV',
        'electrophilicity_index_eV'}
META = {'n_conformers', 'conformer_w_max'}


def main():
    df = pd.read_csv(PANEL)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n]
    keep, _ = pick(df, [c for c in num if c not in ELEC and c not in META],
                   PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)

    kl = df.compound_class.to_numpy()
    i_real = np.flatnonzero(kl == 'real_agent')
    i_pool = np.flatnonzero(kl != 'real_agent')
    i_sim = set(np.flatnonzero(kl == 'simulant').tolist())
    true = tuple(sorted(list(i_pool).index(j) for j in i_sim))
    allc = list(combinations(range(len(i_pool)), 3))
    assert len(allc) == 680
    ti = allc.index(true)
    C = np.asarray(allc, np.int32)
    npool = len(i_pool)
    M = np.zeros((680, npool), bool)
    M[np.arange(680)[:, None], C] = True

    rng = np.random.default_rng(SEED)
    hit = np.zeros(680, int)          # 지정별 p<0.05 인 부분공간 수
    print(f'[부분공간] {N_SUB}회, 좌표 {Z.shape[1]}개, 지정 680개 동시 평가',
          flush=True)
    for t in range(N_SUB):
        cols = rng.choice(Z.shape[1], 3, replace=False)
        Zx = Z[:, cols]
        D = np.sqrt(((Zx[:, None, :] - Zx[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(D, np.inf)
        Ds = D[np.ix_(i_real, i_pool)]
        A3 = Ds[:, None, :]
        ds = np.where(M[None], A3, np.inf).min(2)
        di = np.where(~M[None], A3, np.inf).min(2)
        stat = (ds - di).mean(0)                      # (680,)
        # 각 지정의 일측 p = 자기 통계량 이하인 지정의 비율
        order = stat.argsort(kind='mergesort')
        ranks = np.empty(680, float)
        srt = stat[order]
        # 동률을 포함해 '이하' 개수를 센다
        cnt = np.searchsorted(srt, srt, side='right')
        ranks[order] = cnt / 680.0
        hit += (ranks < 0.05)
        if (t + 1) % 100 == 0:
            print(f'  {t+1}/{N_SUB}  실제 지정 누적 '
                  f'{100*hit[ti]/(t+1):.1f}%', flush=True)

    frac = hit / N_SUB
    actual = float(frac[ti])
    out = dict(
        n_subspaces=N_SUB, n_designations=680, seed=SEED,
        n_coords=int(Z.shape[1]),
        actual_designation_frac=actual,
        null_mean_frac=float(frac.mean()),
        null_median_frac=float(np.median(frac)),
        null_p95_frac=float(np.percentile(frac, 95)),
        null_max_frac=float(frac.max()),
        rank_of_actual=int((frac >= actual).sum()),
        p_of_actual=float((frac >= actual).mean()),
    )
    print('\n실제 지정 (DMMP/DCM/DPGME) p<0.05 비율: '
          f'{100*actual:.1f}%')
    print(f'680개 지정 전체 평균: {100*out["null_mean_frac"]:.1f}%  '
          f'중위 {100*out["null_median_frac"]:.1f}%  '
          f'95백분위 {100*out["null_p95_frac"]:.1f}%  '
          f'최대 {100*out["null_max_frac"]:.1f}%')
    print(f'실제 지정의 순위: {out["rank_of_actual"]}/680  '
          f'(p = {out["p_of_actual"]:.4f})')
    json.dump(out, open('outputs/subspace_null.json', 'w'), indent=1)
    print('\n[done] outputs/subspace_null.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
