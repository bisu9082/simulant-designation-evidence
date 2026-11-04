#!/usr/bin/env python3
"""검정통계량 선택에 대한 순열검정 민감도 (2026-10-12 신설, GATE 8 리뷰어 C M3).

지적. 원고는 평균차 mean(d_sim - d_intf) 하나로 p=0.028 을 보고한다. 그런데
이 통계량은 정규화되지 않은 양이어서 거리 12~14 구간의 원거리 작용제가
거리 4~7 구간과 같은 가중을 받는다. 원고 스스로 2.4절에서 "최근접 모사체가
11.8 거리에 있으면 매칭이 아니다" 라고 쓰는데, 그 VX 가 유의성에 GB 에
버금가게 기여한다. 즉 유의성을 만드는 통계량과 매칭을 판정하는 기준이
서로 다른 양을 재고 있다.

여기서는 귀무(전수 순열)와 거리행렬을 고정하고 통계량만 바꿔 전부 보고한다.
헤드라인을 가장 유리한 통계량 하나로 쓰지 않기 위해서다.

통계량
  mean_diff     mean(d_sim - d_intf)            원고가 쓰던 것
  median_diff   median(d_sim - d_intf)          이상치 가중 제거
  mean_logratio mean(log(d_sim / d_intf))       척도 불변
  mean_dsim     mean(d_sim)                     간섭물 비교 없이
  n_closer      모사체가 더 가까운 작용제 수      순서통계량
  n_valid       유효매칭 수 (저자 자신의 기준)    매칭 판정과 같은 양
"""
import argparse
import json
from itertools import combinations

import numpy as np
import pandas as pd

from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM,
                           OMEGA_DERIVED, PROTECTED, pick)

STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']


def statistics(ds, di, nn_med):
    """작용제별 (d_sim, d_intf) 에서 통계량 6종. 작을수록 '모사체가 좋다'."""
    return dict(
        mean_diff=float(np.mean(ds - di)),
        median_diff=float(np.median(ds - di)),
        mean_logratio=float(np.mean(np.log(ds / di))),
        mean_dsim=float(np.mean(ds)),
        # 아래 둘은 클수록 좋으므로 부호를 뒤집어 방향을 통일한다
        n_closer=float(-np.sum(ds < di)),
        n_valid=float(-np.sum((ds / nn_med <= 2.0) & (ds <= di))),
    )


def run(path, sim_classes, label):
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    SP = {'A_rdkit': [c for c in num if c not in ELECTRONIC_28],
          'B_electronic': [c for c in num if c in ELECTRONIC_DIST],
          'C_merged': [c for c in num if c != OMEGA_DERIVED]}
    kl = df.compound_class.tolist()
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    i_int = [i for i, k in enumerate(kl) if k == 'interferent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(combinations(range(len(i_pool)), len(true_combo)))

    print(f'\n### {label}')
    print(f'  {n}종 (작용제 {len(i_real)}, 비작용제 {len(i_pool)}, '
          f'모사체 {len(i_sim)}, 간섭물 {len(i_int)})  '
          f'귀무 C({len(i_pool)},{len(true_combo)}) = {len(ALL):,}')
    out = {}
    for cond, cols in SP.items():
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(D, np.inf)
        nn_med = float(np.median(D.min(1)))
        Ds = D[np.ix_(i_real, i_pool)]
        # 귀무는 원고와 같은 규약을 쓴다. 한 라벨링에서 선택된 k 종이 모사체,
        # 비작용제 풀의 나머지가 간섭물이다. 간섭물 집합을 실제 집합에 고정하면
        # 다른 귀무가 되고 p 가 달라진다 (28종 구조공간에서 0.028 대 0.085).
        # 벡터화. C(24,8)=735,471 조합을 파이썬 루프로 돌면 분 단위가 된다.
        # 조합 배열을 블록으로 잘라 Ds[:, combo] 의 열 최소를 한 번에 구한다.
        C = np.asarray(ALL, dtype=np.int32)              # (nperm, k)
        null = {s: np.empty(len(ALL)) for s in STATS}
        npool = len(i_pool)
        BLK = 8000
        for a0 in range(0, len(C), BLK):
            cb = C[a0:a0 + BLK]                          # (b, k)
            b = len(cb)
            m = np.zeros((b, npool), bool)
            m[np.arange(b)[:, None], cb] = True
            A3 = Ds[:, None, :]                          # (n_real, 1, npool)
            ds = np.where(m[None, :, :], A3, np.inf).min(2)   # (n_real, b)
            di = np.where(~m[None, :, :], A3, np.inf).min(2)
            null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
            null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
            null['mean_logratio'][a0:a0 + b] = np.log(ds / di).mean(0)
            null['mean_dsim'][a0:a0 + b] = ds.mean(0)
            null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
            null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
        ti = ALL.index(true_combo)
        obs = {s: float(null[s][ti]) for s in STATS}
        row = {}
        for s in STATS:
            p = float((null[s] <= obs[s]).sum() / len(ALL))
            row[s] = dict(observed=obs[s], p_one_sided=p)
        out[cond] = dict(n_kept=len(keep), nn_median=nn_med,
                         n_perm=len(ALL), statistics=row)
        print(f'  {cond}')
        for s in STATS:
            v = obs[s]
            disp = f'{-v:.0f}' if s.startswith('n_') else f'{v:+.4f}'
            print(f'      {s:<14} 관측 {disp:>9}   p = {row[s]["p_one_sided"]:.4f}')
    return out


def jackknife(path, sim_classes, label):
    """작용제 1종씩 빼고 mean_diff 의 p 가 얼마나 흔들리는지."""
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    kl = df.compound_class.tolist()
    name = df.compound_name.tolist()
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    i_int = [i for i, k in enumerate(kl) if k == 'interferent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(combinations(range(len(i_pool)), len(true_combo)))
    print(f'\n### {label} — 작용제 1종 제외 jackknife (mean_diff)')
    res = {}
    for cond, cols in [('A_rdkit', [c for c in num if c not in ELECTRONIC_28]),
                       ('B_electronic', [c for c in num if c in ELECTRONIC_DIST]),
                       ('C_merged', [c for c in num if c != OMEGA_DERIVED])]:
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(D, np.inf)
        ps = []
        for drop in range(len(i_real)):
            keep_r = [r for j, r in enumerate(i_real) if j != drop]
            Ds = D[np.ix_(keep_r, i_pool)]
            C = np.asarray(ALL, dtype=np.int32)
            null = np.empty(len(ALL))
            npool = len(i_pool)
            BLK = 8000
            for a0 in range(0, len(C), BLK):
                cb = C[a0:a0 + BLK]; b = len(cb)
                m = np.zeros((b, npool), bool)
                m[np.arange(b)[:, None], cb] = True
                A3 = Ds[:, None, :]
                dsb = np.where(m[None, :, :], A3, np.inf).min(2)
                dib = np.where(~m[None, :, :], A3, np.inf).min(2)
                null[a0:a0 + b] = (dsb - dib).mean(0)
            obs = null[ALL.index(true_combo)]
            ps.append((name[i_real[drop]], float((null <= obs).sum() / len(ALL))))
        lo, hi = min(p for _, p in ps), max(p for _, p in ps)
        worst = max(ps, key=lambda t: t[1])
        print(f'  {cond:<13} p 범위 {lo:.4f} ~ {hi:.4f}   최악 = {worst[0]} 제외 시 {worst[1]:.4f}')
        res[cond] = dict(range=[lo, hi], worst_agent=worst[0], worst_p=worst[1],
                         all=[dict(dropped=a, p=p) for a, p in ps])
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel35', default='descriptor_cache_xtb_35.csv')
    ap.add_argument('--panel28', default='descriptor_cache_xtb.csv')
    A = ap.parse_args()
    out = {}
    out['panel28_designated'] = run(A.panel28, ('simulant',),
                                    '28종 패널 — 지정 모사체 3종')
    out['panel35_designated'] = run(A.panel35, ('simulant',),
                                    '35종 확장 — 지정 모사체 3종')
    out['panel35_all'] = run(A.panel35, ('simulant', 'literature_simulant'),
                             '35종 확장 — 지정 + 문헌 모사체 8종')
    out['jackknife_panel35_all'] = jackknife(
        A.panel35, ('simulant', 'literature_simulant'),
        '35종 확장 — 지정 + 문헌 모사체 8종')
    json.dump(out, open('outputs/permutation_statistics.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n[done] outputs/permutation_statistics.json')
