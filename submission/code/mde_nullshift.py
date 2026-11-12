#!/usr/bin/env python3
"""MDE 재계산 — 이동된 거리행렬에서 귀무를 재생성한다 (2026-11-12 신설).

왜 필요한가. mde.py 는 귀무분포를 수정하지 않은 거리행렬에서 한 번 만들고
delta 루프 밖에 둔다(mde.py:47-68, 92-96). docstring 이 그 사실을 선언하고
해석을 조건부로 한정하지만, 원고는 그 결과를 조건 없는 검정력 진술로
쓴다 ("The test is therefore sensitive to effects smaller than the one at
issue").

그 절차가 불가능한 값을 낸다. outputs/mde.json 의 panel35 curve 151점 중
94점이 p=0.0 인데, C(24,3)=2024 전수 치환의 p 하한은 1/2024=4.94e-4 다.
지정 3종만 당기고 귀무를 그대로 두면 관측값이 귀무 전체 아래로 빠지기
때문이다.

이 검정의 귀무는 같은 거리행렬에서 라벨을 바꿔 만든다. 2024개 지정 중
C(21,3)=1330개(65.7%)가 기존 3종을 비교풀에 두므로, 3종이 가까워지면
그 1330개의 d_intf 도 같이 줄고 귀무가 같은 방향으로 이동한다. 이동을
반영해야 검정력이 과대추정되지 않는다.

여기서는 delta 마다 이동된 행렬 전체에서 귀무를 다시 만든다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from mde import CASES, ALPHA, GRID, merged_D

BH18 = 0.05 / 18.0


def shift_matrix(D, i_real, i_sim, amount):
    """지정 3종과 작용제 사이 거리만 amount 만큼 줄인 행렬 사본."""
    Ds = D.copy()
    for a in i_real:
        for s in i_sim:
            v = max(Ds[a, s] - amount, 1e-9)
            Ds[a, s] = v
            Ds[s, a] = v
    return Ds


def null_and_obs(D, i_real, i_pool, true_combo, ALL):
    """이 행렬에서 전수 귀무와 실제 지정의 통계량·p 를 낸다."""
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    npool = len(i_pool)
    out = np.empty(len(ALL))
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        out[a0:a0 + b] = (ds - di).mean(0)
    ti = ALL.index(true_combo)
    p = float((out <= out[ti]).sum() / len(out))
    return float(out[ti]), p, len(out)


def run(label, path, sim_cls):
    df, D, nk = merged_D(path)
    kl = df.compound_class.tolist()
    name = df.compound_name.tolist()
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    i_other = [i for i in i_pool if i not in i_sim]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    floor = 1.0 / len(ALL)

    rows, mde_a, mde_b = [], None, None
    for d in GRID:
        Dm = shift_matrix(D, i_real, i_sim, d * nn_med)
        obs, p, n = null_and_obs(Dm, i_real, i_pool, true_combo, ALL)
        rows.append(dict(delta=float(d), observed=obs, p=p))
        if mde_a is None and p < ALPHA:
            mde_a = float(d)
        if mde_b is None and p < BH18:
            mde_b = float(d)

    # 실제 격차
    ds0 = np.array([min(D[a, s] for s in i_sim) for a in i_real])
    di0 = np.array([min(D[a, o] for o in i_other) for a in i_real])
    gap = (ds0 - di0) / nn_med
    gmed = float(np.median(gap))
    per_agent = {name[i_real[k]]: float(gap[k]) for k in range(len(i_real))}

    n0 = sum(1 for r in rows if r['p'] == 0.0)
    print(f'\n### {label}  ({nk}차원, 귀무 {len(ALL):,}, 전형최근접 {nn_med:.3f})')
    print(f'  p 하한 1/{len(ALL)} = {floor:.6f}   curve 에서 p=0.0 인 점 {n0}')
    print(f'  delta=0 에서 p = {rows[0]["p"]:.4f}')
    print(f'  alpha=0.05 MDE : {mde_a}')
    print(f'  BH18 MDE       : {mde_b}')
    print(f'  실제 중위 격차 : {gmed:.4f}')
    for lv, m in (('alpha', mde_a), ('bh18', mde_b)):
        if m:
            print(f'    {lv}: 중위격차/MDE = {gmed/m:.3f} '
                  f'({"검출 가능" if gmed/m >= 1 else "검출 하한 아래"})')
            print(f'    {lv}: MDE 아래인 작용제 '
                  f'{sum(1 for v in per_agent.values() if v < m)}/11')
    return dict(n_dims=nk, n_permutations=len(ALL), p_floor=floor,
                nn_median=nn_med, p_at_delta0=rows[0]['p'],
                mde_alpha=mde_a, mde_bh18=mde_b,
                gap_median=gmed, gap_per_agent=per_agent,
                n_below_mde_alpha=sum(1 for v in per_agent.values()
                                      if mde_a and v < mde_a),
                n_below_mde_bh18=sum(1 for v in per_agent.values()
                                     if mde_b and v < mde_b),
                ratio_gap_to_mde_alpha=gmed / mde_a if mde_a else None,
                ratio_gap_to_mde_bh18=gmed / mde_b if mde_b else None,
                curve=rows)


def main():
    res = {'generated_by': 'mde_nullshift.py',
           'note': 'delta 마다 이동된 거리행렬에서 귀무를 재생성한 MDE. '
                   'mde.py 는 귀무를 고정한다.'}
    for lab, path, cls in CASES:
        res[lab] = run(lab, path, cls)
    with open('outputs/mde_nullshift.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('\n[mde_nullshift] outputs/mde_nullshift.json 기록')


if __name__ == '__main__':
    main()
