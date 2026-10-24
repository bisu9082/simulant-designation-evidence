#!/usr/bin/env python3
"""Mantel 검정의 실제 검정력 (2026-10-24 신설).

왜 필요한가. 실측 어레이 절의 최소탐지효과를 지금까지 Fisher z 근사로 적었다.
처음에는 쌍 15개를 단위로 넣어 |rho| >= 0.669 가 나왔고, 분석 단위가 화합물
이라는 지적을 받아 n = 6 으로 고쳐 0.924 가 나왔다. 그런데 실제로 돌린 검정은
Fisher z 가 아니라 6! = 720 개 라벨 치환을 전수하는 정확 Mantel 이다. 근사의
자유도를 고쳐 봐야 그 근사가 검정 자체와 다르면 의미가 없다.

그래서 검정 자체를 시뮬레이션한다. 화합물 좌표에서 서술자 거리행렬을 만들고,
패턴 거리는 그 단조관계에 잡음을 더해 만든다. 잡음 폭으로 연관의 세기를
조절하고, 각 반복에서 실제 쓰는 절차(쌍 Spearman + 720 전수치환, 양측)를
그대로 돌려 기각률을 센다.

속도. 치환은 화합물 라벨에 걸리므로 그것이 유도하는 '쌍 치환'은 D 와 무관하게
한 번만 계산하면 된다. 순위벡터는 묶음이 없으면 1..15 의 치환이라 Pearson 이
Spearman 과 같고, 주변분포가 고정이므로 상관은 내적의 1차식이 된다.
"""
import json
import numpy as np
from itertools import permutations
from scipy.stats import rankdata

SEED = 0xC0DE
N_REP = 2000
ALPHA = 0.05
N = 6

IU = np.triu_indices(N, 1)
NP_ = len(IU[0])                                   # 15

_pos = -np.ones((N, N), int)
for _k, (_i, _j) in enumerate(zip(*IU)):
    _pos[_i, _j] = _pos[_j, _i] = _k
PAIRPERM = np.array([[_pos[p[i], p[j]] for i, j in zip(*IU)]
                     for p in permutations(range(N))])     # (720, 15)


def spearman_vec(a_rank, b_rank):
    """순위가 1..n 의 치환일 때의 Spearman. 주변분포 고정이므로 내적의 1차식.

    a_rank 는 (..., n) 이다. len() 으로 n 을 잡으면 (720, 15) 에서 720 을
    집어 버린다. 2026-10-24 에 그 상태로 귀무 1종오류가 0.85 로 나왔다.
    """
    n = a_rank.shape[-1]
    m = (n + 1) / 2.0
    s = ((np.arange(1, n + 1) - m) ** 2).sum()
    return ((a_rank - m) * (b_rank - m)).sum(-1) / s


def mantel_p(dr, pr):
    """dr, pr 은 쌍 거리의 순위벡터. 720 치환 전수, 양측(|rho|)."""
    obs = abs(spearman_vec(dr, pr))
    null = np.abs(spearman_vec(dr[PAIRPERM], pr))
    return float((null >= obs - 1e-12).mean()), float(obs)


def run(noise, rng):
    rejects, rhos = 0, np.empty(N_REP)
    for r in range(N_REP):
        X = rng.normal(size=(N, 3))
        D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
        d = D[IU]
        v = d + rng.normal(scale=noise * d.std(), size=NP_)
        p, obs = mantel_p(rankdata(d), rankdata(v))
        rhos[r] = obs
        rejects += (p <= ALPHA)
    return float(rhos.mean()), rejects / N_REP


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    print('잡음폭   평균 관측 |rho|   기각률(검정력)')
    for noise in (0.10, 0.20, 0.30, 0.40, 0.50, 0.65, 0.85, 1.10, 1.50):
        rho, pw = run(noise, rng)
        rows.append(dict(noise=noise, mean_abs_rho=rho, power=pw))
        print(f'  {noise:.2f}        {rho:.3f}            {pw:.3f}')

    rows.sort(key=lambda r: r['mean_abs_rho'])
    mde = None
    for a, b in zip(rows, rows[1:]):
        if a['power'] <= 0.8 <= b['power']:
            t = (0.8 - a['power']) / (b['power'] - a['power'])
            mde = a['mean_abs_rho'] + t * (b['mean_abs_rho'] - a['mean_abs_rho'])
            break
    print()
    if mde is None:
        print('[mantel_power] 80% 구간이 격자 밖이다. 격자를 넓힐 것.')
    else:
        print(f'[mantel_power] 정확 Mantel 80% 검정력 최소탐지 |rho| = {mde:.3f}')
    json.dump(dict(seed=SEED, n_compounds=N, n_perm=len(PAIRPERM),
                   n_rep=N_REP, alpha=ALPHA, grid=rows,
                   mde_rho_mantel=(None if mde is None else round(float(mde), 3))),
              open('outputs/mantel_power.json', 'w'), indent=1)
    print('[done] outputs/mantel_power.json')


if __name__ == '__main__':
    main()
