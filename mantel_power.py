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
N_REP_NULL = 20000
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
    rejects, rejects_bh, rhos = 0, 0, np.empty(N_REP)
    for r in range(N_REP):
        X = rng.normal(size=(N, 3))
        D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
        d = D[IU]
        v = d + rng.normal(scale=noise * d.std(), size=NP_)
        p, obs = mantel_p(rankdata(d), rankdata(v))
        rhos[r] = obs
        rejects += (p <= ALPHA)
        # 2026-11-06 (9차 R3-2). 논문이 실제로 판정에 쓰는 임계는 명목
        # alpha 가 아니라 12셀 BH 다. 지정검정에서는 이미 보정 임계의 MDE
        # 를 싣고 있었는데(mde.json 의 bh18) 헤드라인인 Mantel 에만
        # 적용하지 않았다. 같은 격자에서 두 임계의 검정력을 함께 센다.
        rejects_bh += (p <= ALPHA_BH)
    return float(rhos.mean()), rejects / N_REP, rejects_bh / N_REP


def run_null(n_rep, rng):
    """귀무 조건의 1종오류율 (2026-11-06 신설).

    왜 신설했는가. SI 가 'the type I rate of the simulated test at 0.049'
    이라고 적고 있었는데, 이 스크립트에는 귀무 조건이 없었고 어느 예치물에도
    그 수가 없었다. 수행하지 않은 결과를 적은 것이므로 [C-3] 위반이다.
    문장을 지우는 대신 실제로 돌려서 나온 수로 교체한다.

    귀무의 모양. '두 거리행렬이 무관하다'이므로 잡음을 쌍에 직접 얹으면 안
    된다. 그건 거리행렬이 아니다. 독립한 배치 Y 에서 두 번째 거리행렬을
    만들어야 삼각부등식 같은 거리행렬 고유의 종속구조가 보존된다.

    기대값. 검정이 720 전수 치환이고 항등치환을 귀무분포에 포함하므로
    p >= 1/720 이고 검정은 보수적으로 타당하다. 즉 1종오류율은 0.05 를
    넘지 않아야 하고, 넘으면 절차에 버그가 있는 것이다.
    """
    rejects = 0
    for _ in range(n_rep):
        X = rng.normal(size=(N, 3))
        Y = rng.normal(size=(N, 3))
        d = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))[IU]
        v = np.sqrt(((Y[:, None] - Y[None]) ** 2).sum(-1))[IU]
        p, _ = mantel_p(rankdata(d), rankdata(v))
        rejects += (p <= ALPHA)
    rate = rejects / n_rep
    se = float(np.sqrt(rate * (1 - rate) / n_rep))
    return rate, se


BH_FAMILY = 12                    # 서술자 Mantel 셀 수 (rollback_c.deltaE.mantel_bh)
ALPHA_BH = 0.05 / BH_FAMILY       # BH 에서 최소순위가 통과해야 하는 임계


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    print('잡음폭   평균 관측 |rho|   기각률 @a=0.05   기각률 @a=0.05/12')
    for noise in (0.10, 0.20, 0.30, 0.40, 0.50, 0.65, 0.85, 1.10, 1.50):
        rho, pw, pw_bh = run(noise, rng)
        rows.append(dict(noise=noise, mean_abs_rho=rho, power=pw, power_bh=pw_bh))
        print(f'  {noise:.2f}        {rho:.3f}          {pw:.3f}            {pw_bh:.3f}')

    print()
    t1, t1_se = run_null(N_REP_NULL, rng)
    print(f'[mantel_power] 귀무 1종오류율 = {t1:.4f} +- {t1_se:.4f} '
          f'(명목 {ALPHA}, 반복 {N_REP_NULL})')
    if t1 > ALPHA + 3 * t1_se:
        print('[mantel_power] 경고: 명목 수준을 유의하게 초과한다. 절차 재점검.')
    print()

    rows.sort(key=lambda r: r['mean_abs_rho'])

    def interp(field):
        for a, b in zip(rows, rows[1:]):
            if a[field] <= 0.8 <= b[field]:
                t = (0.8 - a[field]) / (b[field] - a[field])
                return a['mean_abs_rho'] + t * (b['mean_abs_rho']
                                                - a['mean_abs_rho'])
        return None

    mde = interp('power')
    mde_bh = interp('power_bh')
    if mde_bh is None:
        print('[mantel_power] 보정 임계의 80% 구간이 격자 밖이다 — '
              f'최대 격자점에서도 검정력 {rows[-1]["power_bh"]:.3f}. '
              '격자를 넓히거나 상한으로 보고할 것.')
    else:
        print(f'[mantel_power] BH 임계(a={ALPHA_BH:.5f}) 80% 검정력 '
              f'최소탐지 |rho| = {mde_bh:.3f}')
    print()
    if mde is None:
        print('[mantel_power] 80% 구간이 격자 밖이다. 격자를 넓힐 것.')
    else:
        print(f'[mantel_power] 정확 Mantel 80% 검정력 최소탐지 |rho| = {mde:.3f}')
    json.dump(dict(seed=SEED, n_compounds=N, n_perm=len(PAIRPERM),
                   n_rep=N_REP, alpha=ALPHA, grid=rows,
                   mde_rho_mantel=(None if mde is None else round(float(mde), 3)),
                   bh_family=BH_FAMILY, alpha_bh=ALPHA_BH,
                   mde_rho_mantel_bh=(None if mde_bh is None
                                      else round(float(mde_bh), 3)),
                   null=dict(n_rep=N_REP_NULL,
                             type_i_rate=round(float(t1), 5),
                             type_i_se=round(float(t1_se), 5),
                             nominal_alpha=ALPHA,
                             p_min_attainable=1.0 / len(PAIRPERM),
                             note='독립한 두 배치에서 거리행렬 2개를 만들어 '
                                  '실제 절차(720 전수 치환, 양측)를 그대로 돌린 '
                                  '기각률. 항등치환을 귀무분포에 포함하므로 '
                                  'p >= 1/720 이고 검정은 보수적으로 타당하다.')),
              open('outputs/mantel_power.json', 'w'), indent=1)
    print('[done] outputs/mantel_power.json')


if __name__ == '__main__':
    main()
