#!/usr/bin/env python3
"""
simulant / GATE 4 — Sanity Check ①② (③ 교차조건은 analysis_main.py 의 3조건 ablation)

① 단조성  : 서술자를 nested 로 추가하면 표준화 유클리드 거리는 단조 증가해야 한다.
             위반 시 표준화·결측 처리 버그. d/sqrt(k) 는 안정해야 한다.
② Baseline: 모사체 3 vs 간섭물 14 라는 군 크기 불균형 자체가 '간섭물이 더 가깝다'를
             만들어낸다(14번 뽑은 최소값이 3번 뽑은 최소값보다 작기 마련).
             17종 중 3종을 모사체로 두는 모든 라벨링 C(17,3)=680 가지를 전수 열거해
             귀무분포를 만들고, 실제 지정 모사체의 순위로 정확 p값을 낸다.
"""
import json
from itertools import combinations
import numpy as np
import pandas as pd

SEED = 42
df = pd.read_csv('descriptor_cache_xtb.csv')
S = json.load(open('outputs/experiment_summary.json'))
CON = {c['condition']: c for c in S['conditions']}
klass = df.compound_class.to_numpy()
i_real = np.flatnonzero(klass == 'real_agent')
i_pool = np.flatnonzero(klass != 'real_agent')              # 모사체+간섭물 17종
i_sim_true = set(np.flatnonzero(klass == 'simulant').tolist())
print(f'[setup] 진본 {len(i_real)}종 | 비작용제 풀 {len(i_pool)}종 '
      f'(지정 모사체 {len(i_sim_true)}종)')


def zmat(cols):
    X = df[cols].to_numpy(float)
    return (X - X.mean(0)) / X.std(0, ddof=1)


def dmat(Z):
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D


# =============================================================== ① 단조성
print('\n=== Sanity ① 단조성 (nested 서술자 추가) ===')
keep = CON['C_merged']['kept']
ks, means, norm = [], [], []
for k in range(4, len(keep) + 1, 4):
    Z = zmat(keep[:k])
    D = dmat(Z)
    fin = D[np.isfinite(D)]
    ks.append(k); means.append(fin.mean()); norm.append(fin.mean() / np.sqrt(k))
diffs = np.diff(means)
ok1 = bool((diffs > 0).all())
print('  k      : ' + ' '.join(f'{k:>6d}' for k in ks))
print('  mean d : ' + ' '.join(f'{m:6.2f}' for m in means))
print('  d/sqrt(k): ' + ' '.join(f'{m:6.3f}' for m in norm))
print(f'  단조 증가: {ok1}   (증가분 최소 {diffs.min():+.4f})')
print(f'  d/sqrt(k) 범위 {min(norm):.3f}~{max(norm):.3f} '
      f'(변동 {100*(max(norm)-min(norm))/np.mean(norm):.1f}%)')

# ====================================================== ② Baseline 타당성
print('\n=== Sanity ② Baseline 타당성 (C(17,3)=680 전수 순열) ===')
ALL = list(combinations(range(len(i_pool)), 3))
assert len(ALL) == 680
out2 = {}
for cond in ['A_rdkit', 'B_electronic', 'C_merged']:
    D = dmat(zmat(CON[cond]['kept']))
    Dsub = D[np.ix_(i_real, i_pool)]                        # 11 x 17
    cnt_null, diff_null = np.empty(680), np.empty(680)
    true_idx = None
    for t, combo in enumerate(ALL):
        m = np.zeros(len(i_pool), bool); m[list(combo)] = True
        ds = Dsub[:, m].min(1)
        di = Dsub[:, ~m].min(1)
        cnt_null[t] = (di < ds).sum()
        diff_null[t] = (ds - di).mean()
        if {int(i_pool[c]) for c in combo} == i_sim_true:
            true_idx = t
    assert true_idx is not None
    c_obs, d_obs = cnt_null[true_idx], diff_null[true_idx]
    # 실제 모사체가 '무작위 3종보다 더 좋은(=차이가 더 작은)' 쪽의 정확 p값
    p_exact = float((diff_null <= d_obs).sum() / 680)
    print(f'\n  [{cond}]')
    print(f'    귀무 기대 간섭물우세 수 : {cnt_null.mean():.2f} / 11  '
          f'(군크기 14:3 이 만드는 기저값)')
    print(f'    실제 지정 모사체        : {c_obs:.0f} / 11')
    print(f'    귀무 평균차 분포        : {diff_null.mean():+.4f} '
          f'[{np.percentile(diff_null,2.5):+.4f}, {np.percentile(diff_null,97.5):+.4f}]')
    print(f'    실제 평균차             : {d_obs:+.4f}  '
          f'-> 680개 중 {int((diff_null<=d_obs).sum())}번째, 정확 p = {p_exact:.4f}')
    verdict = ('지정 모사체가 무작위 3종보다 유의하게 가깝다'
               if p_exact < 0.05 else
               '지정 모사체가 무작위로 고른 3종과 구별되지 않는다')
    print(f'    판정: {verdict}')
    np.save(f'outputs/_null_{cond}.npy', diff_null)
    out2[cond] = dict(null_mean_count=float(cnt_null.mean()),
                      observed_count=float(c_obs),
                      null_mean_diff=float(diff_null.mean()),
                      null_ci=[float(np.percentile(diff_null, 2.5)),
                               float(np.percentile(diff_null, 97.5))],
                      observed_diff=float(d_obs), p_exact=p_exact,
                      n_permutations=680, verdict=verdict)

res = dict(sanity_monotonicity=dict(passed=ok1, k_grid=ks,
                                    mean_distance=[float(m) for m in means],
                                    d_over_sqrt_k=[float(m) for m in norm],
                                    min_increment=float(diffs.min())),
           sanity_baseline=out2)
with open('outputs/sanity_checks.json', 'w') as f:
    json.dump(res, f, indent=2, ensure_ascii=False)
print('\n[done] outputs/sanity_checks.json')
