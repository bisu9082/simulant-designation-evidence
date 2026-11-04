#!/usr/bin/env python3
"""실측 전제 검정의 검정력과 통합 추정 (2026-10-12).

쌍 15개짜리 상관은 약한 설계다. 결론을 내기 전에 (a) 이 표본에서 탐지
가능한 최소 상관이 얼마인지, (b) 서로 다른 물질군·플레이트인 d1 과 d2 를
Fisher z 로 결합하면 무엇이 남는지를 먼저 본다. uv 와 normal 은 같은
플레이트의 두 채널이므로 결합하지 않고 각각 보고한다.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq

S = pd.read_csv('outputs/deltaE_space_test.csv')
# 2026-10-24. 종전에는 N=15, 즉 쌍의 수를 분석 단위로 썼다. 같은 절의 앞쪽
# 문단이 "쌍은 독립 관측이 아니다, 단위는 화합물이다"라고 선언하고 Mantel 로
# 교정해 놓고, 검정력과 Fisher 통합만 쌍을 단위로 돌리고 있었다. 같은 기준을
# 자기 자료에 적용한다. 화합물 6종이 단위이므로 Fisher z 의 가중은 n-3 = 3 이다.
N_PAIRS = 15
N_COMPOUNDS = 6
N = N_COMPOUNDS


def power(rho, n=N, a=0.05):
    z = np.arctanh(rho) * np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - a / 2)
    return stats.norm.cdf(z - zc) + stats.norm.cdf(-z - zc)


# 2026-10-24 2차 교정. Fisher z 근사는 어느 자유도를 넣어도 실제로 돌린
# 검정이 아니다. 본 분석의 셀 검정은 6! = 720 라벨 치환을 전수하는 정확
# Mantel 이고, 그 치환 귀무는 쌍을 독립으로 본 귀무보다 좁다. 좁은 귀무가
# 작은 표본을 상쇄하므로 n=6 을 Fisher z 에 넣은 0.924 는 과도하게 비관적이다.
# 검정 자체를 시뮬레이션한 값(mantel_power.py)을 정본으로 쓰고, 근사값들은
# 비교용으로만 적는다.
print('[검정력] alpha=.05 양측, power .8 에서 탐지 가능한 최소 |rho|')
try:
    _m = json.load(open('outputs/mantel_power.json'))['mde_rho_mantel']
    print(f'    정확 Mantel 시뮬레이션 (정본) -> |rho| >= {_m:.3f}')
except Exception as _e:
    print(f'    [경고] outputs/mantel_power.json 없음 ({_e}). '
          'mantel_power.py 를 먼저 돌릴 것')
print('    [비교용] Fisher z 근사, 검정 자체가 아님:')
for n, unit in ((N_COMPOUNDS, '화합물'), (N_PAIRS, '쌍')):
    print(f'      n={n:>2} {unit} -> |rho| >= '
          f'{brentq(lambda r: power(r, n) - 0.8, .01, .999):.3f}')

print()
print('[통합] 모드별로 d1·d2 를 Fisher z 결합 (물질군·플레이트가 달라 독립)')
rows = []
for mode in ['uv', 'normal']:
    print(f'  --- {mode}')
    for sp in ['structural', 'electronic', 'merged']:
        s = S[(S['mode'] == mode) & (S.space == sp)]
        z = np.arctanh(s.spearman.to_numpy())
        k, w = len(z), N - 3
        zbar = z.mean()
        se = 1 / np.sqrt(w * k)
        r, lo, hi = np.tanh(zbar), np.tanh(zbar - 1.96 * se), np.tanh(zbar + 1.96 * se)
        p = 2 * (1 - stats.norm.cdf(abs(zbar) / se))
        print(f'      {sp:<11} rho={r:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  p={p:.3f}')
        rows.append(dict(mode=mode, space=sp, pooled_rho=float(r),
                         ci_lo=float(lo), ci_hi=float(hi), p=float(p)))

print()
print('[대조] 전자 공간이 구조 공간을 이기는 셀이 있는가')
wins = 0
for (ds, mode), g in S.groupby(['dataset', 'mode']):
    e = float(g[g.space == 'electronic'].spearman.iloc[0])
    st = float(g[g.space == 'structural'].spearman.iloc[0])
    wins += e > st
    print(f'    {ds} {mode:<7} 전자 {e:+.3f}  구조 {st:+.3f}  '
          f'-> {"전자 우세" if e > st else "구조 우세"}')
print(f'    전자 우세 {wins}/4 셀')

pd.DataFrame(rows).to_csv('outputs/deltaE_pooled.csv', index=False)
print('\n[done] outputs/deltaE_pooled.csv')
