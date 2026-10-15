#!/usr/bin/env python3
"""실측 전제 검정의 검정력과 통합 추정 (2026-10-12).

쌍 15개짜리 상관은 약한 설계다. 결론을 내기 전에 (a) 이 표본에서 탐지
가능한 최소 상관이 얼마인지, (b) 서로 다른 물질군·플레이트인 d1 과 d2 를
Fisher z 로 결합하면 무엇이 남는지를 먼저 본다. uv 와 normal 은 같은
플레이트의 두 채널이므로 결합하지 않고 각각 보고한다.
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq

S = pd.read_csv('outputs/deltaE_space_test.csv')
N = 15


def power(rho, n=N, a=0.05):
    z = np.arctanh(rho) * np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - a / 2)
    return stats.norm.cdf(z - zc) + stats.norm.cdf(-z - zc)


print('[검정력] alpha=.05 양측, power .8 에서 탐지 가능한 최소 |rho|')
for n in (15, 30):
    print(f'    n={n:>2}쌍 -> |rho| >= {brentq(lambda r: power(r, n) - 0.8, .01, .99):.3f}')

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
