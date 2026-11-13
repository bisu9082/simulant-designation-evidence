#!/usr/bin/env python3
"""외부 검증 — Davidson et al. 2020 (ACS Sens. 5, 1102) 의 측정된 군집과 대조.

왜 필요한가. 패널 6인 중 4인이 같은 지점을 짚었다: 분석화학 저널에
'새 측정이 없다'. 우리가 측정을 추가할 수는 없다. 그러나 진본 작용제
14종을 변색 센서 어레이로 측정한 논문이 이미 있고, 우리 패널과 화합물이
겹친다. 우리 서술자 공간이 그 어레이가 실제로 무엇을 묶었는지 예측하는지
검정한다. 새 측정이 아니라 '기존 측정에 대한 검증'이다.

Davidson 정보원. 전문 8면 (ACS Sens. 2020, 5, 1102-1109), Table 1
"Cluster Members" + 본문 "Cluster 6 contains exposures with very weak
color changes and includes all blank exposures. Any compound not listed
as a member of a different cluster in Table 1 is a member of cluster 6."

이원 결과 (Davidson Table 1 에서 직독, 우리 판단 아님):
  진본 작용제와 동일 군집  = DCP(c10: GB,GD,GF) / CEES(c7: HD,HT,T)
                             / malathion(c5: GA)
  약반응·블랭크 군집(c6)   = DMMP(c6 "various; includes ... DMMP")
                             / dichloromethane / methyl salicylate

2026-11-13 단서. 뒤 둘은 Table 1 에 없고 전문 어디에도 이름이 없다. 본문
규칙은 "미등재 = c6" 이지만 그 규칙은 **노출된** 화합물에만 적용된다.
Davidson 은 "The test matrix included 56 analytes" 라 적고 "complete list
of tested analytes can be found in the Supporting Information" 으로 넘긴다.
그 SI 를 얻지 못했다 (pubs.acs.org client_error, 2026-11-13). 따라서
두 화합물이 56종 안에 있었는지는 [수집 불가] 다. ethanol 도 같은 이유로
판정할 수 없어 비교군에 넣지 않았다.

그래서 두 가지를 함께 낸다. primary = 6종(현행), verifiable_only =
Table 1 에서 직접 확인되는 4종(DCP/CEES/malathion/DMMP). 후자는 3대1 이라
도달 가능한 최소 p 가 0.25 이므로 검정력이 없고, 그 사실을 적는다.

사전 지정 검정. 각 후보에서 Davidson 이 측정한 진본 작용제 7종
(GA,GB,GD,GF,VX,HD,L) 까지의 최소거리. 군집 동반 3종이 비동반 3종보다
가까운가. Mann-Whitney U, 3대3 전수 열거 (C(6,3)=20), 단측.
본문과 같은 공간 정의·같은 표준화를 쓴다. 새 규칙을 쓰면 비교가 깨진다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from pesticide_panel import PATH, dmat, spaces_of

# Davidson 과 우리 패널에 공통인 진본 작용제
DAVIDSON_AGENTS = ['tabun (GA)', 'sarin (GB)', 'soman (GD)', 'cyclosarin (GF)',
                   'VX', 'sulfur mustard (HD)', 'lewisite 1 (L)']
# Davidson Table 1 직독. True = 진본 작용제와 동일 군집
CO_CLUSTERED = {'DCP': True, 'CEES': True, 'malathion': True,
                'DMMP': False, 'dichloromethane': False,
                'methyl salicylate': False}
CLUSTER_NOTE = {'DCP': 'c10 with GB,GD,GF', 'CEES': 'c7 with HD,HT,T',
                'malathion': 'c5 with GA', 'DMMP': 'c6 weak/blank',
                'dichloromethane': 'c6 (Table 1 미등재)',
                'methyl salicylate': 'c6 (Table 1 미등재)'}
# Davidson 전문에 이름이 나오는 것만. 나머지 둘은 56종 안에 있었는지
# 확인할 수 없다 (분석물 전체 목록이 SI 에만 있고 SI 를 얻지 못했다).
VERIFIABLE = ['DCP', 'CEES', 'malathion', 'DMMP']
EXPOSURE_UNVERIFIED = ['dichloromethane', 'methyl salicylate']
UNRESOLVED = {'ethanol': ('[수집 불가: Davidson 분석물 전체 목록이 SI 에만 '
                          '있고 SI 를 얻지 못했다. 우리 패널에는 간섭물로 '
                          '있으나 노출 여부를 판정할 수 없다]')}


def mwu_exact(x_yes, x_no):
    """작은 값이 yes 군에 몰릴 확률. 3대3 전수 열거, 단측."""
    allv = list(x_yes) + list(x_no)
    k = len(x_yes)
    r = {v: i for i, v in enumerate(sorted(allv))}
    obs = sum(r[v] for v in x_yes)
    combos = list(itertools.combinations(range(len(allv)), k))
    null = [sum(r[allv[i]] for i in c) for c in combos]
    p = sum(1 for s in null if s <= obs) / len(null)
    return obs, p, len(combos)


def main():
    df = pd.read_csv(PATH)
    sp = spaces_of(df)
    names = df['compound_name'].tolist()
    missing = [c for c in list(CO_CLUSTERED) + DAVIDSON_AGENTS if c not in names]
    if missing:
        raise SystemExit(f'패널에 없음: {missing}')
    ai = [names.index(a) for a in DAVIDSON_AGENTS]

    out = {'source': 'Davidson et al. 2020, ACS Sens. 5, 1102-1109, Table 1',
           'agents_overlapping': DAVIDSON_AGENTS,
           'davidson_outcome': CO_CLUSTERED, 'cluster_note': CLUSTER_NOTE,
           'exposure_verified_from_article': VERIFIABLE,
           'exposure_unverified': EXPOSURE_UNVERIFIED,
           'unresolved': UNRESOLVED,
           'test': 'Mann-Whitney U, exact enumeration, one-sided',
           'spaces': {}, 'verifiable_only': {}}

    for sname, cols in sp.items():
        D = dmat(df, cols)
        d_min = {c: float(D[names.index(c), ai].min()) for c in CO_CLUSTERED}
        nn_med = float(np.median(D.min(1)))
        yes = [d_min[c] for c in CO_CLUSTERED if CO_CLUSTERED[c]]
        no = [d_min[c] for c in CO_CLUSTERED if not CO_CLUSTERED[c]]
        obs, p, n = mwu_exact(yes, no)
        order = sorted(d_min, key=d_min.get)
        perfect = all(CO_CLUSTERED[c] for c in order[:3])
        out['spaces'][sname] = {
            'n_coords': len(cols), 'nn_median': nn_med,
            'd_min_to_davidson_agents': d_min,
            'ranking_closest_first': order,
            'perfect_separation': bool(perfect),
            'rank_sum_coclustered': obs, 'p_one_sided': p, 'n_permutations': n,
            'd_min_normalised': {c: d_min[c] / nn_med for c in d_min}}
        # 노출이 전문에서 확인되는 4종만. 3대1 이므로 최소 도달 p = 0.25.
        vy = [d_min[c] for c in VERIFIABLE if CO_CLUSTERED[c]]
        vn = [d_min[c] for c in VERIFIABLE if not CO_CLUSTERED[c]]
        vobs, vp, vn_ = mwu_exact(vy, vn)
        out['verifiable_only'][sname] = {
            'compounds': VERIFIABLE, 'n_coclustered': len(vy),
            'n_not': len(vn), 'p_one_sided': vp, 'n_permutations': vn_,
            'smallest_attainable_p': round(1.0 / vn_, 4)}
        print(f'[{sname}] n_coords={len(cols)} p={p:.4f} '
              f'perfect={perfect}  order={order}')
        print(f'    검증가능 4종만: p={vp:.4f} '
              f'(도달 가능 최소 {1.0/vn_:.4f}, 검정력 없음)')

    with open('outputs/davidson_external.json', 'w') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print('-> outputs/davidson_external.json')


if __name__ == '__main__':
    main()
