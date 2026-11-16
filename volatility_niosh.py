#!/usr/bin/env python3
"""휘발성 축에 NIOSH 후보 풀을 복원한다 (2026-11-16 신설).

왜 필요한가. niosh_panel.py 가 서술자 공간에 대해 수행한 비교풀 복원
(원 지정문서가 열거한 후보 11종 중 우리 패널에 있는 것은 DMMP 하나뿐이고,
나머지를 되돌려 놓으면 구조공간 0.028 → 0.360, 유효매칭 3 → 0) 이 휘발성
축에는 수행되지 않았다. 심사에서 두 리뷰어가 독립으로 지적했다.

원고가 서술자 공간을 전자/병합으로 확장하지 않은 사유 — xtb 바이너리
버전이 기록되지 않아 두 바이너리 출력을 한 거리공간에 섞게 된다 — 는
휘발성 축에 해당하지 않는다. 이 축은 측정된 1차원 양이고 반경험 계산이
들어가지 않는다.

[C-3]. 증기압을 추정하지 않는다. volatility_fetch.py 와 같은 경로
(PubChem PUG-View, heading=Vapor Pressure) 로 후보를 수집하고,
volatility_space.py 가 사전 선언한 같은 3-tier 규칙으로 채택했다.
수집 결과: 10종 중 5종에 후보가 있고 5종은 0건이다.
"""
import itertools
import json
import os

import numpy as np
import pandas as pd

import volatility_space as V

# ---------------------------------------------------------------------------
# 채택 증기압. (mmHg, 섭씨, tier, 출처 문자열 그대로, 선택 사유)
# 출처 문자열은 PubChem PUG-View 원문이다. 손으로 고치지 말 것.
# ---------------------------------------------------------------------------
ADDED_VP = {
 'DIMP': (0.17, 25, 'T1',
          'Bennett SR et al; Environmental hazards of chemical agent '
          'simylants CRDC-TR-84055. Chemical Research and Development Ctr, '
          'Aberdeen Proving Ground, MD 1984 | Vapor pressure = 0.17 mmHg '
          'at 25 °C',
          'CRDC(ECBC 전신) 기술보고서의 25 °C 값이므로 T1. Haz-Map 과 '
          'DOE PAC 은 0.28 @25 °C 로 더 높으나 T2 이므로 규칙상 뒤다. '
          '낮은 값을 채택하는 쪽이 이 화합물을 DMMP 에 더 가깝게 놓지 '
          '않는다'),
 'TEP': (0.39, 25, 'T1',
         'Daubert, T.E., R.P. Danner. Physical and Thermodynamic Properties '
         'of Pure Chemicals Data Compilation. Washington, D.C.: Taylor and '
         'Francis, 1989. | 0.39 mm Hg at 25 °C',
         'T1 목록에 이름이 있는 편람의 25 °C 값. Haz-Map 과 DOE PAC 도 '
         '0.39 로 일치한다'),
 '2,4-dichlorophenol': (0.09, 25, 'T1',
         'Shiu WY et al; Chemosphere 29: 1155-224 (1994) | '
         '0.09 mm Hg at 25 °C',
         'HSDB 가 인용하는 1차 문헌의 25 °C 값. DOE PAC 도 0.09 로 '
         '일치한다. Haz-Map 만 0.11'),
 # 아래 둘은 Haz-Map 단독이고 출처가 온도를 적지 않는다. T2 는 '25 °C 값'
 # 을 요구하므로 규칙을 엄히 읽으면 [수집 불가] 다. 다만 물을 Haz-Map
 # 23.75 mmHg 로 T2 채택한 선례가 있고, 풀을 키우는 것은 우리 결과에
 # 불리한 방향이므로 민감도에만 넣고 1차 결과에서는 뺀다.
 'di-n-butyl sulfide': (1.2, None, 'T2?',
         'Haz-Map, Information on Hazardous Chemicals and Occupational '
         'Diseases | 1.2 [mmHg]',
         '출처가 온도를 적지 않는다. 민감도 전용'),
 'DEEP': (0.31, None, 'T2?',
          'Haz-Map, Information on Hazardous Chemicals and Occupational '
          'Diseases | 0.31 [mmHg]',
          '출처가 온도를 적지 않는다. 민감도 전용'),
}

# PubChem PUG-View 에 Vapor Pressure 항목이 0건인 것. 추정하지 않고 뺀다.
NOT_OBTAINED = {
 '2-chloroethyl phenyl sulfide': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 (CID 21702)]',
 '2-chloroethyl cyclohexyl sulfide': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 (CID 22469748)]',
 '1,6-dichlorohexane': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 (CID 16551)]',
 'bis(4-chlorobutyl) ether': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 (CID 22800)]',
 'DEMP': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 (CID 12685)]',
}

DATED = ['DIMP', 'TEP', '2,4-dichlorophenol']
UNDATED = ['di-n-butyl sulfide', 'DEEP']


def run(added):
    """28종 패널에 added 를 간섭물로 넣고 같은 전수열거를 돌린다."""
    df = pd.read_csv('descriptor_cache_xtb.csv')
    names = list(df.compound_name) + list(added)
    kl = list(df.compound_class) + ['interferent'] * len(added)
    vp = dict(V.VP)
    for n in added:
        vp[n] = ADDED_VP[n]
    x = np.array([np.log10(vp[n][0]) for n in names], float)
    z = (x - x.mean()) / x.std(ddof=1)
    D = np.abs(z[:, None] - z[None, :])
    res, nperm = V.exact_perm(D, kl, ['simulant'])
    n_nonagent = sum(1 for k in kl if k != 'real_agent')
    return {'added': list(added), 'n_compounds': len(names),
            'n_nonagent': n_nonagent, 'n_permutations': nperm,
            'statistics': res,
            'prespecified': {'statistic': V.PRESPECIFIED,
                             **res[V.PRESPECIFIED]},
            'log10_vp_added': {n: float(np.log10(ADDED_VP[n][0]))
                               for n in added}}


def main():
    np.random.seed(V.SEED)
    print('=' * 78)
    print('[V-NIOSH] 휘발성 축 비교풀 복원')
    print('=' * 78)
    out = {'generated_by': 'volatility_niosh.py',
           'source_document': 'NIOSH CBRN APR simulant programme (31 July 2002)',
           'adoption_rule': V.VP and {
               'T1': '25 °C 증기압 전용 측정·편람값',
               'T2': 'T1 없을 때 정부 편람 25 °C 값',
               'T3': '25 °C 값이 없으면 최근접 온도 + flag'},
           'adopted': {k: {'mmHg': v[0], 'temp_C': v[1], 'tier': v[2],
                           'source': v[3], 'reason': v[4]}
                       for k, v in ADDED_VP.items()},
           'not_obtained': NOT_OBTAINED,
           'n_candidates_in_source': 11,
           'n_on_panel_already': 1,
           'runs': {}}
    for label, added in (('baseline', []),
                         ('plus3_dated', DATED),
                         ('plus5_with_undated', DATED + UNDATED)):
        r = run(added)
        out['runs'][label] = r
        p = r['prespecified']
        print(f"  {label:<20} n={r['n_compounds']:>3} "
              f"비작용제={r['n_nonagent']:>3} 열거={r['n_permutations']:>5} "
              f"p={p['p']:.5f}  n_valid={r['statistics']['n_valid']['observed']:.0f}")
    os.makedirs('outputs', exist_ok=True)
    json.dump(out, open('outputs/volatility_niosh.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n-> outputs/volatility_niosh.json')


if __name__ == '__main__':
    main()
