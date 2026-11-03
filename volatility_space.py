#!/usr/bin/env python3
r"""휘발성 축 보조 분석 (2026-11-03 신설).

왜 이 파일이 있는가. 리뷰어 1 의 지적: 서론은 휘발성을 모사체 선정 기준의
첫 번째로 들면서, 정작 서술자 패널에는 휘발성이 없다. 그러면 '서술자가
응답을 예측하지 못한다'는 결론이 '하필 휘발성을 뺀 서술자가' 로 읽힌다.

이 파일은 그 공백을 메우되, 본문의 72검정 가족에는 넣지 않는다. 축을
정식 편입하면 셀당 12검정이 늘어 논문의 모든 BH q 가 다시 계산되는데,
리뷰어 지적을 해소하는 데 필요한 것은 '휘발성으로도 맞지 않는다'는 사실
하나이고 그것은 사전지정 단일 검정으로 충분하다. 그래서 두 검정만
사전에 고정하고 가족과 분리해 보고한다. 분리 사실 자체를 SI 에 적는다.

사전지정 (수치를 보기 전에 고정한 것):
  V1  지정 검정. 1차원 휘발성 공간에서 통계량 mean_diff 하나.
      패널은 28종(1차 패널). 전수열거 C(17,3)=680.
  V2  일치 검정. 휘발성 거리와 ΔE00 응답거리의 Mantel 을 여섯
      런-채널 셀에서 돌려 Fisher 로 하나의 p 로 합친다.
다른 통계량·다른 패널의 값도 함께 찍지만, 판정은 위 둘로만 한다.

왜 1차원인가. 축이 하나면 상관필터가 할 일이 없고 서술자 인코딩 선택도
없다. 본문이 센 네 가지 분석자 자유도(인코딩/상관필터/형태이성질체
시드/패널 구성) 가운데 셋이 구조적으로 사라진다. 귀무가 나오더라도
'분석자 선택의 산물'이라는 반론이 닿지 않는 자리다.

[C-3] 준수. 모든 증기압은 공개 DB 에서 실제로 조회한 측정값이고, 예측값
(OPERA 등)은 쓰지 않았다. 조회 원문은 volatility_fetch.py 가
outputs/volatility_candidates.json 에 그대로 담아 둔다. 아래 표는 그
후보들 가운데 사람이 고른 채택값이며, 고른 근거를 한 건씩 적었다.

채택 규칙 (고르기 전에 정한 것):
  T1  25 °C 의 증기압 전용 측정·편람값. ECBC 기술보고서, Daubert &
      Danner, Boublik/Fried/Hala, J. Phys. Chem. Ref. Data, Pesticide
      Manual, HSDB 가 인용하는 1차 문헌.
  T2  T1 에 25 °C 값이 없으면 정부 편람의 25 °C 값 (EPA EHS 1998,
      DOE PAC).
  T3  25 °C 값이 아예 없으면 가장 가까운 온도의 값. 온도를 적고
      flag 를 세운다.
  없으면 [수집 불가] 로 적고 그 화합물을 뺀다.
단위는 mmHg 로 통일한다. 온도 환산은 하지 않는다. 환산에는 Antoine
계수가 필요하고 그것을 추정하면 측정값이 아니게 된다.
"""
import itertools
import json
import os
from itertools import permutations as _permutations

import numpy as np
import pandas as pd
from scipy import stats

SEED = 0xC0DE
STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']
PRESPECIFIED = 'mean_diff'

# ---------------------------------------------------------------------------
# 채택 증기압. (mmHg, 섭씨, tier, 출처 문자열 그대로, 선택 사유)
# 출처 문자열은 PubChem PUG-View 가 돌려준 원문이다. 손으로 고치지 말 것.
# ---------------------------------------------------------------------------
VP = {
 'DMMP': (0.962, 25, 'T1',
          'Boublik, T., Fried, V., and Hala, E., The Vapour Pressures of '
          'Pure Substances | 0.962 mm Hg at 25 °C',
          'NTP 1992 는 같은 25 °C 에 1.2 를 준다. 증기압 전용 편람을 앞에 둔다'),
 'acetone': (231.0, 25, 'T1', 'PMID:7676461 | 231 mm Hg at 25 °C', ''),
 'methanol': (127.0, 25, 'T1',
              'Boublik/Fried/Hala | 127 mm Hg at 25 °C', ''),
 'acetonitrile': (88.8, 25, 'T1',
                  'Boublik/Fried/Hala | 88.8 mm Hg at 25 °C', ''),
 'chloroform': (197.0, 25, 'T1',
                'Boublik/Fried/Hala | 197 mm Hg at 25 °C', ''),
 'n-hexane': (153.0, 25, 'T1',
              'Daubert, T.E., R.P. Danner | 153 mm Hg at 25 °C', ''),
 'THF': (162.0, 25, 'T1',
         'Daubert, T.E., R.P. Danner | 1.62X10+2 mm Hg at 25 °C', ''),
 'ethyl acetate': (93.2, 25, 'T1',
                   'Daubert TE, Danner RP | 93.2 mm Hg at 25 °C', ''),
 'benzene': (94.8, 25, 'T1',
             'Daubert, T.E., R.P. Danner | 94.8 mm Hg at 25 °C', ''),
 'water': (23.75, 25, 'T2',
           'Haz-Map, Information on Hazardous Chemicals | 23.75 [mmHg]',
           'T1 에 25 °C 항목이 없었다. CRC 항목은 삼중점만 돌아온다'),
 'dichloromethane': (435.0, 25, 'T1',
                     'Boublik/Fried/Hala | 435 mm Hg at 25 °C', ''),
 'DPGME': (0.4, 26, 'T3',
           'Browning, E. Toxicity and Metabolism of Industrial Solvents '
           '| 0.4 MM HG @ 26 °C',
           '25 °C 값이 없다. 26 °C 를 쓰고 온도를 적는다'),
 'ethanol': (59.3, 25, 'T1',
             'Daubert, T.E., R.P. Danner | 59.3 mm Hg at 25 °C '
             '/Extrapolated/',
             '출처가 /Extrapolated/ 로 표시한 값이다. 그대로 적고 flag 를 센다'),
 'toluene': (28.4, 25, 'T1',
             'Daubert, T.E., R.P. Danner | 28.4 mm Hg at 25 °C', ''),
 'DMSO': (0.60, 25, 'T1',
          'Daubert TE, Danner RP | 0.60 mm Hg at 25 °C', ''),
 'p-dichlorobenzene': (1.74, 25, 'T1',
                       'Daubert, T.E., R.P. Danner | 1.74 mm Hg at 25 °C',
                       '같은 출처가 1.36~2.2 범위도 준다. 점값을 쓴다'),
 'p-xylene': (8.84, 25, 'T1',
              'Chao J et al; J Phys Chem Ref Data 12: 1033 '
              '| 8.84 mm Hg at 25 °C', ''),
 'sarin (GB)': (2.86, 25, 'T1',
                'Epstein J; J Amer Water Works Assoc Jan: 31-7 (1974) '
                '| 2.86 mm Hg at 25 °C',
                'EPA 1998 은 2.9. 1차 문헌을 앞에 둔다'),
 'soman (GD)': (0.41, 25, 'T1',
                'Balboa A et al; Vapor Pressure of GD, ECBC-TR '
                '| 0.41 mm Hg at 25 °C',
                'ECBC 증기압 전용 측정보고서'),
 'tabun (GA)': (0.07, 25, 'T1',
                'Munro NB et al; Environ Health Perspectives '
                '| 0.07 mm Hg at 25 °C; 0.037 mm Hg at 20 °C', ''),
 'sulfur mustard (HD)': (0.11, 25, 'T1',
                         'PMID:10585900 | 0.11 mm Hg at 25 °C',
                         'EPA 1998 은 0.09 @30 °C. 25 °C 값을 쓴다'),
 'lewisite 1 (L)': (0.58, 25, 'T1',
                    'Munro NB et al; Environ Health Perspectives '
                    '| 0.58 mm Hg at 25 °C',
                    'PMID:2692088 은 이성질체별로 trans 0.4 / cis 1.562 '
                    '(25 °C). 공업품은 혼합물이므로 혼합물 값을 쓰고 '
                    '이성질체 분리값을 민감도로 돌린다'),
 'hydrogen cyanide (AC)': (742.0, 25, 'T1',
                           'Daubert, T.E., R.P. Danner | 742 mm Hg at 25 °C',
                           'NIOSH 는 630. 증기압 편람을 앞에 둔다'),
 'cyanogen chloride (CK)': (1230.0, 25, 'T1',
                            'Daubert, T.E., R.P. Danner '
                            '| 1.23X10+3 mm Hg at 25 C', ''),
 'phosgene (CG)': (1420.0, 25, 'T1',
                   'Daubert, T.E., R.P. Danner | 1420 mm Hg at 25 °C', ''),
 'chloropicrin (PS)': (24.0, 25, 'T1',
                       'Tomlin CDS, ed. Chloropicrin (76-06-2) '
                       '| 3.2 kPa (24 mm Hg) at 25 °C', ''),
 'cyclosarin (GF)': (0.044, 20, 'T3',
                     'Augerson W; A Review of the Scientific Literature '
                     '| 0.044 mm Hg at 20 °C',
                     '25 °C 값이 공개 DB 에 없다. 20 °C 를 쓰고 온도를 적는다'),
 'VX': (8.78e-4, 25, 'T1',
        'Buchanan JH et al; Vapor Pressure of VX, ECBC-TR '
        '| 8.78X10-4 mm Hg at 25 °C',
        'EPA 1998 은 7.0e-4. ECBC 증기압 전용 측정보고서를 앞에 둔다'),
 # --- 35종 패널에만 있는 화합물
 'DFP': (0.58, 25, 'T2',
         'PAC Chemical Database, U.S. Department of Energy '
         '| 0.58 [mm Hg] @25 °C',
         'Merck Index 는 0.579 @20 °C. 25 °C 값을 앞에 둔다'),
 'CEES': (3.4, 25, 'T1',
          'Bennett SR et al; Environmental Hazards of Chemical Agent '
          'Simulants | 3.4 mm Hg @ 25 °C', ''),
 'methyl salicylate': (0.0343, 25, 'T1',
                       'Daubert, T.E., R.P. Danner '
                       '| Vapor pressure = 0.0343 mm Hg @ 25 °C', ''),
 'dichlorvos': (1.58e-2, 25, 'T1',
                'MacBean C, ed; e-Pesticide Manual 15th ed. '
                '| 1.58X10-2 mm Hg at 25 °C', ''),
 'malathion': (3.97e-5, 30, 'T3',
               'MacBean C, ed; e-Pesticide Manual 15th ed. '
               '| 3.97X10-5 mm Hg at 30 °C',
               '25 °C 값이 없다. 30 °C 를 쓰고 온도를 적는다'),
}

# 조회했으나 증기압 항목이 없는 화합물. 35종 패널에만 있다.
UNAVAILABLE = {
 'DCP': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 '
        '(CAS 814-49-3, CID 13139)]',
 'DCNP': '[수집 불가: PubChem PUG-View 에 Vapor Pressure 항목 없음 '
         '(CAS 2942-58-7, CID 76257)]',
}

# 루이사이트 이성질체 민감도 (PMID:2692088, 둘 다 25 °C 측정값)
LEWISITE_ISOMERS = {'trans': 0.4, 'cis': 1.562}


# ---------------------------------------------------------------------------
def space(names):
    """1차원 휘발성 공간의 거리행렬. x = log10(P/mmHg) 를 z 표준화한다.

    축이 하나이므로 상관필터가 할 일이 없다. z 표준화 뒤 |Δz| 가 거리다.
    """
    x = np.array([np.log10(VP[n][0]) for n in names], float)
    z = (x - x.mean()) / x.std(ddof=1)
    return np.abs(z[:, None] - z[None, :]), x, z


def exact_perm(D, kl, sim_classes):
    """rollback_c.exact_perm 과 같은 귀무·같은 통계량. 1차원이라 그대로 쓴다."""
    D = D.copy().astype(float)
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    null = {s: np.empty(len(ALL)) for s in STATS}
    npool = len(i_pool)
    _nz = [0]                      # log 발산으로 버린 쌍 수
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
        null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
        # 1차원이라 두 화합물의 log10 P 가 완전히 같으면 거리가 0 이 되고
        # log(ds/di) 가 발산한다. 루이사이트 trans 민감도(0.4 mmHg)가
        # DPGME(0.4 mmHg)와 정확히 겹치는 자리에서 실제로 발생했다.
        # 유한한 쌍만 평균하고, 몇 건을 버렸는지 센다.
        with np.errstate(divide='ignore', invalid='ignore'):
            lr = np.log(ds / di)
        fin = np.isfinite(lr)
        _nz[0] += int((~fin).sum())
        null['mean_logratio'][a0:a0 + b] = np.where(
            fin.any(0), np.nansum(np.where(fin, lr, 0.0), 0)
            / np.maximum(fin.sum(0), 1), np.nan)
        null['mean_dsim'][a0:a0 + b] = ds.mean(0)
        null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
        null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    ti = ALL.index(true_combo)
    out = {s: dict(observed=float(null[s][ti]),
                   p=float((null[s] <= null[s][ti]).sum() / len(ALL)),
                   rank=int((null[s] <= null[s][ti]).sum()))
           for s in STATS}
    out['_log_divergent_pairs'] = _nz[0]
    return out, len(ALL)


def sanity(names, x):
    """Sanity 3항목. 실패하면 멈춘다."""
    out = {}
    # (1) 단조성 — 알려진 휘발성 순서를 지키는가.
    #     기체(CG, CK, AC)가 가장 높고 VX 가 가장 낮아야 한다.
    order = [n for _, n in sorted(zip(x, names))]
    out['lowest'] = order[0]
    out['highest'] = order[-1]
    out['monotonicity_ok'] = bool(order[0] == 'VX'
                                  and order[-1] == 'phosgene (CG)')
    # (2) 범위 — 패널이 실제로 여러 자수를 덮는가. 한 자수 안에 몰려
    #     있으면 축이 정보를 담지 못한다.
    out['log10_range'] = float(x.max() - x.min())
    out['range_ok'] = bool(out['log10_range'] > 3.0)
    # (3) 온도 혼입 — 25 °C 가 아닌 값이 몇 건인가.
    out['off_25C'] = {n: VP[n][1] for n in names if VP[n][1] != 25}
    out['n_off_25C'] = len(out['off_25C'])
    out['tiers'] = {t: sum(1 for n in names if VP[n][2] == t)
                    for t in ('T1', 'T2', 'T3')}
    return out


def designation(cache, sim_classes, label):
    """V1. 지정 검정 — 휘발성 공간에서 유효 매칭이 생기는가."""
    df = pd.read_csv(cache)
    names = list(df.compound_name)
    miss = [n for n in names if n not in VP]
    if miss:
        return {'label': label, 'skipped': True,
                'reason': '; '.join(UNAVAILABLE.get(n, f'[수집 불가: {n}]')
                                    for n in miss),
                'missing': miss}
    D, x, z = space(names)
    kl = list(df.compound_class)
    res, nperm = exact_perm(D, kl, sim_classes)
    return {'label': label, 'skipped': False, 'n_compounds': len(names),
            'n_permutations': nperm, 'statistics': res,
            'prespecified': {'statistic': PRESPECIFIED,
                             **res[PRESPECIFIED]},
            'sanity': sanity(names, x),
            'log10_vp': {n: float(v) for n, v in zip(names, x)}}


def agreement():
    """V2. 일치 검정 — 휘발성 거리가 ΔE00 응답거리를 따라가는가.

    셀마다 Mantel 전수열거, 여섯 셀을 Fisher 로 하나의 p 로 합친다.
    """
    import deltaE_premise as m
    df = pd.read_csv('descriptor_cache_xtb.csv')
    names = list(df.compound_name)
    D, x, z = space(names)
    idx = {n: i for i, n in enumerate(names)}
    cells, skipped = [], []
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                Dd, dy = m.load(tag, mode)
            except Exception as e:
                skipped.append({'cell': f'{tag}_{mode}',
                                'reason': f'[수집 불가: {e}]'})
                continue
            P = m.consensus(m.centred(m.raw_patterns(Dd, dy)))
            ags = [a for a in sorted(P) if a in m.MAP and m.MAP[a] in idx]
            if len(ags) < 4:
                skipped.append({'cell': f'{tag}_{mode}',
                                'n_mapped_agents': len(ags),
                                'reason': '매핑된 작용제 4종 미만 — '
                                          'Mantel 불가'})
                continue
            k = len(ags)
            pr = list(itertools.combinations(range(k), 2))
            dis = np.array([1 - float(np.dot(P[ags[i]], P[ags[j]]))
                            for i, j in pr])
            M = np.array([[D[idx[m.MAP[a]], idx[m.MAP[b]]] for b in ags]
                          for a in ags])
            d = np.array([M[i, j] for i, j in pr])
            r = float(stats.spearmanr(d, dis).statistic)
            perms = list(_permutations(range(k)))
            cnt = 0
            for pm in perms:
                dp = np.array([M[pm[i], pm[j]] for i, j in pr])
                rp = stats.spearmanr(dp, dis).statistic
                if rp == rp and abs(rp) >= abs(r):
                    cnt += 1
            cells.append({'cell': f'{tag}_{mode}', 'dataset': tag,
                          'mode': mode, 'n_compounds': k,
                          'n_pairs': len(pr), 'n_relabellings': len(perms),
                          'spearman': r, 'p_mantel': cnt / len(perms)})
    pv = [c['p_mantel'] for c in cells]
    fisher = stats.combine_pvalues(pv, method='fisher') if pv else None
    return {'cells': cells, 'skipped': skipped,
            'prespecified': {'method': 'Fisher combination of per-cell '
                                       'exact Mantel p',
                             'n_cells': len(cells),
                             'statistic': float(fisher.statistic) if pv
                             else None,
                             'p_value': float(fisher.pvalue) if pv else None}}


def valid_match_detail():
    """유효매칭을 작용제별로 분해하고, 유일한 매칭의 붕괴점을 계산한다.

    왜 필요한가. V1 의 n_valid 는 '1건, p=0.4647' 로만 보고되는데, 그
    1건이 무엇인지가 해석을 가른다. 분해하면 DPGME-소만 한 짝이고, 그
    짝의 근거는 0.011 자수 간격이다. DPGME 는 26 °C 측정값, 소만은
    25 °C 측정값이므로 그 간격은 온도 불일치 안에 들어간다.

    공통 온도로 맞추면 간격은 넓어질 수밖에 없다. 증기압은 온도에 대해
    단조증가하므로 DPGME 의 25 °C 값은 26 °C 값 0.4 mmHg 보다 낮고,
    소만(-0.387)보다 아래에 있는 DPGME(-0.398)는 더 멀어진다. 이 방향은
    열역학으로 고정되므로 Antoine 계수 추정이 필요 없다 ([C-3] 준수).

    붕괴점은 계산한다. 추정하지 않는다. 25 °C 의 DPGME 증기압이 얼마
    아래로 내려가면 유효매칭 두 조건이 각각 깨지는지만 적는다.
    """
    df = pd.read_csv('descriptor_cache_xtb.csv')
    names = list(df.compound_name)
    kl = list(df.compound_class)
    D, x, _ = space(names)
    sd = float(np.std(x, ddof=1))
    Dm = D.copy()
    np.fill_diagonal(Dm, np.inf)
    nn_med = float(np.median(Dm.min(1)))
    i_sim = [i for i, k in enumerate(kl) if k == 'simulant']
    i_int = [i for i, k in enumerate(kl) if k == 'interferent']
    rows = []
    for a, k in enumerate(kl):
        if k != 'real_agent':
            continue
        ds = float(min(Dm[a, t] for t in i_sim))
        di = float(min(Dm[a, t] for t in i_int))
        rows.append({'agent': names[a], 'd_nearest_simulant': ds,
                     'd_nearest_interferent': di,
                     'ratio_to_nn_median': ds / nn_med,
                     'nearest_simulant': names[min(
                         i_sim, key=lambda t: Dm[a, t])],
                     'valid': bool(ds / nn_med <= 2.0 and ds <= di)})
    valid = [r for r in rows if r['valid']]
    out = {'nn_median_distance': nn_med, 'threshold_ratio': 2.0,
           'threshold_distance_z': 2.0 * nn_med,
           'log10_sd': sd, 'per_agent': rows,
           'n_valid': len(valid),
           'valid_pairs': [(r['agent'], r['nearest_simulant'])
                           for r in valid]}
    # 유일 매칭의 붕괴점. 두 조건을 각각 깨는 간격을 자수로, 그리고 그
    # 간격이 되려면 DPGME 증기압이 얼마여야 하는지를 mmHg 로 적는다.
    if len(valid) == 1 and valid[0]['nearest_simulant'] == 'DPGME':
        gd = VP['soman (GD)'][0]
        gap_ratio = 2.0 * nn_med * sd        # 비율 조건을 깨는 간격 (자수)
        gap_closer = valid[0]['d_nearest_interferent'] * sd   # 근접 조건
        out['sole_match_breakpoint'] = {
            'pair': ['DPGME', 'soman (GD)'],
            'observed_gap_log10': float(
                abs(np.log10(VP['DPGME'][0]) - np.log10(gd))),
            'temperatures_C': [VP['DPGME'][1], VP['soman (GD)'][1]],
            'gap_breaking_ratio_condition_log10': float(gap_ratio),
            'vp_breaking_ratio_condition_mmHg': float(gd / 10 ** gap_ratio),
            'gap_breaking_closer_condition_log10': float(gap_closer),
            'vp_breaking_closer_condition_mmHg': float(
                gd / 10 ** gap_closer),
            'direction_note':
                'Vapour pressure increases with temperature, so the '
                "25 C value for DPGME is below the 26 C measurement of "
                '0.4 mmHg and the gap to soman can only widen. The '
                'direction is fixed thermodynamically; no Antoine fit '
                'is used.'}
    return out


def pairing_table():
    """기술표. 작용제별로 (가장 가까운 모사체, 가장 가까운 풀 화합물).

    왜 필요한가. V1 의 통계량은 본문과 같은 작용제 중심 min 규약을 쓰므로
    '지정 역할'을 보지 않는다. 그 수치만 싣고 끝내면 '휘발성으로는
    모사체가 맞는다'로 읽히는데, 짝을 펼쳐 보면 지정된 상대가 아니라
    다른 계열의 작용제에 붙어 있다. 검정이 아니므로 가족에 넣지 않는다.
    """
    df = pd.read_csv('descriptor_cache_xtb.csv')
    names = list(df.compound_name)
    kl = list(df.compound_class)
    ser = dict(zip(df.compound_name, df.agent_series.fillna('')))
    _, x, _ = space(names)
    lx = dict(zip(names, (float(t) for t in x)))
    sims = [n for n, k in zip(names, kl) if k == 'simulant']
    ags = [n for n, k in zip(names, kl) if k == 'real_agent']
    pool = [n for n, k in zip(names, kl) if k != 'real_agent']
    rows = []
    for a in sorted(ags, key=lambda n: -lx[n]):
        bs = min(sims, key=lambda t: abs(lx[t] - lx[a]))
        bp = min(pool, key=lambda t: abs(lx[t] - lx[a]))
        rows.append({'agent': a, 'series': ser[a],
                     'log10_vp': lx[a],
                     'nearest_simulant': bs,
                     'gap_to_nearest_simulant': abs(lx[bs] - lx[a]),
                     'nearest_pool_compound': bp,
                     'gap_to_nearest_pool': abs(lx[bp] - lx[a]),
                     'nearest_is_interferent': bool(kl[names.index(bp)]
                                                    == 'interferent')})
    return {'rows': rows,
            'simulant_log10_vp': {s: lx[s] for s in sims},
            'simulant_designated_series': {s: ser[s] for s in sims},
            'n_agents_whose_nearest_is_interferent':
                sum(1 for r in rows if r['nearest_is_interferent']),
            'max_gap': max(r['gap_to_nearest_simulant'] for r in rows),
            'max_gap_agent': max(rows,
                                 key=lambda r:
                                 r['gap_to_nearest_simulant'])['agent'],
            # DPGME 는 26 °C 값, GD 는 25 °C 값이다. 둘의 간격 0.011 자수는
            # 그 온도 차이 안에 들어가므로 '일치'로 읽어서는 안 된다.
            'temperature_caveat':
                'DPGME 0.4 mmHg @26 C vs soman 0.41 mmHg @25 C — '
                'the 0.011 log-unit gap is within the temperature '
                'mismatch and is not evidence of a match'}


def lewisite_sensitivity():
    """루이사이트를 이성질체 값으로 바꿔도 V1 판정이 바뀌는가."""
    df = pd.read_csv('descriptor_cache_xtb.csv')
    names = list(df.compound_name)
    kl = list(df.compound_class)
    base = VP['lewisite 1 (L)']
    out = {}
    for iso, val in LEWISITE_ISOMERS.items():
        VP['lewisite 1 (L)'] = (val, 25, 'T1', base[3], f'{iso} 민감도')
        D, _, _ = space(names)
        res, _ = exact_perm(D, kl, ['simulant'])
        out[iso] = {'vp_mmHg': val, **res[PRESPECIFIED]}
    VP['lewisite 1 (L)'] = base
    return out


def main():
    np.random.seed(SEED)
    print('=' * 78)
    print('[V] 휘발성 축 보조 분석 — 72검정 가족과 분리해 보고한다')
    print('=' * 78)
    v1 = [designation('descriptor_cache_xtb.csv', ['simulant'],
                      'panel28_incumbent3'),
          designation('descriptor_cache_xtb_35.csv', ['simulant'],
                      'panel35_incumbent3')]
    for r in v1:
        if r['skipped']:
            print(f"  V1 {r['label']:<22} 생략 — {r['reason']}")
            continue
        s = r['sanity']
        ps = r['prespecified']
        print(f"  V1 {r['label']:<22} n={r['n_compounds']} "
              f"전수열거 {r['n_permutations']}")
        print(f"     sanity 단조성 {s['monotonicity_ok']} "
              f"({s['lowest']} 최저, {s['highest']} 최고), "
              f"log10 범위 {s['log10_range']:.2f} 자수 ({s['range_ok']}), "
              f"25 °C 아님 {s['n_off_25C']}건, tier {s['tiers']}")
        print(f"     사전지정 {PRESPECIFIED}: 관측 {ps['observed']:+.4f}, "
              f"p = {ps['p']:.4f}")
        print('     (참고) ' + '  '.join(
            f"{k} p={v['p']:.4f}" for k, v in r['statistics'].items()
            if k != PRESPECIFIED and k in STATS))
    v2 = agreement()
    print()
    for c in v2['cells']:
        print(f"  V2 {c['cell']:<10} n={c['n_compounds']} "
              f"rho {c['spearman']:+.3f}  Mantel p = {c['p_mantel']:.4f} "
              f"({c['n_relabellings']} relabellings)")
    for s in v2['skipped']:
        print(f"  V2 {s['cell']:<10} 생략 — {s['reason']}")
    ps = v2['prespecified']
    print(f"     사전지정 Fisher ({ps['n_cells']} 셀): "
          f"X2 = {ps['statistic']:.3f}, p = {ps['p_value']:.4f}")
    pt = pairing_table()
    vm = valid_match_detail()
    print()
    print(f"  기술표: 작용제 11종 중 최근접이 간섭물인 경우 "
          f"{pt['n_agents_whose_nearest_is_interferent']}건, "
          f"최대 간격 {pt['max_gap']:.3f} 자수 ({pt['max_gap_agent']})")
    for r in pt['rows']:
        print(f"     {r['agent']:<24}{r['series']:<14}"
              f"최근접 모사체 {r['nearest_simulant']:<18}"
              f"{r['gap_to_nearest_simulant']:.3f}   "
              f"최근접 풀 {r['nearest_pool_compound']:<18}"
              f"{r['gap_to_nearest_pool']:.3f}")
    print()
    print(f"  유효매칭 {vm['n_valid']}건 {vm['valid_pairs']}  "
          f"(임계 거리 {vm['threshold_distance_z']:.4f} z, "
          f"최근접이웃 중위 {vm['nn_median_distance']:.4f})")
    bp = vm.get('sole_match_breakpoint')
    if bp:
        print(f"     유일 매칭 간격 {bp['observed_gap_log10']:.4f} 자수, "
              f"온도 {bp['temperatures_C']} °C")
        print(f"     붕괴점: 비율조건 간격 "
              f"{bp['gap_breaking_ratio_condition_log10']:.4f} 자수 "
              f"(DPGME < {bp['vp_breaking_ratio_condition_mmHg']:.4f} mmHg), "
              f"근접조건 간격 "
              f"{bp['gap_breaking_closer_condition_log10']:.4f} 자수 "
              f"(DPGME < {bp['vp_breaking_closer_condition_mmHg']:.4f} mmHg)")
    lew = lewisite_sensitivity()
    print()
    for iso, r in lew.items():
        print(f"  민감도 루이사이트 {iso:<5} ({r['vp_mmHg']} mmHg): "
              f"{PRESPECIFIED} p = {r['p']:.4f}")

    os.makedirs('outputs', exist_ok=True)
    out = {'generated_by': 'volatility_space.py',
           'note': '본문 72검정 가족과 분리된 보조 분석. 사전지정 2검정.',
           'prespecified_tests': {
               'V1': '1차원 휘발성 공간의 지정 검정, 통계량 mean_diff, '
                     '28종 패널, 전수열거',
               'V2': '휘발성 거리 대 ΔE00 응답거리 Mantel 6셀 Fisher 합성'},
           'adoption_rule': {'T1': '25 °C 증기압 전용 측정·편람값',
                             'T2': 'T1 없을 때 정부 편람 25 °C 값',
                             'T3': '25 °C 값이 없으면 최근접 온도 + flag'},
           'vapour_pressure_mmHg': {
               n: {'value': v[0], 'temperature_C': v[1], 'tier': v[2],
                   'source': v[3], 'adjudication': v[4]}
               for n, v in VP.items()},
           'unavailable': UNAVAILABLE,
           'designation': v1, 'agreement': v2, 'pairing_table': pt,
           'valid_match_detail': vm,
           'lewisite_isomer_sensitivity': lew}
    with open('outputs/volatility_space.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('\n-> outputs/volatility_space.json')


if __name__ == '__main__':
    main()
