#!/usr/bin/env python3
"""루이사이트 순도 67% 민감도 (2026-11-02 신설, T2-6).

왜 필요한가. 2차 런의 원논문은 순도를 이렇게 적는다. "The purity of each
agent was confirmed by GC-MS analysis; AC, CK, CG, PS, HD, and HN exhibited
purities > 95%, while L was used with a purity of 67% (determined by
quantitative analysis)." 원고는 이 수치를 인용하면서 그것이 결과에 미치는
영향은 다루지 않았다. 67% 순도의 시료에서 나온 응답은 3분의 1이 무엇인지
모르는 물질의 응답이며, 루이사이트는 작용제 11종 중 하나로 매칭 분석과
2차 런 전제 검정에 모두 들어간다.

이 스크립트는 루이사이트를 빼고 같은 분석을 다시 돌려 결론이 바뀌는지 센다.
빼는 것이 옳다고 주장하지 않는다. 두 수를 나란히 두는 것이 목적이다.

1차 런과 3차 런의 순도는 원논문에 없다. 세 PDF 의 Methods 를 직접 검색해
확인했고 ('purit' 0건), 추정하지 않고 [수집 불가] 로 적는다 ([C-3]).

사용: python3 purity_sensitivity.py
"""
import itertools
import json

import numpy as np
import pandas as pd

import deltaE_premise as m
from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, pick

LEWISITE_DESCRIPTOR = 'lewisite 1 (L)'
LEWISITE_RESPONSE = 'L'
THRESHOLD = 2.0
CASES = [
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
]

PURITY_SOURCE = {
    'second_run': {
        'quote': 'The purity of each agent was confirmed by GC-MS analysis; '
                 'AC, CK, CG, PS, HD, and HN exhibited purities > 95%, while '
                 'L was used with a purity of 67% (determined by quantitative '
                 'analysis).',
        'source': 'multiplexref, Methods',
    },
    'first_run': '[수집 불가: 원논문에 순도 서술이 없다 (Methods 전문에서 '
                 "'purit' 0건).]",
    'third_run': '[수집 불가: 원논문에 순도 서술이 없다 (Methods 전문에서 '
                 "'purit' 0건). 반복측정은 triplicate 로 명시되어 있다.]",
}


def merged_space(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return [c for c in num if c != OMEGA_DERIVED]


def match_counts(df, sim_cls):
    keep, _ = pick(df, merged_space(df), PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    kl, name = df.compound_class.tolist(), df.compound_name.tolist()
    i_agent = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    i_cmp = [i for i, k in enumerate(kl)
             if k != 'real_agent' and k not in sim_cls]
    matched = []
    for i in i_agent:
        js = min(i_sim, key=lambda j: D[i, j])
        if D[i, js] / nn_med <= THRESHOLD and all(D[i, js] < D[i, j]
                                                  for j in i_cmp):
            matched.append(name[i])
    return dict(n_dims=len(keep), n_agents=len(i_agent),
                nn_median=round(nn_med, 4), n_match=len(matched),
                matched=sorted(matched))


def premise_cell(tag, mode, drop=None):
    D, dy = m.load(tag, mode)
    if drop:
        D = D[D.Agent != drop]
    cen = m.centred(m.raw_patterns(D, dy))
    if not cen:
        return None
    P = m.consensus(cen)
    names = sorted(cen)
    d = []
    for a in names:
        ks = sorted(cen[a])
        wa = [float(np.dot(cen[a][i], cen[a][j]))
              for i, j in itertools.combinations(ks, 2)]
        ba = [float(np.dot(P[a], P[o])) for o in names if o != a]
        if wa and ba:
            d.append(np.mean(wa) - np.mean(ba))
    d = np.array(d); n = len(d)
    if n < 2:
        return None
    obs = float(d.mean())
    cnt = sum(1 for s in itertools.product([1, -1], repeat=n)
              if float((np.array(s) * d).mean()) >= obs - 1e-15)
    return dict(n_compounds=n, mean_diff=obs, p_exact=cnt / 2 ** n,
                p_min_attainable=1 / 2 ** n)


def main():
    print('=' * 78)
    print('[T2-6] 루이사이트 순도 67% 민감도')
    print('=' * 78)
    rep = {'purity_source': PURITY_SOURCE, 'match': {}, 'premise': {}}

    print('\n(1) 서술자 매칭 — 루이사이트 포함 대 제외')
    for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        with_l = match_counts(df, sim_cls)
        d2 = df[df.compound_name != LEWISITE_DESCRIPTOR].reset_index(drop=True)
        without = match_counts(d2, sim_cls)
        rep['match'][label] = {'with_lewisite': with_l,
                               'without_lewisite': without}
        print(f'  {label:<30} 포함 작용제 {with_l["n_agents"]:>2} 차원 '
              f'{with_l["n_dims"]:>2} 매칭 {with_l["n_match"]}   |   '
              f'제외 작용제 {without["n_agents"]:>2} 차원 '
              f'{without["n_dims"]:>2} 매칭 {without["n_match"]}')
        if with_l['matched'] != without['matched']:
            print(f'      매칭 명단 변화: {with_l["matched"]} -> '
                  f'{without["matched"]}')

    print('\n(2) 2차 런 전제 검정 — 루이사이트 포함 대 제외')
    for mode in ['uv', 'normal']:
        a = premise_cell('d2', mode)
        b = premise_cell('d2', mode, drop=LEWISITE_RESPONSE)
        rep['premise'][f'd2.{mode}'] = {'with_lewisite': a,
                                        'without_lewisite': b}
        print(f'  d2.{mode:<7} 포함 n={a["n_compounds"]} 차 {a["mean_diff"]:+.3f} '
              f'정확p {a["p_exact"]:.4f} (최소 {a["p_min_attainable"]:.4f})   |   '
              f'제외 n={b["n_compounds"]} 차 {b["mean_diff"]:+.3f} '
              f'정확p {b["p_exact"]:.4f} (최소 {b["p_min_attainable"]:.4f})')

    print('\n(3) 순도 출처')
    print('  2차 런: 원문 인용 확보 (multiplexref Methods)')
    print(f'  1차 런: {PURITY_SOURCE["first_run"]}')
    print(f'  3차 런: {PURITY_SOURCE["third_run"]}')

    json.dump(rep, open('outputs/purity_sensitivity.json', 'w'), indent=1,
              ensure_ascii=False)
    print('\n-> outputs/purity_sensitivity.json')


if __name__ == '__main__':
    main()
