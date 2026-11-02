#!/usr/bin/env python3
"""매칭 규칙의 두 판독 대조 (2026-11-02 신설, T2-2).

왜 필요한가. 본문 실험절은 유효매칭을 '정규화거리 2 이하이고 지정집합 밖의
모든 화합물보다 가깝다'로 정의한다. 코드에는 그 둘째 조건의 비교군이 두 가지로
구현되어 있었다.

  넓은 판독 (match_by_panel.py): 작용제도 지정모사체도 아닌 모든 화합물
  좁은 판독 (analysis_main.py): compound_class == 'interferent' 인 화합물

28종 패널에서는 두 판독이 같은 집합을 가리킨다. 그 패널에는 interferent 14 +
real_agent 11 + simulant 3 만 있어 '지정집합 밖'이 곧 interferent 다. 35종
패널에는 literature_simulant 5 종이 더 있으므로 두 판독이 갈라진다. 본문이
'주 패널'이라 부르며 매칭 0 건을 보고하는 쪽이 바로 그 35종 패널이다.

이 스크립트는 어느 쪽이 옳은지 고르지 않고 두 수를 모두 센다. 문장과 코드를
맞추는 일은 둘을 본 다음이다.

사용: python3 match_rule_audit.py
"""
import json

import numpy as np
import pandas as pd

from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, _pick_legacy, pick

THRESHOLD = 2.0
CASES = [
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
]


def merged_space(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return [c for c in num if c != OMEGA_DERIVED]


def run(df, sim_cls, reading, rule):
    keep, _ = rule(df, merged_space(df), PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))

    kl, name = df.compound_class.tolist(), df.compound_name.tolist()
    i_agent = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    if reading == 'wide':
        i_cmp = [i for i, k in enumerate(kl)
                 if k != 'real_agent' and k not in sim_cls]
    else:
        i_cmp = [i for i, k in enumerate(kl) if k == 'interferent']

    out = []
    for i in i_agent:
        js = min(i_sim, key=lambda j: D[i, j])
        d_norm = float(D[i, js] / nn_med)
        beats = all(D[i, js] < D[i, j] for j in i_cmp)
        blocker = None
        if not beats:
            jb = min(i_cmp, key=lambda j: D[i, j])
            blocker = name[jb]
        out.append({'agent': name[i], 'nearest_simulant': name[js],
                    'd_norm': round(d_norm, 4),
                    'within_threshold': d_norm <= THRESHOLD,
                    'beats_comparison_set': beats, 'blocker': blocker,
                    'match': bool(d_norm <= THRESHOLD and beats)})
    return {'n_dims': len(keep), 'nn_median': round(nn_med, 4),
            'n_comparison': len(i_cmp),
            'n_match': sum(a['match'] for a in out), 'agents': out}


def main():
    rep = {}
    for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        rep[label] = {}
        for rname, rule in (('deterministic', pick), ('legacy', _pick_legacy)):
            for reading in ('wide', 'narrow'):
                r = run(df, sim_cls, reading, rule)
                rep[label][f'{rname}/{reading}'] = r
                print(f'{label:<30} {rname:<14} {reading:<7} '
                      f'dims {r["n_dims"]:>3}  비교군 {r["n_comparison"]:>2}  '
                      f'매칭 {r["n_match"]}')
        # 두 판독이 갈라지는 작용제를 짚는다.
        for rname in ('deterministic', 'legacy'):
            w = {a['agent']: a['match'] for a in rep[label][f'{rname}/wide']['agents']}
            n = {a['agent']: a['match'] for a in rep[label][f'{rname}/narrow']['agents']}
            diff = [k for k in w if w[k] != n[k]]
            if diff:
                print(f'    [{rname}] 판독에 따라 갈리는 작용제: {", ".join(diff)}')
        print()
    json.dump(rep, open('outputs/match_rule_audit.json', 'w'), indent=1,
              ensure_ascii=False)
    print('-> outputs/match_rule_audit.json')


if __name__ == '__main__':
    main()
