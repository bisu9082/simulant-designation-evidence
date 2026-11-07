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
    elif reading == 'literal':
        # 2026-11-07 (12차 A-2). 세 번째 판독이다. 본문 실험절은 '지정집합
        # 밖의 모든 화합물' 이라고 적고, 작용제도 지정집합 밖이다. 종전의
        # wide/narrow 는 둘 다 작용제를 비교군에서 뺐으므로 문자 그대로의
        # 판독은 어느 쪽도 아니었고 감사되지 않았다. 자기 자신은 뺀다.
        i_cmp = [i for i, k in enumerate(kl) if k not in sim_cls]
    else:
        i_cmp = [i for i, k in enumerate(kl) if k == 'interferent']

    out = []
    for i in i_agent:
        js = min(i_sim, key=lambda j: D[i, j])
        d_norm = float(D[i, js] / nn_med)
        cmp_i = [j for j in i_cmp if j != i]
        beats = all(D[i, js] < D[i, j] for j in cmp_i)
        blocker = None
        if not beats:
            jb = min(cmp_i, key=lambda j: D[i, j])
            blocker = name[jb]
        out.append({'agent': name[i], 'nearest_simulant': name[js],
                    'd_norm': round(d_norm, 4),
                    'within_threshold': d_norm <= THRESHOLD,
                    'beats_comparison_set': beats, 'blocker': blocker,
                    'match': bool(d_norm <= THRESHOLD and beats)})
    return {'n_dims': len(keep), 'nn_median': round(nn_med, 4),
            'n_comparison': len(i_cmp),
            'n_match': sum(a['match'] for a in out), 'agents': out}


def panel_swap():
    """패널 교체가 무엇을 바꾸는가 (2026-11-07 신설, 12차 R2-F3).

    왜 필요한가. 본문은 패널 자유도를 "Adding back the five compounds the
    sensing literature also treats as simulants" 로 적고 "The incumbents do
    not move; what moves is whether anything closer is allowed into the
    comparison" 으로 닫았다. 둘 다 실측과 다르다.

      (1) 28종 -> 35종은 일곱 종을 더한다. 문헌 모사체 다섯에 더해
          간섭물 둘(dichlorvos, malathion)이 들어가고, 간섭물은 매칭
          2조건의 비교군이므로 후보와 무관하게 기준을 바꾼다. SI 의
          생성표는 이미 'omits the seven' 으로 적고 있었다.
      (2) 표준화는 패널 자신의 평균·표준편차로 하므로 추가된 화합물과
          아무 관계 없는 쌍의 거리도 바뀐다. 유지 좌표 집합도 바뀐다.
          따라서 지정 모사체도 움직인다.
    """
    import itertools as _it
    out = {}
    d28 = pd.read_csv('descriptor_cache_xtb.csv')
    d35 = pd.read_csv('descriptor_cache_xtb_35.csv')

    def prep(df):
        keep, _ = pick(df, merged_space(df), PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        Dx = D.copy()
        np.fill_diagonal(Dx, np.inf)
        return D, keep, float(np.median(Dx.min(1)))

    D28, k28, m28 = prep(d28)
    D35, k35, m35 = prep(d35)
    n28, n35 = list(d28.compound_name), list(d35.compound_name)
    i28 = {c: i for i, c in enumerate(n28)}
    i35 = {c: i for i, c in enumerate(n35)}
    added = d35[~d35.compound_name.isin(set(n28))]
    out['n_added'] = int(len(added))
    out['added'] = [{'name': r.compound_name, 'klass': r.compound_class}
                    for r in added.itertuples()]
    out['n_added_by_class'] = (added.compound_class.value_counts()
                               .to_dict())
    out['n_kept_28'] = len(k28)
    out['n_kept_35'] = len(k35)
    out['n_kept_common'] = len(set(k28) & set(k35))
    out['nn_median_28'] = m28
    out['nn_median_35'] = m35

    common = [c for c in n28 if c in i35]
    addnames = [r.compound_name for r in added.itertuples()]
    rel, far = [], []
    for a_, b_ in _it.combinations(common, 2):
        x = D28[i28[a_], i28[b_]] / m28
        y = D35[i35[a_], i35[b_]] / m35
        if x <= 0:
            continue
        rel.append(abs(y - x) / x)
        # 두 성분 각각에서 추가된 일곱 종까지의 최소 거리 중 작은 쪽
        far.append(min(min(D35[i35[a_], i35[c]] for c in addnames),
                       min(D35[i35[b_], i35[c]] for c in addnames)) / m35)
    rel = np.asarray(rel)
    far = np.asarray(far)
    out['common_pairs'] = int(len(rel))
    out['common_pair_rel_change_median'] = float(np.median(rel))
    out['common_pair_rel_change_max'] = float(rel.max())
    # 추가된 일곱 종에서 가장 먼 사분위. 효과가 국소적인지 본다.
    _cut = float(np.quantile(far, 0.75))
    _fq = rel[far >= _cut]
    out['far_quartile_n'] = int(len(_fq))
    out['far_quartile_rel_change_median'] = float(np.median(_fq))

    def nearest_inc(df, D, med):
        kl, nm = list(df.compound_class), list(df.compound_name)
        ia = [i for i, k in enumerate(kl) if k == 'real_agent']
        isim = [i for i, k in enumerate(kl) if k == 'simulant']
        res = {}
        for i in ia:
            j = min(isim, key=lambda j: D[i, j])
            res[nm[i]] = (nm[j], float(D[i, j] / med))
        return res

    a28, a35 = nearest_inc(d28, D28, m28), nearest_inc(d35, D35, m35)
    rows, rels = [], []
    for ag in sorted(a28):
        s0, x = a28[ag]
        s1, y = a35[ag]
        rels.append(abs(y - x) / x)
        rows.append({'agent': ag, 'nearest_28': s0, 'd_28': x,
                     'nearest_35': s1, 'd_35': y,
                     'identity_changed': s0 != s1,
                     'rel_change': abs(y - x) / x})
    out['incumbent_distance'] = rows
    out['incumbent_rel_change_median'] = float(np.median(rels))
    out['incumbent_rel_change_max'] = float(max(rels))
    out['incumbent_rel_change_max_agent'] = rows[
        int(np.argmax(rels))]['agent']
    out['n_identity_changed'] = int(sum(r['identity_changed'] for r in rows))

    print('\n[패널 교체] 28종 -> 35종')
    print(f"  추가 {out['n_added']}종 {out['n_added_by_class']}")
    print(f"  유지 좌표 {out['n_kept_28']} -> {out['n_kept_35']} "
          f"(공통 {out['n_kept_common']})")
    print(f"  공통 {out['common_pairs']}쌍 정규화 거리 상대변화: 중위 "
          f"{out['common_pair_rel_change_median'] * 100:.1f}%  최대 "
          f"{out['common_pair_rel_change_max'] * 100:.1f}%")
    print(f"  추가종에서 가장 먼 사분위 {out['far_quartile_n']}쌍: 중위 "
          f"{out['far_quartile_rel_change_median'] * 100:.1f}%")
    print(f"  작용제->최근접 지정모사체 거리 상대변화: 중위 "
          f"{out['incumbent_rel_change_median'] * 100:.1f}%  최대 "
          f"{out['incumbent_rel_change_max'] * 100:.1f}% "
          f"({out['incumbent_rel_change_max_agent']})")
    print(f"  최근접 상대가 바뀐 작용제 {out['n_identity_changed']}종")
    return out


def main():
    rep = {}
    for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        rep[label] = {}
        for rname, rule in (('deterministic', pick), ('legacy', _pick_legacy)):
            for reading in ('wide', 'narrow', 'literal'):
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
    rep['panel_swap'] = panel_swap()
    json.dump(rep, open('outputs/match_rule_audit.json', 'w'), indent=1,
              ensure_ascii=False)
    print('-> outputs/match_rule_audit.json')


if __name__ == '__main__':
    main()
