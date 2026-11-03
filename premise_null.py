#!/usr/bin/env python3
r"""전제검정의 귀무 교체 — 라벨 순열 (2026-11-03 신설, 작업 A-2/A-5).

왜 교체하는가. 두 가지가 겹쳐 있었다.

(1) 중심화 '바닥' 논증이 계산 대상에 적용되지 않는다. deltaE_premise.centred()
    는 물질축 평균을 뺀 뒤 unit() 으로 다시 정규화한다. 재정규화가
    Sum_i u_i = 0 을 깨므로, -1/(n-1) 을 고정하는 항등식은 논문이 코사인을
    계산하는 벡터에 성립하지 않는다. 덧붙여 -1/(n-1) 은 영합 벡터의
    노름가중 평균 코사인에 대한 상한이고 등노름에서만 등호다
    (Sum_{i<j}||u_i||||u_j|| <= (n-1)/2 Sum||u_i||^2 에서 바로 나온다).
    즉 '바닥' 이라는 이름부터 방향이 반대다. 따라서 "관측 평균이 바닥을
    0.0032~0.1576 넘는다" 는 적용되지 않는 경계와의 비교였다.

(2) 부호뒤집기 전수열거가 여섯 셀 전부에서 자기 바닥을 돌려준다
    (p == p_min_attainable, 2^-6/2^-7/2^-4). d_a = w_a - b_a 의 부호는
    b_a 가 중심화 때문에 음수 상수 근처에 묶여 있어 거의 모든 비퇴화
    데이터에서 전부 양수가 된다. 그러면 p 는 데이터가 아니라 설계의 함수다.
    바닥 위 여유도(0.0032~0.1576)가 주장을 담을 양인데 거기에는 추론이
    붙어 있지 않았다.

교체 설계. 통계량은 그대로 둔다 — 화합물 a 마다 w_a(농도간 코사인 평균)와
b_a(합의패턴 간 코사인 평균), 검정통계량은 mean_a (w_a - b_a). 바꾸는 것은
귀무다. 각 농도 안에서 화합물 라벨을 순열한다. 중심화는 농도별로 화합물축
평균을 빼므로, 농도 안의 벡터 집합은 순열해도 불변이고 평균도 불변이다.
따라서 중심화 벡터는 그대로이고 배정만 바뀐다. 이 귀무가 묻는 것은
정확히 "한 화합물의 어떤 농도 응답이 그 화합물의 다른 농도 응답과
관계가 있는가" 이며, 그것이 전제가 주장하는 바다.

전수열거 가능성. 첫 농도의 라벨을 고정하면 서로 다른 상대 배정은
(n!)^(C-1) 가지다. d3 은 4종 4농도이므로 24^3 = 13824 로 전수 가능하고,
d1 은 720^3 = 3.7억이라 표집한다. 전수인 셀과 표집인 셀을 구분해 적는다
(논문이 '전수열거' 를 범위 없이 주장했던 결함을 되풀이하지 않는다).

산출: outputs/premise_null.json
"""
import itertools
import json
import math
import os

import numpy as np

import deltaE_premise as m

SEED = 49374
N_SAMPLE = 10000
EXACT_MAX = 200000


def cell_vectors(tag, mode):
    """농도별 중심화 단위벡터 행렬. (농도 -> (화합물목록, V))"""
    cen = m.centred(m.raw_patterns(*m.load(tag, mode)))
    names = sorted(cen)
    concs = sorted({c for a in names for c in cen[a]})
    out = {}
    for c in concs:
        here = [a for a in names if c in cen[a]]
        if len(here) < 2:
            continue
        out[c] = (here, np.array([cen[a][c] for a in here], float))
    return names, out


def statistic(names, per_conc, assign):
    """mean_a (w_a - b_a). assign[c] 는 슬롯 -> 화합물 이름."""
    vecs = {a: [] for a in names}
    for c, (here, V) in per_conc.items():
        lab = assign[c]
        for s in range(len(here)):
            vecs[lab[s]].append(V[s])
    cons, within = {}, {}
    for a, vs in vecs.items():
        if len(vs) < 2:
            continue
        A = np.array(vs)
        cons[a] = A.mean(0) / np.linalg.norm(A.mean(0))
        within[a] = float(np.mean([float(A[i] @ A[j])
                                   for i, j in itertools.combinations(
                                       range(len(A)), 2)]))
    ks = sorted(cons)
    if len(ks) < 2:
        return float('nan'), {}, {}
    between = {a: float(np.mean([float(cons[a] @ cons[o])
                                 for o in ks if o != a])) for a in ks}
    diff = [within[a] - between[a] for a in ks]
    return float(np.mean(diff)), within, between


def run_cell(tag, mode, rng):
    names, per_conc = cell_vectors(tag, mode)
    concs = sorted(per_conc)
    if len(concs) < 2:
        return None
    base = {c: list(per_conc[c][0]) for c in concs}
    obs, within, between = statistic(names, per_conc, base)

    sizes = [len(per_conc[c][0]) for c in concs[1:]]
    n_total = 1
    for s in sizes:
        n_total *= math.factorial(s)
    exact = n_total <= EXACT_MAX

    null = []
    if exact:
        pools = [list(itertools.permutations(range(s))) for s in sizes]
        for combo in itertools.product(*pools):
            a = {concs[0]: list(per_conc[concs[0]][0])}
            for c, pm in zip(concs[1:], combo):
                here = per_conc[c][0]
                a[c] = [here[i] for i in pm]
            null.append(statistic(names, per_conc, a)[0])
    else:
        for _ in range(N_SAMPLE):
            a = {concs[0]: list(per_conc[concs[0]][0])}
            for c in concs[1:]:
                here = list(per_conc[c][0])
                rng.shuffle(here)
                a[c] = here
            null.append(statistic(names, per_conc, a)[0])
    null = np.array([x for x in null if x == x])
    # 단측. 통계량이 클수록 변별력이 크다.
    if exact:
        p = float((null >= obs).sum() / len(null))
    else:
        p = float(((null >= obs).sum() + 1) / (len(null) + 1))

    n = len(names)
    mean_between = float(np.mean(list(between.values())))
    return {'dataset': tag, 'mode': mode, 'n_compounds': n,
            'n_concentrations': len(concs),
            'observed': obs,
            'null_mean': float(null.mean()), 'null_sd': float(null.std(ddof=1)),
            'p_one_sided': p,
            'enumerated': bool(exact),
            'n_null': int(len(null)),
            'n_distinct_labelings': int(n_total),
            'p_min_attainable': float(1.0 / len(null)) if exact else None,
            'mean_within': float(np.mean(list(within.values()))),
            'mean_between': mean_between,
            # 기록용. 이 경계는 재정규화 때문에 계산 대상에 적용되지 않는다.
            'minus_one_over_n_minus_1': float(-1.0 / (n - 1)),
            'note_bound': ('centred() 가 unit(v-mu) 로 재정규화하므로 '
                           'Sum u = 0 이 깨진다. -1/(n-1) 은 영합 벡터의 '
                           '노름가중 평균에 대한 상한이며 이 양의 바닥이 '
                           '아니다. 비교를 위해 적어 두기만 한다.')}


def main():
    rng = np.random.default_rng(SEED)
    print('=' * 78)
    print('[A-2/A-5] 전제검정 귀무 교체 — 농도 내 라벨 순열')
    print('=' * 78)
    rows = []
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                r = run_cell(tag, mode, rng)
            except Exception as e:
                print(f'  {tag}_{mode}: [수집 불가: {e}]')
                continue
            if r is None:
                continue
            rows.append(r)
            kind = '전수' if r['enumerated'] else f"표집 {r['n_null']}"
            print(f"  {tag}_{mode:<7} n={r['n_compounds']} "
                  f"관측 {r['observed']:+.4f}  귀무 {r['null_mean']:+.4f}"
                  f"±{r['null_sd']:.4f}  p={r['p_one_sided']:.4f}  ({kind}, "
                  f"가능 배정 {r['n_distinct_labelings']:,})")
    # 여섯 셀은 하나의 family 다. 본문이 18검정·72검정을 보정하면서 이쪽을
    # 보정하지 않으면 같은 오류다.
    o = np.sort(np.array([r['p_one_sided'] for r in rows]))
    K = len(o)
    bh = np.minimum.accumulate((o * K / np.arange(1, K + 1))[::-1])[::-1]
    print(f'\n  {K}셀 BH 최소 q = {bh[0]:.4f}  '
          f'유의 {int((bh < 0.05).sum())}건')
    os.makedirs('outputs', exist_ok=True)
    with open('outputs/premise_null.json', 'w') as fh:
        json.dump({'generated_by': 'premise_null.py', 'seed': SEED,
                   'statistic': 'mean over compounds of (within - between)',
                   'null': 'compound labels permuted within each '
                           'concentration',
                   'cells': rows,
                   'bh': {'n_cells': K, 'min_q': float(bh[0]),
                          'n_sig_q': int((bh < 0.05).sum())}},
                  fh, indent=1, ensure_ascii=False)
    print('-> outputs/premise_null.json')


if __name__ == '__main__':
    main()
