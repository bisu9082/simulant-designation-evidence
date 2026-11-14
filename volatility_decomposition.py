#!/usr/bin/env python3
r"""휘발성 축의 유의성이 무엇에서 오는가 (2026-11-14 신설).

왜 필요한가. 본 논문에서 우연과 구별되는 유일한 지정 검정이 휘발성 축의
V1 (사전지정 mean_diff, p=0.0206) 이다. 그 하나가 초록에 들어간다. 그런데
이 축에 대해 두 가지를 보고하지 않고 있었다.

(1) 측면 규약. 서술자 공간 세 곳에는 단측·양측 두 규약을 모두 싣는다
    (perm_sidedness.py, SI res:paneldetail). 휘발성 축에는 단측 하나만
    싣는다. 빠진 쪽이 우리에게 유리한 방향이다.

(2) 무엇이 통계량을 움직이는가. mean_diff 는 작용제별 '최근접 지정까지의
    거리 - 최근접 비지정까지의 거리' 의 평균이다. 1차원 축에서 지정 3종이
    축의 양끝과 중간을 잡으면 어느 작용제든 가까운 지정이 있으므로 이 양이
    작아진다. 즉 통계량이 충실도가 아니라 축 피복(span)을 재고 있을 수
    있다. 그 가설을 검정한다.

방법. 귀무·거리·통계량 정의는 volatility_space.py 를 그대로 임포트해 쓴다.
새 규칙을 쓰면 본문 값과 비교가 깨진다. 추가로 k=1,2 의 전수 열거를 돌려
현행 3종이 더 작은 지정 규모와 어떻게 비교되는지 적는다.

[C-3] 증기압은 volatility_space.VP 의 채택값 그대로이고, 새로 조회하거나
추정한 값은 없다.
"""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

from volatility_space import VP, space

PATH = 'descriptor_cache_xtb.csv'
SPAN_QUANTILES = (0.95, 0.90, 0.80)


def main():
    df = pd.read_csv(PATH)
    names = df.compound_name.tolist()
    kl = df.compound_class.tolist()
    D, x, z = space(names)
    D = D.astype(float).copy()
    np.fill_diagonal(D, np.inf)
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    pos = {g: j for j, g in enumerate(i_pool)}
    i_sim = [i for i, k in enumerate(kl) if k == 'simulant']
    true = tuple(sorted(pos[s] for s in i_sim))
    Ds = D[np.ix_(i_real, i_pool)]
    pool_names = [names[i] for i in i_pool]
    zp = z[i_pool]

    def mean_diff(combo):
        m = np.zeros(len(i_pool), bool)
        m[list(combo)] = True
        ds = np.where(m, Ds, np.inf).min(1)
        di = np.where(~m, Ds, np.inf).min(1)
        return float((ds - di).mean())

    out = {'generated_by': 'volatility_decomposition.py',
           'panel': 'panel28_incumbent3',
           'statistic': 'mean_diff (사전지정)',
           'note': ('휘발성 축 V1 의 측면 규약과 피복 분해. 본문 값을 '
                    '바꾸지 않고 그 값이 무엇에서 오는지 적는다.')}

    # ---- 측면 규약 세 가지 -------------------------------------------
    ALL3 = list(itertools.combinations(range(len(i_pool)), 3))
    null = np.array([mean_diff(c) for c in ALL3])
    obs = mean_diff(true)
    lo = float((null <= obs).mean())
    hi = float((null >= obs).mean())
    out['sidedness'] = dict(
        n_permutations=len(ALL3),
        observed=obs, null_mean=float(null.mean()),
        one_sided_low=lo,
        two_sided_2min=float(min(1.0, 2 * min(lo, hi))),
        two_sided_abs=float((np.abs(null - null.mean())
                             >= abs(obs - null.mean())).mean()),
        convention_source='perm_sidedness.py 와 같은 세 정의')
    print('[측면] 단측 {one_sided_low:.4f}  양측(2x) {two_sided_2min:.4f}  '
          '양측(편차) {two_sided_abs:.4f}'.format(**out['sidedness']))

    # ---- 지정 규모별 전수 열거 ---------------------------------------
    out['by_designation_size'] = {}
    for k in (1, 2, 3):
        ALL = list(itertools.combinations(range(len(i_pool)), k))
        vals = np.array([mean_diff(c) for c in ALL])
        order = np.argsort(vals)
        rows = []
        for j in order:
            c = ALL[j]
            rows.append(dict(
                compounds=sorted(pool_names[t] for t in c),
                mean_diff=float(vals[j]),
                rank=int(1 + (vals < vals[j]).sum()),
                p=float((vals <= vals[j]).mean()),
                span=float(zp[list(c)].max() - zp[list(c)].min())))
        rec = dict(n_combinations=len(ALL), ranking=rows)
        if k == 3:
            ti = ALL.index(true)
            rec['incumbent'] = dict(
                compounds=sorted(pool_names[t] for t in true),
                mean_diff=obs, rank=int(1 + (vals < obs).sum()), p=lo)
        out['by_designation_size'][f'k{k}'] = rec
        best = rows[0]
        print(f'[k={k}] C(17,{k})={len(ALL)}  최상 '
              f'{"+".join(best["compounds"])} ({best["mean_diff"]:+.4f})')

    # ---- 피복(span) 효과 ---------------------------------------------
    spans = np.array([zp[list(c)].max() - zp[list(c)].min() for c in ALL3])
    rho, p_rho = stats.spearmanr(spans, null)
    s_true = float(zp[list(true)].max() - zp[list(true)].min())
    lad = []
    for q in SPAN_QUANTILES:
        thr = float(np.quantile(spans, q))
        sel = spans >= thr
        lad.append(dict(quantile=q, threshold=thr, n=int(sel.sum()),
                        conditional_p=float((null[sel] <= obs).mean()),
                        rank=int(1 + (null[sel] < obs).sum())))
    sel_max = spans >= s_true
    out['span'] = dict(
        spearman_rho=float(rho), spearman_p=float(p_rho),
        n=len(ALL3),
        incumbent_span=s_true,
        incumbent_span_percentile=float(100 * (spans <= s_true).mean()),
        max_span=float(spans.max()),
        conditional_at_incumbent_span=dict(
            n=int(sel_max.sum()),
            conditional_p=float((null[sel_max] <= obs).mean()),
            rank=int(1 + (null[sel_max] < obs).sum()),
            caveat=('현행 span 이 최대값이므로 이 조건화는 표본 15개로 '
                    '거의 축퇴한다. 판정은 아래 사다리로 한다.')),
        ladder=lad,
        interpretation=('통계량이 축 피복의 단조함수에 가깝다. 같은 정도로 '
                        '퍼진 지정들 중에서 보면 현행 3종의 p 는 0.0206 이 '
                        '아니라 0.10~0.33 이다.'))
    print(f'[피복] Spearman {rho:.4f} (n={len(ALL3)})  '
          f'현행 span {s_true:.4f} = {out["span"]["incumbent_span_percentile"]:.0f}퍼센타일')
    for r in lad:
        print(f'       상위 {100 * (1 - r["quantile"]):.0f}% (span>={r["threshold"]:.4f}, '
              f'n={r["n"]}) 조건부 p = {r["conditional_p"]:.4f}')

    with open('outputs/volatility_decomposition.json', 'w',
              encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('  -> outputs/volatility_decomposition.json')


if __name__ == '__main__':
    main()
