#!/usr/bin/env python3
"""실측 응답자료에서 더 끌어낼 수 있는 것 (2026-10-27 신설).

왜. deltaE_premise.py 가 전제(판별력 6/6 유의)와 전이검정(12셀, BH min_q
0.10, 유의 0건)과 DMMP 코사인을 낸다. 그런데 자료에 아직 쓰지 않은 축이
셋 남아 있고, 셋 다 모사체의 용도에 직접 걸린다. 새 측정 없이 여기서 낸다.

E1 농도 의존성. consensus() 가 네 농도를 평균해 하나로 접는다. 그런데
   현장 점탐지기는 낮은 농도에서 돈다. DMMP 가 500 에서만 작용제를 따라가고
   10 에서 안 따라가면 그 유사도는 용도에 대해 아무 말도 하지 않는다.
   접기 전 값을 농도별로 낸다. 추세는 4! = 24 전수치환으로 검정한다.

E2 응답공간 매칭. 본문 2.5절 매칭 규칙은 '모든 비모사체보다 가까움'을
   요구하는데 d1 에는 DMMP 말고 비작용제가 없어 그 규칙을 적용할 수 없다.
   대신 가환성 가정을 안 쓰는 질문으로 바꾼다. DMMP 가 어떤 작용제의
   대체물이라면, 그 작용제에서 DMMP 까지의 거리가 그 작용제에서 다른
   작용제까지의 거리들 사이에서 앞쪽에 앉아야 한다. 순위를 그대로 낸다.

E3 염료 수준 분해. 코사인은 29개 채널 기여의 합이다. DMMP-GA 의 +0.764 가
   몇 개 염료에서 나오는지 센다. 두세 개에 걸려 있으면 그 유사도는 약하다.

E4 데이터셋 간 결합. 기존 보고는 12셀 BH 로 '어느 셀이 유의한가'를 묻는다.
   '전이 증거가 전체적으로 있는가'는 다른 질문이고 결합검정이 답한다.
   d1 과 d2 는 물질도 플레이트도 겹치지 않아 독립이다. d3 는 서술자가 없어
   Mantel 셀에 없다. 따라서 결합은 2개 p 의 Fisher 이고, 공간 3 x 채널 2 의
   여섯 조합 전부를 낸 뒤 BH 를 건다. 사후 선택을 피하려고 전부 낸다.

E5 구조적 사실. 측정된 17종 중 지정 모사체가 몇 종인지, 어느 작용제군에
   모사체가 한 번도 함께 측정되지 않았는지.

출력: outputs/deltaE_extended.json
"""
import itertools
import json
import sys

import numpy as np
from scipy import stats

from deltaE_premise import (CONC, load, raw_patterns, centred, consensus)

INCUMBENT = {'DMMP', 'DCM', 'DPGME'}          # 본문의 현행 모사체 3종
D1_AGENTS = ['GA', 'GB', 'GD', 'GF', 'VX']
MODES = ['uv', 'normal']


def cos(u, v):
    return float(np.dot(u, v))


# ---------------------------------------------------------------- E1
def concentration_profile():
    """DMMP 대 각 작용제 코사인을 농도별로. 추세는 4! 전수치환."""
    out = []
    for mode in MODES:
        D, dy = load('d1', mode)
        cen = centred(raw_patterns(D, dy))
        if 'DMMP' not in cen:
            raise SystemExit('[중단] d1 에 DMMP 가 없다')
        for ag in D1_AGENTS:
            if ag not in cen:
                continue
            cs = [c for c in CONC if c in cen['DMMP'] and c in cen[ag]]
            vals = [cos(cen['DMMP'][c], cen[ag][c]) for c in cs]
            rec = dict(mode=mode, agent=ag, conc=cs,
                       cosines=[round(v, 4) for v in vals])
            if len(cs) >= 3:
                # 추세. 농도 순서와 코사인의 Spearman, 전수치환으로 양측 p.
                x = np.log10(cs)
                obs = stats.spearmanr(x, vals).statistic
                perms = list(itertools.permutations(range(len(cs))))
                null = [stats.spearmanr(x, [vals[i] for i in p]).statistic
                        for p in perms]
                null = np.array(null)
                rec['trend_rho'] = round(float(obs), 4)
                rec['trend_p_exact'] = round(
                    float((np.abs(null) >= abs(obs) - 1e-12).mean()), 4)
                rec['n_perm'] = len(perms)
                rec['cos_low'] = round(float(vals[0]), 4)
                rec['cos_high'] = round(float(vals[-1]), 4)
                rec['sign_flips'] = int(sum(
                    1 for a, b in zip(vals, vals[1:]) if a * b < 0))
            out.append(rec)
    return out


# ---------------------------------------------------------------- E2
def response_matching():
    """작용제 시점에서 본 DMMP 의 순위. 거리는 단위벡터 유클리드."""
    out = []
    for mode in MODES:
        D, dy = load('d1', mode)
        con = consensus(centred(raw_patterns(D, dy)))
        names = [n for n in con]
        def dist(a, b):
            return float(np.linalg.norm(con[a] - con[b]))
        for ag in D1_AGENTS:
            if ag not in con:
                continue
            others = [a for a in D1_AGENTS if a != ag and a in con]
            d_sim = dist(ag, 'DMMP')
            d_oth = sorted(dist(ag, o) for o in others)
            rank = int(sum(1 for d in d_oth if d < d_sim)) + 1
            out.append(dict(
                mode=mode, agent=ag,
                d_to_DMMP=round(d_sim, 4),
                d_to_nearest_agent=round(d_oth[0], 4),
                d_to_agents=[round(d, 4) for d in d_oth],
                rank_of_DMMP=rank, n_candidates=len(others) + 1,
                nearer_than_every_agent=bool(rank == 1)))
    # 귀무 대조. 후보가 k 개면 무작위 배정에서 1위 확률은 1/k 다. 10개 셀의
    # 1위 횟수를 기대값과 비교한다. 이것 없이 '두 셀에서 1위'를 쓰면
    # 우연과 구별되지 않는다.
    k = [r['n_candidates'] for r in out]
    exp = float(sum(1.0 / x for x in k))
    obs = int(sum(r['nearer_than_every_agent'] for r in out))
    # 이항 꼬리 (셀마다 k 가 같으므로 단일 p)
    n, p0 = len(out), 1.0 / k[0]
    tail = float(sum(stats.binom.pmf(i, n, p0) for i in range(obs, n + 1)))
    return dict(cells=out, n_cells=n, n_rank1_observed=obs,
                n_rank1_expected=round(exp, 2), p_binom_one_sided=round(tail, 4),
                channels_agree_on_best=sorted(
                    {r['agent'] for r in out
                     if r['mode'] == 'uv' and r['nearer_than_every_agent']}
                    & {r['agent'] for r in out
                       if r['mode'] == 'normal' and r['nearer_than_every_agent']}))


# ---------------------------------------------------------------- E3
def _concentration(u, v):
    """양의 기여 질량의 절반·80% 를 지는 채널 수와 최대 몫."""
    c = u * v
    pos = c[c > 0]
    if not len(pos):
        return 0, 0, None, 0
    srt = np.sort(pos)[::-1]
    cum, s = np.cumsum(srt), float(srt.sum())
    return (int(np.searchsorted(cum, 0.5 * s) + 1),
            int(np.searchsorted(cum, 0.8 * s) + 1),
            float(srt[0] / s), int((c > 0).sum()))


def dye_decomposition():
    """코사인을 채널 기여로 분해. 양의 합의 50% / 80% 를 몇 채널이 지는가.

    귀무 대조가 둘 필요하다 (2026-10-27). 첫째, 작용제-작용제 쌍도 같은 식으로
    소수 채널에 몰려 있으면 집중도는 어레이의 성질이고 DMMP 의 성질이 아니다.
    둘째, 중심화한 단위벡터 두 개를 곱하면 기여가 얼마나 몰리는 것이 보통인지
    모르므로, 채널 라벨을 섞은 치환 귀무를 둔다.
    """
    rows, controls = [], []
    rng = np.random.default_rng(20261027)
    for mode in MODES:
        D, dy = load('d1', mode)
        con = consensus(centred(raw_patterns(D, dy)))
        n_ch = len(con['DMMP'])
        ags = [x for x in D1_AGENTS if x in con]
        for ag in ags:
            n50, n80, top, npos = _concentration(con['DMMP'], con[ag])
            rows.append(dict(
                mode=mode, agent=ag, n_channels=int(n_ch),
                cosine=round(cos(con['DMMP'], con[ag]), 4), n_positive=npos,
                top1_share=round(top, 4) if top else None,
                n_channels_for_half=n50, n_channels_for_80pct=n80))
        # 대조 1. 작용제-작용제 쌍 전부. 집중도가 어레이의 성질인지 가른다.
        aa = []
        for a, b in itertools.combinations(ags, 2):
            n50, n80, top, npos = _concentration(con[a], con[b])
            aa.append(dict(pair=f'{a}-{b}', n_channels_for_half=n50,
                           n_channels_for_80pct=n80,
                           top1_share=round(top, 4) if top else None,
                           cosine=round(cos(con[a], con[b]), 4)))
        # 대조 2. 채널 라벨 치환 귀무. DMMP 벡터의 채널 순서만 섞는다.
        null50 = [_concentration(con['DMMP'][rng.permutation(n_ch)],
                                 con[ag])[0]
                  for _ in range(2000) for ag in ags]
        obs50 = [r['n_channels_for_half'] for r in rows if r['mode'] == mode]
        controls.append(dict(
            mode=mode, agent_pairs=aa,
            agent_pair_half_median=float(np.median(
                [x['n_channels_for_half'] for x in aa])),
            dmmp_half_median=float(np.median(obs50)),
            null_half_mean=round(float(np.mean(null50)), 2),
            null_half_median=float(np.median(null50)),
            null_half_p05=int(np.percentile(null50, 5)),
            n_null=len(null50)))
    return dict(cells=rows, controls=controls)


# ---------------------------------------------------------------- E4
def combine_datasets():
    """d1 과 d2 의 전수 Mantel p 를 Fisher 로 결합. 여섯 조합 전부 + BH."""
    cells = json.load(open('outputs/rollback_c.json'))['deltaE']['mantel']
    idx = {(c['dataset'], c['mode'], c['space']): c['p_mantel'] for c in cells}
    spaces = sorted({c['space'] for c in cells})
    rows = []
    for sp in spaces:
        for mode in MODES:
            ps = [idx.get(('d1', mode, sp)), idx.get(('d2', mode, sp))]
            if any(p is None for p in ps):
                continue
            # 전수치환 p 는 이산이고 0 이 될 수 있다. 여기서는 0 이 없다.
            stat, pf = stats.combine_pvalues(ps, method='fisher')
            rows.append(dict(space=sp, mode=mode,
                             p_d1=round(ps[0], 4), p_d2=round(ps[1], 4),
                             fisher_stat=round(float(stat), 4),
                             fisher_p=round(float(pf), 4)))
    # BH over the six
    pv = np.array([r['fisher_p'] for r in rows])
    order = pv.argsort()
    m = len(pv)
    q = np.empty(m)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = m - rank
        prev = min(prev, pv[i] * m / k)
        q[i] = prev
    for r, qq in zip(rows, q):
        r['bh_q'] = round(float(qq), 4)
    return dict(cells=rows, n_tests=m, min_q=round(float(q.min()), 4),
                n_sig_q05=int((q < 0.05).sum()))


# ---------------------------------------------------------------- E5
def structural_facts():
    measured, per = set(), {}
    for tag in ('d1', 'd2', 'd3'):
        D, _ = load(tag, 'uv')
        cs = sorted(D.Agent.unique())
        per[tag] = cs
        measured |= set(cs)
    sims = sorted(measured & INCUMBENT)
    return dict(
        n_measured=len(measured), compounds=sorted(measured),
        per_dataset=per,
        incumbents_measured=sims,
        n_incumbents_measured=len(sims),
        n_incumbents_total=len(INCUMBENT),
        datasets_with_an_incumbent=[t for t, cs in per.items()
                                    if set(cs) & INCUMBENT],
        datasets_without_an_incumbent=[t for t, cs in per.items()
                                       if not set(cs) & INCUMBENT],
        # d1 에 DMMP 외 비작용제가 없으므로 지정 치환 귀무가 비어 있다.
        d1_alternative_candidates=[c for c in per['d1']
                                   if c not in D1_AGENTS and c != 'DMMP'],
        n_designation_permutations_d1=0)


def main():
    res = dict(
        concentration=concentration_profile(),
        matching=response_matching(),
        dyes=dye_decomposition(),
        combined=combine_datasets(),
        facts=structural_facts())

    print('=== E1 농도 의존성 (DMMP 대 작용제 코사인) ===')
    for r in res['concentration']:
        cs = '  '.join(f'{c}:{v:+.3f}' for c, v in
                       zip(r['conc'], r['cosines']))
        t = (f"  추세 rho {r['trend_rho']:+.2f} p {r['trend_p_exact']:.3f}"
             if 'trend_rho' in r else '')
        print(f"  {r['mode']:6s} {r['agent']:3s}  {cs}{t}")

    M = res['matching']
    print('\n=== E2 응답공간에서 본 DMMP 순위 ===')
    for r in M['cells']:
        print(f"  {r['mode']:6s} {r['agent']:3s}  DMMP까지 {r['d_to_DMMP']:.3f}  "
              f"최근접작용제 {r['d_to_nearest_agent']:.3f}  "
              f"순위 {r['rank_of_DMMP']}/{r['n_candidates']}"
              f"{'  <- 1위' if r['nearer_than_every_agent'] else ''}")
    print(f"  1위 {M['n_rank1_observed']}/{M['n_cells']} 셀, "
          f"무작위 기대 {M['n_rank1_expected']}, "
          f"이항 일측 p {M['p_binom_one_sided']:.3f}")
    print(f"  두 채널이 같은 작용제에서 1위라고 하는 경우: "
          f"{M['channels_agree_on_best'] or '없음'}")

    print('\n=== E3 염료 수준 분해 ===')
    for r in res['dyes']['cells']:
        print(f"  {r['mode']:6s} {r['agent']:3s}  cos {r['cosine']:+.3f}  "
              f"양의채널 {r['n_positive']:2d}/{r['n_channels']}  "
              f"절반을 지는 채널 {r['n_channels_for_half']}  "
              f"80%를 지는 채널 {r['n_channels_for_80pct']}  "
              f"최대채널몫 {r['top1_share']}")
    for c in res['dyes']['controls']:
        print(f"  [대조 {c['mode']:6s}] DMMP 쌍 중위 {c['dmmp_half_median']:.1f}  "
              f"작용제-작용제 쌍 중위 {c['agent_pair_half_median']:.1f}  "
              f"채널치환 귀무 평균 {c['null_half_mean']:.2f} "
              f"중위 {c['null_half_median']:.1f} 5백분위 {c['null_half_p05']}")

    print('\n=== E4 d1+d2 Fisher 결합 ===')
    for r in res['combined']['cells']:
        print(f"  {r['space']:12s} {r['mode']:6s}  d1 {r['p_d1']:.3f}  "
              f"d2 {r['p_d2']:.3f}  ->  Fisher p {r['fisher_p']:.4f}  "
              f"BH q {r['bh_q']:.4f}")
    print(f"  최소 q {res['combined']['min_q']:.4f}, "
          f"q<0.05 셀 {res['combined']['n_sig_q05']}/"
          f"{res['combined']['n_tests']}")

    f = res['facts']
    print('\n=== E5 구조적 사실 ===')
    print(f"  측정된 화합물 {f['n_measured']}종")
    print(f"  그중 현행 모사체 {f['n_incumbents_measured']}/"
          f"{f['n_incumbents_total']}종: {f['incumbents_measured']}")
    print(f"  모사체가 함께 측정된 데이터셋 {f['datasets_with_an_incumbent']}")
    print(f"  모사체 없는 데이터셋 {f['datasets_without_an_incumbent']}")
    print(f"  d1 의 대안 후보 모사체 {f['d1_alternative_candidates']} "
          f"-> 지정 치환 귀무 {f['n_designation_permutations_d1']}개")

    json.dump(res, open('outputs/deltaE_extended.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n[done] outputs/deltaE_extended.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
