#!/usr/bin/env python3
"""자기보다 남이 가까운 사례의 계수 (2026-11-07 신설, 12차 GATE 8).

왜 필요한가. 원고 main.tex:354 와 supporting.tex:1606 이 "13 of the 34
compound-and-cell instances" 를 인용하는데, 어느 산출물에도 그 수가 없고
어느 스크립트도 그것을 계산하지 않았다. 수치 자체는 재계산으로 맞았으나
재현 경로가 없었다. 여기서 그 경로를 만든다.

무엇을 세는가. 한 (화합물, 셀) 사례에서

  자기일치  = 그 화합물의 농도별 응답 패턴끼리의 평균 코사인
  최근접 남 = 그 화합물의 합의 패턴과 다른 화합물 합의 패턴의 최대 코사인

이고, 최근접 남 > 자기일치 이면 '역전' 1건으로 센다. 역전은 그 표현이
화합물 정체성을 자기 반복측정만큼도 보존하지 못한다는 뜻이다.

두 표현에서 센다. 이것이 이 스크립트의 핵심이다.
  (1) 중심화 전 : 단위정규화만 한 응답. 원천 논문 분류기가 보는 것에 가깝다.
  (2) 중심화 후 : 농도별로 물질축 평균을 뺀 응답. 원고의 모든 분석이 쓴다.
중심화가 역전을 만드는지 줄이는지를 수치로 답해야 한다. 방향 A 가 중심화
인공물로 무너진 뒤이므로, 이 통계도 같은 검사를 통과해야 쓸 수 있다.

귀무. 라벨순열 - 농도마다 화합물 라벨을 뒤섞는다. 중심화 벡터 자체는
건드리지 않고 어느 응답이 어느 화합물에 배정되는지만 바꾼다. 원고
deltaE_premise 가 쓰는 것과 같은 귀무다. 화합물 수와 농도 수가 작은
셀은 전수열거하고, 그렇지 않은 셀은 시드 5개로 각 10,000회 표집해
흩어짐을 함께 보고한다.

출력: outputs/reversal_audit.json
"""
import itertools
import json
import os
from itertools import permutations as _perm

import numpy as np

import deltaE_premise as m

SEEDS = [0xC0DE, 1, 2, 7, 1234]   # 시드 5개 (통계 엄밀성 최소 기준)
NDRAW = 10000                      # 표집 셀의 1시드당 추출 수
EXACT_CAP = 200000                 # 전수열거를 허용하는 상한

CELLS = [(t, md) for t in ('d1', 'd2', 'd3') for md in ('uv', 'normal')]


def consensus_of(pat):
    return {a: m.unit(np.mean([pat[a][c] for c in sorted(pat[a])], axis=0))
            for a in pat}


def reversals(pat):
    """(역전수, 사례수, 화합물별 상세)."""
    P = consensus_of(pat)
    names = sorted(P)
    hit, det = 0, []
    for a in names:
        cs = sorted(pat[a])
        ws = [float(np.dot(pat[a][i], pat[a][j]))
              for i, j in itertools.combinations(cs, 2)]
        if not ws:
            continue
        within = float(np.mean(ws))
        best_v, best_o = -9.0, None
        for b in names:
            if b == a:
                continue
            c = float(np.dot(P[a], P[b]))
            if c > best_v:
                best_v, best_o = c, b
        rev = bool(best_v > within)
        hit += int(rev)
        det.append(dict(compound=a, within=within,
                        nearest_other=best_o, nearest_cos=best_v,
                        reversed=rev))
    return hit, len(det), det


def relabel(pat, perms_by_conc):
    """농도마다 주어진 순열로 라벨을 재배정한다. 벡터는 그대로다."""
    names = sorted(pat)
    out = {a: {} for a in names}
    for c, pmx in perms_by_conc.items():
        present = [a for a in names if c in pat[a]]
        for i, a in enumerate(present):
            out[a][c] = pat[present[pmx[i]]][c]
    return {a: d for a, d in out.items() if d}


def null_counts(pat, rng=None, exact=False):
    """귀무에서의 역전수 분포."""
    names = sorted(pat)
    concs = sorted({c for a in names for c in pat[a]})
    sizes = {c: sum(1 for a in names if c in pat[a]) for c in concs}
    if exact:
        pools = {c: list(_perm(range(sizes[c]))) for c in concs}
        out = []
        for combo in itertools.product(*[pools[c] for c in concs]):
            pb = {c: combo[i] for i, c in enumerate(concs)}
            out.append(reversals(relabel(pat, pb))[0])
        return out
    out = []
    for _ in range(NDRAW):
        pb = {c: tuple(rng.permutation(sizes[c])) for c in concs}
        out.append(reversals(relabel(pat, pb))[0])
    return out


def _fact(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def exact_space(pat):
    """(재배정 가지수, 전수열거 가능 여부)."""
    names = sorted(pat)
    concs = sorted({c for a in names for c in pat[a]})
    n = 1
    for c in concs:
        n *= _fact(sum(1 for a in names if c in pat[a]))
    return n, n <= EXACT_CAP


# 비작용제는 DMMP 하나다. 2차 런의 일곱 종은 원천 논문이 모두 작용제로
# 세고 (main.tex:285 "all seven agents"), 그 가운데 HN 은 서술자 벡터가
# 없어 Mantel 에서만 빠진다 (main.tex:453). 3차 런의 A-230/232/234/242 도
# 작용제다. 따라서 "지정 모사체가 작용제와 함께 계측된" 셀은 1차 런
# 두 채널뿐이고, 이것이 main.tex:301 의 서술과 일치한다.
SIMULANTS = {'DMMP'}


def pair_baseline():
    """셀 기저선. 작용제끼리의 쌍 분포 대 비작용제 각자의 대작용제 분포.

    왜 필요한가. 중심화 코사인에서 음의 부호는 중심화가 만든다. 그래서
    "이 화합물이 저 화합물에 반대로 응답한다" 는 진술은 같은 셀의 임의
    쌍이 어디에 있는지를 함께 적지 않으면 뜻이 없다. 11차 판본의 헤드라인
    이 그 대조 없이 세워져 있었다.
    """
    out = {}
    for tag, mode in CELLS:
        try:
            D, dy = m.load(tag, mode)
        except Exception as e:
            out[f'{tag}.{mode}'] = dict(unavailable=str(e))
            continue
        P = consensus_of(m.centred(m.raw_patterns(D, dy)))
        nm = sorted(P)
        sims = [a for a in nm if a in SIMULANTS]
        ags = [a for a in nm if a not in SIMULANTS]
        if len(ags) < 3 or not sims:
            out[f'{tag}.{mode}'] = dict(n_agents=len(ags),
                                        non_agents=sims,
                                        note='작용제 쌍 또는 비작용제 부족')
            continue
        aa = {}
        for a, b in itertools.combinations(ags, 2):
            aa[(a, b)] = float(np.dot(P[a], P[b]))
        av = list(aa.values())
        worst = min(aa.items(), key=lambda kv: kv[1])
        blk = dict(
            n_agents=len(ags), agents=ags, n_agent_pairs=len(av),
            agent_pair_mean=float(np.mean(av)),
            agent_pair_median=float(np.median(av)),
            agent_pair_min=float(min(av)), agent_pair_max=float(max(av)),
            most_negative_pair=list(worst[0]),
            most_negative_value=float(worst[1]),
            non_agents={})
        for s in sims:
            sv = [float(np.dot(P[s], P[a])) for a in ags]
            blk['non_agents'][s] = dict(
                n_pairs=len(sv), mean=float(np.mean(sv)),
                median=float(np.median(sv)),
                per_agent={a: float(np.dot(P[s], P[a])) for a in ags},
                n_agent_pairs_below_its_mean=int(
                    sum(1 for x in av if x < np.mean(sv))),
                above_agent_pair_mean=bool(np.mean(sv) >
                                           np.mean(av)))
        out[f'{tag}.{mode}'] = blk
    return out


def inter_channel():
    """같은 플레이트 두 채널의 평균 상호 코사인.

    왜 여기에 있는가. SI 가 이 양을 0.584 로 두 곳에서 인용하는데
    (채널 두 개를 한 Fisher 결합에 넣지 않는 근거다), 그 값은
    volatility_space.py 의 주석과 rationale 문자열에 리터럴로만 있었고
    계산하는 코드가 없었다. 실제로 재보니 0.584 는 중심화 전 값이고,
    검정이 쓰는 중심화 표현에서는 0.353 이다. 두 값을 모두 낸다.

    두 채널의 유효 채널 수가 다를 수 있으므로 (축약은 채널별로 걸린다)
    공통 채널로 제한한 뒤 비교한다.
    """
    cen, raw = [], []
    per_plate = {}
    for tag in ('d1', 'd2', 'd3'):
        try:
            Du, du = m.load(tag, 'uv')
            Dn, dn = m.load(tag, 'normal')
        except Exception as e:
            per_plate[tag] = dict(unavailable=str(e))
            continue
        common = [c for c in du if c in dn]
        Pu = consensus_of(m.centred(m.raw_patterns(Du, common)))
        Pn = consensus_of(m.centred(m.raw_patterns(Dn, common)))
        Ru = consensus_of(m.raw_patterns(Du, common))
        Rn = consensus_of(m.raw_patterns(Dn, common))
        nm = sorted(set(Pu) & set(Pn))
        cv = [float(np.dot(Pu[a], Pn[a])) for a in nm]
        rv = [float(np.dot(Ru[a], Rn[a])) for a in nm]
        cen += cv
        raw += rv
        per_plate[tag] = dict(n_common_channels=len(common),
                              n_compounds=len(nm),
                              centred_mean=float(np.mean(cv)),
                              uncentred_mean=float(np.mean(rv)))
    return dict(n_pairs=len(cen),
                centred_mean=float(np.mean(cen)) if cen else None,
                uncentred_mean=float(np.mean(raw)) if raw else None,
                per_plate=per_plate)


def collapse_free_counts():
    """중복 채널 축약을 끈 상태에서 같은 계수를 낸다.

    이 논문의 다른 통계는 전부 분석자 자유도에서 움직였다. 응답공간
    통계에 영향을 주는 자유도는 축약 하나이므로 (패널 구성·서술자 집합·
    임계값·부호화는 서술자 쪽이고 응답 패턴을 건드리지 않는다) 그것만
    끈다. collapse_sensitivity.py 와 같은 자리바꿈을 쓴다.
    """
    _orig = m.collapse_duplicate_channels
    m.collapse_duplicate_channels = lambda D, dy, key: list(dy)
    out = {}
    try:
        for tag, mode in CELLS:
            try:
                D, dy = m.load(tag, mode)
            except Exception as e:
                out[f'{tag}.{mode}'] = dict(unavailable=str(e))
                continue
            raw = m.raw_patterns(D, dy)
            hu, nu, _ = reversals(raw)
            hc, nc, _ = reversals(m.centred(raw))
            out[f'{tag}.{mode}'] = dict(n_channels=len(dy), n_instances=nu,
                                        reversals_uncentred=hu,
                                        reversals_centred=hc)
    finally:
        m.collapse_duplicate_channels = _orig
    return out


def main():
    os.makedirs('outputs', exist_ok=True)
    # n_draws 는 검증기의 표집 탐지가 보는 이름이다 (2026-11-07,
    # 12차). n_draws_per_seed 로만 적으면 이 산출물이 전수열거로
    # 오인되고, '전수' 주장 가드가 그만큼 헐거워진다.
    report = dict(generated_by='reversal_audit.py', seeds=SEEDS,
                  n_seeds=len(SEEDS), n_draws=NDRAW,
                  n_draws_per_seed=NDRAW, cells=[])
    tot = {'uncentred': 0, 'centred': 0}
    inst = 0

    print('=' * 78)
    print('[역전 감사] 자기 반복측정보다 남이 가까운 (화합물, 셀) 사례')
    print('=' * 78)

    for tag, mode in CELLS:
        try:
            D, dy = m.load(tag, mode)
        except Exception as e:
            print(f'  [수집 불가: {tag} {mode} {e}]')
            report['cells'].append(dict(dataset=tag, mode=mode,
                                        unavailable=str(e)))
            continue
        raw = m.raw_patterns(D, dy)
        cen = m.centred(raw)
        hu, nu, du = reversals(raw)
        hc, nc, dc = reversals(cen)
        tot['uncentred'] += hu
        tot['centred'] += hc
        inst += nu

        # 본문이 인용하는 두 범위를 여기서 낸다.
        #  (1) 중심화 전 쌍평균: 모든 화합물이 서로 닮는 수준
        #  (2) 중심화 후 자기일치 셀평균 / 쌍평균 셀평균
        Pr, Pc = consensus_of(raw), consensus_of(cen)
        nm = sorted(Pc)
        pre = [float(np.dot(Pr[a], Pr[b]))
               for a, b in itertools.combinations(nm, 2)]
        post = [float(np.dot(Pc[a], Pc[b]))
                for a, b in itertools.combinations(nm, 2)]
        wi = []
        for a in nm:
            cs = sorted(cen[a])
            wi += [float(np.dot(cen[a][i], cen[a][j]))
                   for i, j in itertools.combinations(cs, 2)]
        cellstat = dict(pre_centring_pair_mean=float(np.mean(pre)),
                        post_centring_pair_mean=float(np.mean(post)),
                        post_centring_within_mean=float(np.mean(wi)))

        # 귀무를 두 표현에서 각각 구한다. 중심화 전 표현의 귀무가
        # 핵심이다 - 원천 논문 분류기가 보는 표현에서 화합물 정체성이
        # 라벨을 뒤섞은 것보다 나은지를 그 수가 답한다.
        nulls = {}
        for label, pat, obs in (('uncentred', raw, hu), ('centred', cen, hc)):
            nsp, do_exact = exact_space(pat)
            if do_exact:
                nn = null_counts(pat, exact=True)
                nulls[label] = dict(
                    mode='exact', n_relabellings=nsp,
                    mean=float(np.mean(nn)),
                    p5=float(np.percentile(nn, 5)),
                    p95=float(np.percentile(nn, 95)),
                    p_le_observed=float(np.mean(np.asarray(nn) <= obs)))
            else:
                per_seed, floors = [], []
                for s in SEEDS:
                    nn = null_counts(pat, rng=np.random.default_rng(s))
                    arr = np.asarray(nn)
                    per_seed.append(dict(
                        seed=s, mean=float(np.mean(arr)),
                        p_le_observed=float(np.mean(arr <= obs))))
                    # 도달 가능한 최소 p. 역전수는 0..n 의 이산값만 갖는다.
                    # 관측이 0 이었다 해도 p 는 P(귀무 = 0) 아래로 못 간다.
                    floors.append(float(np.mean(arr <= 0)))
                nulls[label] = dict(
                    mode='sampled', n_relabellings=nsp, per_seed=per_seed,
                    mean=float(np.mean([p['mean'] for p in per_seed])),
                    p_le_observed_min=float(
                        min(p['p_le_observed'] for p in per_seed)),
                    p_le_observed_max=float(
                        max(p['p_le_observed'] for p in per_seed)),
                    p_floor_max=float(max(floors)))
        null_block = nulls['centred']

        print(f'\n  {tag} {mode}  사례 {nu}건')
        print(f'    중심화 전 역전 {hu}/{nu}   귀무 평균 '
              f'{nulls["uncentred"]["mean"]:.2f}/{nu}  '
              f'({nulls["uncentred"]["mode"]})')
        print(f'    중심화 후 역전 {hc}/{nc}   귀무 평균 '
              f'{nulls["centred"]["mean"]:.2f}/{nc}  '
              f'({nulls["centred"]["mode"]})')
        for d in du:
            if d['reversed']:
                print(f'      [전] {d["compound"]:<6} 자기 {d["within"]:+.3f} '
                      f'< {d["nearest_other"]} {d["nearest_cos"]:+.3f}')

        report['cells'].append(dict(
            dataset=tag, mode=mode, n_instances=nu,
            reversals_uncentred=hu, reversals_centred=hc,
            detail_uncentred=du, detail_centred=dc,
            null_uncentred=nulls['uncentred'], null_centred=null_block,
            **cellstat))

    # 6셀 BH 보정. 표집 셀은 시드 간 최대 p 를 쓴다 (보수적).
    def _p(block):
        if block['mode'] == 'exact':
            return block['p_le_observed']
        return block['p_le_observed_max']

    def _bh(pv):
        pv = np.sort(np.asarray(pv, float))
        k = len(pv)
        adj = np.minimum.accumulate(
            (pv * k / np.arange(1, k + 1))[::-1])[::-1]
        return [float(x) for x in adj], float(adj[0])

    bh = {}
    for label in ('uncentred', 'centred'):
        rows = [(c['dataset'], c['mode'], _p(c['null_' + label]),
                 c['null_' + label].get('p_floor_max'))
                for c in report['cells'] if 'unavailable' not in c]
        rows.sort(key=lambda r: r[2])
        adj, mn = _bh([r[2] for r in rows])
        bh[label] = dict(
            n_cells=len(rows), min_q=mn,
            cells=[dict(dataset=d, mode=md, p=p, q=q, p_floor=fl)
                   for (d, md, p, fl), q in zip(rows, adj)],
            n_cells_q_below_05=int(sum(1 for q in adj if q < 0.05)),
            n_cells_floor_above_05=int(
                sum(1 for r in rows if r[3] is not None and r[3] >= 0.05)))
    report['bh_over_cells'] = bh
    print('\n  [BH 보정, 6셀]')
    for label in ('uncentred', 'centred'):
        b = bh[label]
        print(f'    중심화 {"전" if label == "uncentred" else "후"}  '
              f'최소 q {b["min_q"]:.4f}   q<0.05 인 셀 '
              f'{b["n_cells_q_below_05"]}/{b["n_cells"]}')
        for r in b['cells']:
            fl = r.get('p_floor')
            fs = '[수집 불가]' if fl is None else f'{fl:.4f}'
            reach = '' if fl is None else (
                '  도달가능' if fl < 0.05 else '  **0.05 도달 불가**')
            print(f'      {r["dataset"]} {r["mode"]:<7} '
                  f'p {r["p"]:.4f}  q {r["q"]:.4f}  최소달성p {fs}{reach}')

    ic = inter_channel()
    report['inter_channel'] = ic
    print(f'\n  [채널 간] 같은 플레이트 두 채널, {ic["n_pairs"]}건')
    print(f'    중심화 후 평균 {ic["centred_mean"]:+.4f}   '
          f'중심화 전 평균 {ic["uncentred_mean"]:+.4f}')

    pb = pair_baseline()
    report['pair_baseline'] = pb
    print('\n  [셀 기저선] 작용제끼리의 쌍 대 비작용제의 대작용제 쌍')
    for k, v in pb.items():
        if 'unavailable' in v:
            print(f'    {k:<12} [수집 불가: {v["unavailable"]}]')
            continue
        if 'note' in v:
            print(f'    {k:<12} {v["note"]} (작용제 {v["n_agents"]}종)')
            continue
        print(f'    {k:<12} 작용제 {v["n_agent_pairs"]}쌍 평균 '
              f'{v["agent_pair_mean"]:+.4f} 중위 {v["agent_pair_median"]:+.4f}'
              f'  최저쌍 {"-".join(v["most_negative_pair"])} '
              f'{v["most_negative_value"]:+.4f}')
        for s2, d2 in v['non_agents'].items():
            print(f'      {s2:<8} 평균 {d2["mean"]:+.4f} '
                  f'중위 {d2["median"]:+.4f}  작용제쌍 '
                  f'{d2["n_agent_pairs_below_its_mean"]}/{v["n_agent_pairs"]}'
                  f' 이 더 낮다  기저선 위 {d2["above_agent_pair_mean"]}')

    cf = collapse_free_counts()
    report['collapse_free'] = cf
    cfu = sum(v.get('reversals_uncentred', 0) for v in cf.values())
    cfc = sum(v.get('reversals_centred', 0) for v in cf.values())
    cfn = sum(v.get('n_instances', 0) for v in cf.values())
    report['collapse_free_totals'] = dict(n_instances=cfn,
                                          reversals_uncentred=cfu,
                                          reversals_centred=cfc)

    ok = [c for c in report['cells'] if 'unavailable' not in c]
    report['totals'] = dict(
        n_instances=inst,
        reversals_uncentred=tot['uncentred'],
        reversals_centred=tot['centred'],
        null_mean_sum_uncentred=float(
            sum(c['null_uncentred']['mean'] for c in ok)),
        null_mean_sum_centred=float(
            sum(c['null_centred']['mean'] for c in ok)),
        pre_centring_pair_mean_range=[
            float(min(c['pre_centring_pair_mean'] for c in ok)),
            float(max(c['pre_centring_pair_mean'] for c in ok))],
        post_centring_pair_mean_range=[
            float(min(c['post_centring_pair_mean'] for c in ok)),
            float(max(c['post_centring_pair_mean'] for c in ok))],
        post_centring_within_mean_range=[
            float(min(c['post_centring_within_mean'] for c in ok)),
            float(max(c['post_centring_within_mean'] for c in ok))])
    print('\n' + '-' * 78)
    print(f'  합계 사례 {inst}건')
    print(f'    중심화 전 역전 {tot["uncentred"]}/{inst}')
    print(f'    중심화 후 역전 {tot["centred"]}/{inst}')
    print('  중심화는 역전을 줄인다 (늘리지 않는다). 이 통계는 중심화의'
          ' 산물이 아니다.')
    print(f'\n  [자유도] 중복 채널 축약을 끈 경우  사례 {cfn}건')
    print(f'    중심화 전 역전 {cfu}/{cfn}   중심화 후 역전 {cfc}/{cfn}')
    for k, v in cf.items():
        if 'unavailable' in v:
            print(f'    {k:<12} [수집 불가: {v["unavailable"]}]')
            continue
        print(f'    {k:<12} 채널 {v["n_channels"]:>2}  '
              f'전 {v["reversals_uncentred"]}/{v["n_instances"]}  '
              f'후 {v["reversals_centred"]}/{v["n_instances"]}')

    with open('outputs/reversal_audit.json', 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print('\n[done] outputs/reversal_audit.json')


if __name__ == '__main__':
    main()
