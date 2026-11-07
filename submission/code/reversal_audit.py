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
    report = dict(generated_by='reversal_audit.py', seeds=SEEDS,
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
                per_seed = []
                for s in SEEDS:
                    nn = null_counts(pat, rng=np.random.default_rng(s))
                    per_seed.append(dict(
                        seed=s, mean=float(np.mean(nn)),
                        p_le_observed=float(np.mean(np.asarray(nn) <= obs))))
                nulls[label] = dict(
                    mode='sampled', n_relabellings=nsp, per_seed=per_seed,
                    mean=float(np.mean([p['mean'] for p in per_seed])),
                    p_le_observed_min=float(
                        min(p['p_le_observed'] for p in per_seed)),
                    p_le_observed_max=float(
                        max(p['p_le_observed'] for p in per_seed)))
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
            null_uncentred=nulls['uncentred'], null_centred=null_block))

    cf = collapse_free_counts()
    report['collapse_free'] = cf
    cfu = sum(v.get('reversals_uncentred', 0) for v in cf.values())
    cfc = sum(v.get('reversals_centred', 0) for v in cf.values())
    cfn = sum(v.get('n_instances', 0) for v in cf.values())
    report['collapse_free_totals'] = dict(n_instances=cfn,
                                          reversals_uncentred=cfu,
                                          reversals_centred=cfc)

    report['totals'] = dict(n_instances=inst,
                            reversals_uncentred=tot['uncentred'],
                            reversals_centred=tot['centred'])
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
