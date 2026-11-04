#!/usr/bin/env python3
"""ΔE 응답 크기와 농도 추세의 SNR 교란 검사 (2026-11-02 신설, T2-5).

왜 이 모양인가. GATE 8 7차는 두 가지를 요구했다. 하나는 ΔE 잡음 바닥을
보고하라는 것, 다른 하나는 농도가 내려가며 DMMP-사린 유사도가 떨어지는
추세가 화학이 아니라 신호대잡음비의 함수일 가능성을 다루라는 것이다.

첫째는 이 저장소의 데이터로 할 수 없다. 예치된 피벗은 (물질, 농도) 한 행이고
반복측정은 이미 평균되어 있다. 0 µM 행은 ΔE 가 공시료 기준이라 정의상
0 이다 (d3 uv 의 최대 8.03e-02 만 예외이며 이는 반올림 잔차다). 반복측정
수준 값도 공시료 대 공시료 변동도 예치물 안에 없으므로, 잡음 바닥은
[수집 불가] 로 적고 무엇이 있어야 하는지를 명시한다. 가상 수치를 만들지
않는다 ([C-3]).

둘째는 할 수 있다. 농도 대신 실측 응답 크기를 설명변수로 놓으면 된다.
유사도가 농도가 아니라 응답 크기를 따라간다면 그 추세는 SNR 효과다. 농도와
크기는 단조로 묶여 있으므로 두 설명을 완전히 분리할 수는 없고, 이 스크립트는
분리했다고 주장하지 않는다. 다만 어느 쪽과 더 붙어 있는지는 센다.

함께 내는 것: 염료별 응답 크기 분포와, 약한 염료를 떨어뜨렸을 때 전제 검정과
DMMP 결론이 버티는지의 민감도.

사용: python3 deltaE_noise.py
"""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

import deltaE_premise as m

CONC_NZ = [c for c in m.CONC if c != 0]


def cells():
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            yield tag, mode


def blank_audit():
    """0 µM 행이 정보를 가지고 있는지 확인한다. 가지고 있지 않음을 기록한다."""
    rows = []
    for tag, mode in cells():
        D = pd.read_csv(f'deltaE/{tag}.{mode}_deltaE_pivot.csv')
        dy = [c for c in D.columns if c.lower().startswith('dye')]
        z = D[D.Concentration == 0][dy].to_numpy(float)
        rows.append(dict(cell=f'{tag}.{mode}', n_blank_rows=int(len(z)),
                         max_abs=float(np.abs(z).max()) if z.size else None,
                         mean_abs=float(np.abs(z).mean()) if z.size else None,
                         carries_variation=bool(z.size and np.abs(z).max() > 1e-6)))
    return rows


def dye_magnitudes():
    """염료별 |ΔE| 분포. 유효 채널(중복 축약 후)만 본다."""
    rows = []
    for tag, mode in cells():
        D, dy = m.load(tag, mode)
        A = D[D.Concentration.isin(CONC_NZ)]
        X = np.abs(A[dy].to_numpy(float))
        med = X.mean(0)
        tot = float(med.sum())
        order = np.argsort(med)
        for rank, j in enumerate(order):
            rows.append(dict(cell=f'{tag}.{mode}', dye=dy[j],
                             rank_weakest=rank + 1,
                             mean_abs_deltaE=float(med[j]),
                             share_of_cell=float(med[j] / tot) if tot else None))
    return rows


def drop_weak_sensitivity(k_list=(0, 3, 6, 9)):
    """약한 염료를 k 개 떨어뜨리며 전제 검정과 DMMP 결론을 다시 본다.

    전제 검정은 deltaE_premise.discrimination 과 같은 화합물 단위 부호뒤집기
    전수열거다. 여기서 다시 구현하지 않고 같은 규칙을 적용한다.
    """
    rows = []
    for tag, mode in cells():
        D, dy = m.load(tag, mode)
        A = D[D.Concentration.isin(CONC_NZ)]
        mean_abs = np.abs(A[dy].to_numpy(float)).mean(0)
        order = [dy[j] for j in np.argsort(mean_abs)]       # 약한 것부터
        for k in k_list:
            if k >= len(dy) - 3:
                continue
            keep = [c for c in dy if c not in set(order[:k])]
            cen = m.centred(m.raw_patterns(D, keep))
            if not cen:
                continue
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
            obs = float(d.mean())
            cnt = sum(1 for s in itertools.product([1, -1], repeat=n)
                      if float((np.array(s) * d).mean()) >= obs - 1e-15)
            row = dict(cell=f'{tag}.{mode}', n_dropped=k, n_channels=len(keep),
                       n_compounds=n, mean_diff=obs, p_exact=cnt / 2 ** n,
                       p_min_attainable=1 / 2 ** n)
            # DMMP 가 있는 셀에서는 헤드라인 비교도 같이 본다.
            if 'DMMP' in P:
                s = {a: float(np.dot(P['DMMP'], P[a])) for a in P if a != 'DMMP'}
                g = [s[a] for a in ['GB', 'GD', 'GF'] if a in s]
                if g and 'VX' in s:
                    row['dmmp_G_mean'] = float(np.mean(g))
                    row['dmmp_VX'] = float(s['VX'])
                    row['dmmp_G_minus_VX'] = float(np.mean(g) - s['VX'])
            rows.append(row)
    return rows


def _exact_spearman_p(x, y):
    """Spearman 의 양측 정확 p. n! 전수열거.

    왜 필요한가. 농도점은 네 개 또는 다섯 개다. scipy 의 spearmanr 은 t 근사
    p 를 주는데, n=4 에서 완전 단조이면 그 값이 0.000 으로 인쇄된다. 전수로
    세면 양측 최소 p 는 2/4! = 0.0833 이고 0 이 될 수 없다. 근사값을 본문에
    적으면 없는 유의성을 주장하게 된다.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    r = float(stats.spearmanr(x, y).statistic)
    if not (r == r):
        return float('nan')
    cnt = 0
    for pm in itertools.permutations(range(n)):
        rp = stats.spearmanr(x, y[list(pm)]).statistic
        if rp == rp and abs(rp) >= abs(r) - 1e-12:
            cnt += 1
    import math
    return cnt / math.factorial(n)


def concentration_vs_magnitude():
    """농도별 응답 크기와 DMMP-작용제 유사도. 추세가 어느 쪽을 따르는가."""
    rows = []
    for tag, mode in cells():
        D, dy = m.load(tag, mode)
        raw = m.raw_patterns(D, dy)             # 단위정규화된 패턴
        if 'DMMP' not in raw:
            continue
        # 크기는 정규화 전 원시 벡터에서 재다.
        norms = {}
        for c in CONC_NZ:
            sub = D[D.Concentration == c]
            if len(sub):
                V = sub[dy].to_numpy(float)
                norms[c] = float(np.linalg.norm(V, axis=1).mean())
        # 농도별 중심화. discrimination 과 같은 경로를 쓴다.
        cen = m.centred(raw)
        for target in ['GB', 'GD', 'GF', 'GA', 'VX']:
            if target not in cen:
                continue
            cs = sorted(set(cen['DMMP']) & set(cen[target]) & set(norms))
            if len(cs) < 3:
                continue
            cos = [float(np.dot(cen['DMMP'][c], cen[target][c])) for c in cs]
            nm = [norms[c] for c in cs]
            r_conc = float(stats.spearmanr(cs, cos).statistic)
            r_norm = float(stats.spearmanr(nm, cos).statistic)
            rows.append(dict(cell=f'{tag}.{mode}', target=target,
                             n_points=len(cs),
                             concentrations=[int(c) for c in cs],
                             mean_response_norm=[round(x, 3) for x in nm],
                             cosine=[round(x, 4) for x in cos],
                             rho_vs_concentration=r_conc,
                             p_vs_concentration=_exact_spearman_p(cs, cos),
                             p_vs_concentration_asymptotic=float(
                                 stats.spearmanr(cs, cos).pvalue),
                             rho_vs_response_norm=r_norm,
                             p_vs_response_norm=_exact_spearman_p(nm, cos),
                             rho_conc_vs_norm=float(
                                 stats.spearmanr(cs, nm).statistic),
                             separable=bool(abs(float(
                                 stats.spearmanr(cs, nm).statistic)) < 1.0)))
    return rows


def main():
    print('=' * 78)
    print('[T2-5] ΔE 응답 크기, 잡음 바닥 가용성, 농도 추세의 SNR 교란')
    print('=' * 78)

    ba = blank_audit()
    print('\n(1) 0 µM 행 감사 — 잡음 바닥을 여기서 뽑을 수 있는가')
    for r in ba:
        print(f'  {r["cell"]:<12} 행 {r["n_blank_rows"]}  '
              f'|값| 최대 {r["max_abs"]:.4g}  변동 보유 {r["carries_variation"]}')
    print('  -> [수집 불가: ΔE 잡음 바닥]. 예치 피벗은 반복측정 평균이고 '
          '0 µM 행은 정의상 0 이다.')
    print('     필요한 것: 반복측정 수준 ΔE 또는 동일 조건 공시료 반복 촬상.')

    dm = dye_magnitudes()
    print('\n(2) 셀별 가장 약한 염료 다섯 (평균 |ΔE|)')
    for tag, mode in cells():
        cl = f'{tag}.{mode}'
        sub = [r for r in dm if r['cell'] == cl][:5]
        s = ', '.join(f'{r["dye"]} {r["mean_abs_deltaE"]:.3f}'
                      f'({100*r["share_of_cell"]:.1f}%)' for r in sub)
        print(f'  {cl:<12} {s}')

    ds = drop_weak_sensitivity()
    print('\n(3) 약한 염료 제거 민감도 — 전제 검정과 DMMP 결론')
    for r in ds:
        extra = ''
        if 'dmmp_G_minus_VX' in r:
            extra = (f'   DMMP: G {r["dmmp_G_mean"]:+.3f} vs VX '
                     f'{r["dmmp_VX"]:+.3f}  차 {r["dmmp_G_minus_VX"]:+.3f}')
        print(f'  {r["cell"]:<12} 제거 {r["n_dropped"]:<2} 채널 '
              f'{r["n_channels"]:<3} 차 {r["mean_diff"]:+.3f}  '
              f'정확p {r["p_exact"]:.4f}{extra}')

    cv = concentration_vs_magnitude()
    print('\n(4) 농도 대 응답크기 — 유사도는 어느 쪽을 따르는가')
    for r in cv:
        print(f'  {r["cell"]:<12} DMMP-{r["target"]:<3} n={r["n_points"]}  '
              f'rho(농도) {r["rho_vs_concentration"]:+.3f} '
              f'(정확p={r["p_vs_concentration"]:.4f}, '
              f'근사p={r["p_vs_concentration_asymptotic"]:.4f})   '
              f'rho(응답크기) {r["rho_vs_response_norm"]:+.3f}   '
              f'rho(농도,크기) {r["rho_conc_vs_norm"]:+.3f}  '
              f'분리가능 {r["separable"]}')

    json.dump(dict(noise_floor='[수집 불가: 반복측정 수준 ΔE 및 공시료 반복 '
                               '촬상이 예치물에 없다. 0 µM 행은 정의상 0.]',
                   blank_audit=ba, dye_magnitudes=dm,
                   drop_weak_sensitivity=ds,
                   concentration_vs_magnitude=cv),
              open('outputs/deltaE_noise.json', 'w'), indent=1,
              ensure_ascii=False)
    print('\n-> outputs/deltaE_noise.json')


if __name__ == '__main__':
    main()
