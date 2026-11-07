#!/usr/bin/env python3
"""희소 지시자 인코딩이 거리에 넣는 받침대 (2026-11-03 신설, T1-14).

왜 필요한가. GATE 8 7차에서 리뷰어 R4 가 서술자 인코딩이 결과를 움직이는데
본문의 analyst freedom 목록에 빠져 있다고 지적했다. 그 지적은 구판 상관필터
위에서 측정된 것이므로(타분 1.995 -> 1.444, 매칭 4 -> 6), 결정적 필터에서
다시 센다.

무엇이 문제인가. 병합공간 38좌표에는 28종 중 정확히 1종만 1 이고 나머지가
0 인 이진 지시자가 넷 있다. n_As(루이사이트), n_P_CN(타분), n_P_S(VX),
n_nitro(클로로피크린) 이다. 1종만 1 인 0/1 열은 평균이 1/n 이고
분산이 정확히 1/n 이므로(표본분산, ddof=1), 표준편차가 1/sqrt(n) 이고
1 을 가진 화합물과 0 을 가진 화합물의 z 차는 정확히 sqrt(n) 이다. 따라서 그
한 화합물에서 다른 어느 화합물로 가는 제곱거리에 **정확히 n = 28** 이
고정으로 더해진다. 두 0 사이에는 0 이 더해진다. n=28 에서 z 는 각각
+5.103 과 -0.189 이고 차의 제곱이 28.0 이다.

그러므로 단일값 지시자는 그것을 가진 한 화합물을 나머지 전부에서 일정하게
밀어낸다. 화학적 근접성과 무관한 받침대다. 타분은 n_P_CN 을 가진 유일한
화합물이므로 DMMP 와의 거리가 이 받침대를 안고 있고, 역치 2 바로 위에 앉은
그 거리가 받침대에서 얼마를 받는지 세는 것이 이 스크립트의 일이다.

받침대가 '틀렸다'고 주장하지 않는다. 지시자를 넣는 것도 빼는 것도 분석자의
선택이며, 어느 쪽이 옳은지는 이 데이터가 정하지 않는다. 그것이 바로 요점이다.

사용: python3 sparse_indicators.py
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, pick

THRESHOLD = 2.0
N_PERM_PANEL = 17          # C(17,3) = 680, 28종 패널의 전수 열거 모집단


def merged_space(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return [c for c in num if c != OMEGA_DERIVED]


def classify(df, cols):
    """열마다 (고유값 수, 비0 개수) 를 낸다. 지시자 판정의 근거다."""
    out = {}
    for c in cols:
        v = df[c].to_numpy(float)
        out[c] = dict(n_unique=int(len(np.unique(v))),
                      n_nonzero=int((v != 0).sum()),
                      vmin=float(v.min()), vmax=float(v.max()),
                      is_binary=bool(set(np.unique(v)) <= {0.0, 1.0}))
    return out


def zmat(df, keep):
    X = df[keep].to_numpy(float)
    return (X - X.mean(0)) / X.std(0, ddof=1)


def panel(df, keep, sim_cls=('simulant',)):
    """match_by_panel 과 같은 거리·매칭 규칙 (넓은 판독)."""
    Z = zmat(df, keep)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    kl, name = df.compound_class.tolist(), df.compound_name.tolist()
    i_ag = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_si = [i for i, k in enumerate(kl) if k in sim_cls]
    i_cm = [i for i, k in enumerate(kl)
            if k != 'real_agent' and k not in sim_cls]
    res, matched = {}, []
    for i in i_ag:
        js = min(i_si, key=lambda j: D[i, j])
        dn = float(D[i, js] / nn_med)
        ok = bool(dn <= THRESHOLD and all(D[i, js] < D[i, j] for j in i_cm))
        res[name[i]] = dict(nearest=name[js], d_norm=round(dn, 4), match=ok)
        if ok:
            matched.append(name[i])
    return dict(n_dims=len(keep), nn_median=round(nn_med, 4),
                n_match=len(matched), matched=sorted(matched), agents=res)


def perm_p(df, keep):
    """지정 3종의 평균차에 대한 단측 p. C(17,3)=680 전수 열거.

    analysis_main 과 같은 통계량(지정 모사체 평균거리 - 비지정 평균거리)을
    쓴다. 모집단은 작용제가 아닌 17종에서 3종을 뽑는 모든 조합이다.
    """
    Z = zmat(df, keep)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    kl = df.compound_class.tolist()
    i_ag = [i for i, k in enumerate(kl) if k == 'real_agent']
    pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    obs_set = tuple(i for i, k in enumerate(kl) if k == 'simulant')

    def stat(sel):
        sel = set(sel)
        rest = [j for j in pool if j not in sel]
        a = np.mean([min(D[i, j] for j in sel) for i in i_ag])
        b = np.mean([min(D[i, j] for j in rest) for i in i_ag])
        return float(a - b)

    obs = stat(obs_set)
    allc = list(itertools.combinations(pool, 3))
    vals = np.array([stat(c) for c in allc])
    return dict(n_perm=len(allc), observed=obs,
                p_one_sided=float((vals <= obs).sum() / len(allc)))


def pedestal_share(df, keep, sparse, pairs):
    """쌍마다 제곱거리에서 희소 지시자 좌표가 차지하는 비율."""
    Z = zmat(df, keep)
    idx = {n: i for i, n in enumerate(df.compound_name)}
    jj = [keep.index(c) for c in sparse if c in keep]
    out = {}
    for a, b in pairs:
        if a not in idx or b not in idx:
            continue
        d = Z[idx[a]] - Z[idx[b]]
        tot = float((d ** 2).sum())
        part = float((d[jj] ** 2).sum()) if jj else 0.0
        out[f'{a} -> {b}'] = dict(sq_total=round(tot, 3),
                                  sq_sparse=round(part, 3),
                                  share=round(part / tot, 4) if tot else None)
    return out


def mantel_stability(df, keep0, variants):
    """인코딩 변종에서 Mantel 12셀 BH 가 움직이는가. 논지 자체의 안정성.

    왜 이것까지 재는가. 희소 지시자를 빼면 28종 패널의 매칭 수와 순열 p 가
    크게 움직인다. 그러면 논문의 본 주장(서술자 거리가 응답 거리를 예측하지
    못한다)도 같이 움직이는지가 당연한 질문이 된다. 그쪽은 Mantel 12셀의
    BH 보정값이 결정하므로 변종마다 다시 계산한다. 움직이지 않는다면 어느
    결과가 논문을 떠받치고 어느 결과가 분석자 선택의 함수인지를 독자가
    구별할 수 있다.
    """
    import deltaE_premise as m
    from scipy import stats
    from spaces import ELECTRONIC_28

    ELEC3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]
    idx = {nm: i for i, nm in enumerate(df.compound_name)}
    out = []
    for label, drop in variants:
        keep = [c for c in keep0 if c not in set(drop)]
        st = [c for c in keep if c not in ELECTRONIC_28]
        el = [c for c in keep if c in ELEC3]
        SP = {'structural': st, 'electronic3': el, 'merged': st + el}
        DIST = {}
        for k, v in SP.items():
            Z = zmat(df, v)
            DIST[k] = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
        ps = []
        for tag in ['d1', 'd2', 'd3']:
            for mode in ['uv', 'normal']:
                D, dy = m.load(tag, mode)
                P = m.consensus(m.centred(m.raw_patterns(D, dy)))
                ags = [a for a in sorted(P) if a in m.MAP]
                if len(ags) < 4:
                    continue
                k = len(ags)
                pr = list(itertools.combinations(range(k), 2))
                dis = np.array([1 - float(np.dot(P[ags[i]], P[ags[j]]))
                                for i, j in pr])
                for sp in SP:
                    M = np.array([[DIST[sp][idx[m.MAP[a]], idx[m.MAP[b]]]
                                   for b in ags] for a in ags])
                    d = np.array([M[i, j] for i, j in pr])
                    r = float(stats.spearmanr(d, dis).statistic)
                    cnt, perms = 0, list(itertools.permutations(range(k)))
                    for pm in perms:
                        dp = np.array([M[pm[i], pm[j]] for i, j in pr])
                        rp = stats.spearmanr(dp, dis).statistic
                        if rp == rp and abs(rp) >= abs(r):
                            cnt += 1
                    ps.append(cnt / len(perms))
        o = np.sort(np.array(ps))
        K = len(o)
        bh = np.minimum.accumulate((o * K / np.arange(1, K + 1))[::-1])[::-1]
        out.append(dict(variant=label, n_dropped=len(drop), n_dims=len(keep),
                        n_cells=int(K), min_p=float(o[0]),
                        min_q=float(bh[0]), n_sig_q=int((bh < 0.05).sum())))
    return out


def main():
    df = pd.read_csv('descriptor_cache_xtb.csv')
    base = merged_space(df)
    keep0, _ = pick(df, base, PROTECTED)
    info = classify(df, keep0)

    print('=' * 78)
    print('[T1-14] 희소 지시자 인코딩 민감도 — 결정적 필터 위에서')
    print('=' * 78)
    print(f'\n기준 병합공간 {len(keep0)}좌표. 이진 지시자의 비0 개수 분포:')
    binz = sorted((v['n_nonzero'], c) for c, v in info.items()
                  if v['is_binary'])
    for nz, c in binz:
        print(f'   {c:<22} 비0 {nz}/28')

    singles = sorted(c for c, v in info.items()
                     if v['is_binary'] and v['n_nonzero'] == 1)
    print(f'\n단일값 지시자 {len(singles)}개: {", ".join(singles)}')
    carrier = {c: df.compound_name[int(np.argmax(df[c].to_numpy(float)))]
               for c in singles}
    for c, who in carrier.items():
        print(f'   {c:<14} 보유 화합물 {who}')

    # 받침대 크기를 식과 실측 양쪽으로 확인한다. 1종만 1 인 0/1 열의
    # 표본분산은 정확히 1/n 이므로 z 차는 sqrt(n), 제곱 기여는 n 이다.
    n = len(df)
    _c = singles[0]
    _v = df[_c].to_numpy(float)
    _z = (_v - _v.mean()) / _v.std(ddof=1)
    _emp = float((_z.max() - _z.min()) ** 2)
    assert abs(_emp - n) < 1e-9, (_emp, n)
    print(f'\n단일값 지시자 하나가 보유 화합물과 나머지 사이의 제곱거리에 '
          f'더하는 양 = n = {n}  (식 sqrt(n)^2, 실측 {_emp:.6f})')

    # 민감도 변종. 비0 개수 역치를 올려가며 지시자를 떨어뜨린다.
    variants = [('baseline', [])]
    for k in (1, 2, 3, 4, 5):
        drop = sorted(c for c, v in info.items()
                      if v['is_binary'] and v['n_nonzero'] <= k)
        variants.append((f'drop_binary_nonzero_le_{k}', drop))
    low_card = sorted(c for c, v in info.items() if v['n_unique'] <= 5)
    variants.append(('drop_cardinality_le_5', low_card))

    rows = []
    print('\n민감도 (넓은 판독, 28종 패널)')
    print(f'{"변종":<28}{"제거":>4}{"차원":>5}{"매칭":>5}  타분 d_norm   단측 p')
    for label, drop in variants:
        keep = [c for c in keep0 if c not in set(drop)]
        if len(keep) < 5:
            continue
        r = panel(df, keep)
        pp = perm_p(df, keep)
        tb = r['agents']['tabun (GA)']['d_norm']
        rows.append(dict(variant=label, dropped=drop, n_dropped=len(drop),
                         **{k: v for k, v in r.items() if k != 'agents'},
                         tabun_d_norm=tb,
                         tabun_match=r['agents']['tabun (GA)']['match'],
                         agents=r['agents'],
                         p_one_sided=pp['p_one_sided'],
                         observed=pp['observed']))
        print(f'{label:<28}{len(drop):>4}{r["n_dims"]:>5}{r["n_match"]:>5}'
              f'   {tb:<10}  {pp["p_one_sided"]:.4f}   '
              f'{"".join(sorted(a[0] for a in r["matched"]))}')

    pairs = [('tabun (GA)', 'DMMP'), ('sarin (GB)', 'DMMP'),
             ('soman (GD)', 'DMMP'), ('cyclosarin (GF)', 'DMMP'),
             ('VX', 'DPGME'), ('lewisite 1 (L)', 'dichloromethane')]
    share_singles = pedestal_share(df, keep0, singles, pairs)
    all_bin = sorted(c for c, v in info.items() if v['is_binary'])
    share_allbin = pedestal_share(df, keep0, all_bin, pairs)

    print('\n제곱거리에서 지시자가 차지하는 비율 (기준 38좌표)')
    print(f'{"쌍":<34}{"제곱거리":>10}{"단일값4":>9}{"이진전체":>9}')
    for k in share_singles:
        a, b = share_singles[k], share_allbin[k]
        print(f'{k:<34}{a["sq_total"]:>10.2f}'
              f'{100*a["share"]:>8.1f}%{100*b["share"]:>8.1f}%')

    # 논지 안정성.
    # 2026-11-07 (10차 R3). 종전에는 '전수 열거가 비싸므로' 세 변종만 돌았고,
    # 그 셋에 매칭을 가장 크게 움직이는 drop_cardinality_le_5 (매칭 6건) 가
    # 빠져 있었다. 그런데 원고는 그 상태로 "survives every choice we
    # enumerated" 라고 적었다. 즉 가장 큰 자유도에서 Mantel 을 돌린 적이
    # 없는데 전부 견딘다고 주장하고 있었다. Mantel 은 6!=720 전수이고 12셀
    # 이라 변종 하나당 8640 회, 비싸지 않다. 전 변종을 돈다.
    mv = list(variants)
    ms = mantel_stability(df, keep0, mv)
    print('\n논지 안정성 — Mantel 12셀 BH (서술자 거리 대 응답 거리)')
    print(f'{"변종":<28}{"차원":>5}{"최소 p":>9}{"최소 q":>9}{"q<0.05":>8}')
    for r in ms:
        print(f'{r["variant"]:<28}{r["n_dims"]:>5}{r["min_p"]:>9.4f}'
              f'{r["min_q"]:>9.4f}{r["n_sig_q"]:>8}건')

    json.dump(dict(n_compounds=n, baseline_dims=len(keep0),
                   mantel_stability=ms,
                   column_info=info, singleton_indicators=singles,
                   singleton_carriers=carrier,
                   pedestal_per_singleton=n,
                   pedestal_formula='sqrt(n)^2 = n; 1종만 1 인 '
                                    '0/1 열의 표본분산은 1/n',
                   variants=rows,
                   share_singletons=share_singles,
                   share_all_binary=share_allbin),
              open('outputs/sparse_indicators.json', 'w'), indent=1,
              ensure_ascii=False)
    print('\n-> outputs/sparse_indicators.json')


if __name__ == '__main__':
    main()
