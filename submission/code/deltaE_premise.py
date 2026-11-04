#!/usr/bin/env python3
"""실측 비색 어레이 응답으로 논문의 중심 전제를 검정한다 (2026-10-12 신설).

배경. GATE 8 4차에서 심사자 두 명이 독립적으로 같은 지적을 했다. 원고는
"전자 서술자 공간이 기상 센서 응답을 지배한다"에 제목과 초록을 걸면서
센서 응답 자료를 한 건도 쓰지 않는다. 그 상태에서 "물리화학적 근거로
검증된 모사체는 탐지 화학에 대해 아무 타당성도 함의하지 않는다"고 결론
내리는 것은 전제의 재서술이다. 검정에 필요한 자료가 이미 있다.

자료. 29염료 비색 어레이의 ΔE 응답.
  d1: DMMP, GA, GB, GD, GF, VX   (신경작용제 + 모사체)
  d2: AC, CG, CK, HD, HN, L, PS  (혈액·질식·수포작용제)
  d3: A, B, C, D                 (A 계열)
각 물질 × 29염료 × 4농도 × 2촬영모드(uv 여기 / normal 반사).

설계 원칙.
  1) d1/d2/d3 는 같은 물질의 반복이 아니라 서로 다른 물질군이다. 서로 다른
     플레이트의 ΔE 를 직접 비교하려면 가환성 가정이 필요하므로 분석은
     데이터셋 안에서만 하고 결과를 따로 보고한다.
  2) uv 와 normal 은 반복이 아니라 물리적으로 다른 채널이다(여기 형광 대
     반사). 같은 물질에서도 두 패턴의 코사인이 평균 0.58 에 그친다. 따라서
     독립 모달리티로 각각 분석하고, 두 채널에서 같은 결론이 나오는지를
     재현의 근거로 쓴다.
  3) 크기가 아니라 패턴을 본다. 각 농도의 29차원 벡터를 단위길이로 정규화한다.
  4) **염료별 중심화가 필수다.** ΔE 는 전부 음이 아니므로 모든 물질이
     "무언가 반응했다"는 공통 성분을 공유하고, 정규화만 하면 코사인이 1 쪽으로
     압축된다. 실제로 중심화 전에는 물질간 유사도 0.77~0.80 이 물질내 잡음
     0.79~0.83 과 구별되지 않았다(p=0.06~0.94). 각 농도에서 물질 축 평균을
     빼면 물질내 +0.35~+0.54 대 물질간 -0.10~-0.33 으로 갈라진다
     (여섯 조합 전부 p<0.03). 공통 성분 제거는 비색 어레이 분석의 표준이다.
"""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

CONC = [10, 50, 100, 500]
SEED, N_BOOT = 42, 10000

# 실측 라벨 -> 28종 패널 이름. HN(질소겨자)과 d3 의 A~D 는 패널에 없다.
MAP = {'DMMP': 'DMMP', 'GA': 'tabun (GA)', 'GB': 'sarin (GB)', 'GD': 'soman (GD)',
       'GF': 'cyclosarin (GF)', 'VX': 'VX',
       'AC': 'hydrogen cyanide (AC)', 'CG': 'phosgene (CG)',
       'CK': 'cyanogen chloride (CK)', 'HD': 'sulfur mustard (HD)',
       'L': 'lewisite 1 (L)', 'PS': 'chloropicrin (PS)'}


# 2026-10-24. 중복 채널 축약. d2 의 자외선 피벗에서 Dye22~27 여섯 열이 35행
# 전체에 걸쳐 비트 단위로 같은 값을 갖는다. 원본 RGB 까지 거슬러 올라가 확인한
# 결과 파일명은 서로 다른데(B22AC10.jpg, B23AC10.jpg, ...) 읽힌 화소값이 같다.
# 저자 확인 결과 이 여섯 염료는 자외선에서 형광을 내지 않아 ROI 가 플레이트
# 배경을 읽었다. 즉 배경 1채널이 여섯 번 기록된 것이고, 단위 정규화한 패턴에서
# 그 한 채널이 6배 가중을 갖는다. 하필 12셀 중 계수가 가장 큰 두 셀이 여기서
# 나오므로, 접지 않으면 분석자 선택이 결과를 키우는 바로 그 구조가 된다.
# 규칙은 사전 선언한 것 하나뿐이다: 전 행에서 완전히 같은 열은 한 번만 센다.
CHANNEL_LOG = {}


def collapse_duplicate_channels(D, dy, key):
    keep, drop = [], {}
    for c in dy:
        v = D[c].to_numpy(float)
        same = next((k for k in keep if np.array_equal(D[k].to_numpy(float), v)),
                    None)
        if same is None:
            keep.append(c)
        else:
            drop.setdefault(same, []).append(c)
    CHANNEL_LOG[key] = dict(n_declared=len(dy), n_effective=len(keep),
                            collapsed={k: v for k, v in drop.items()},
                            all_zero=[c for c in dy
                                      if not D[c].to_numpy(float).any()])
    return keep


def load(tag, mode):
    D = pd.read_csv(f'deltaE/{tag}.{mode}_deltaE_pivot.csv')
    dy = [c for c in D.columns if c.lower().startswith('dye')]
    dy = collapse_duplicate_channels(D, dy, f'{tag}.{mode}')
    return D, dy


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def raw_patterns(D, dy):
    """물질 -> 농도 -> 단위 정규화한 응답. 차원은 중복 축약 후의 유효 채널 수."""
    out = {}
    for a in sorted(D.Agent.unique()):
        g = D[D.Agent == a]
        d = {}
        for c in CONC:
            r = g[g.Concentration == c]
            if len(r):
                v = r[dy].to_numpy(float).ravel()
                if np.linalg.norm(v) > 0:
                    d[c] = unit(v)
        if d:
            out[a] = d
    return out


def centred(raw):
    """각 농도에서 물질 축 평균을 빼 공통 응답 성분을 제거한다."""
    ags = sorted(raw)
    out = {a: {} for a in ags}
    for c in CONC:
        names = [a for a in ags if c in raw[a]]
        if len(names) < 2:
            continue
        V = np.array([raw[a][c] for a in names])
        mu = V.mean(0)
        for a, v in zip(names, V):
            out[a][c] = unit(v - mu)
    return {a: d for a, d in out.items() if d}


def consensus(cen):
    return {a: unit(np.mean([cen[a][c] for c in sorted(cen[a])], axis=0))
            for a in cen}


def discrimination():
    """변별력 게이트. 분석단위는 화합물이다 (2026-11-02, T2-3 에서 교체).

    왜 교체하는가. 종전 판은 물질내 코사인 쌍 전부(d1 이면 6종 x C(4,2) = 36)와
    물질간 합의패턴 쌍 전부(C(6,2) = 15)를 두 표본으로 보고 Mann-Whitney 를
    돌렸다. 그 36 과 15 는 독립 관측이 아니다. 같은 화합물에서 나온 쌍 여섯
    개가 한 화합물의 잡음을 공유하고, 물질간 쌍 15 개도 여섯 화합물을 돌려
    쓴다. 게다가 두 집단이 같은 중심화 벡터에서 만들어지므로 서로 독립도
    아니다. 그 결과 자유도가 실제보다 여섯 배 가까이 부풀어 p 가 10^-8 까지
    내려갔다. 이 논문은 다른 곳에서 바로 이 오류(쌍을 단위로 세는 것)를
    교정했으므로 여기에 남겨둘 수 없다.

    교체한 설계. 화합물 a 마다 두 수를 만든다.
      w_a = a 의 농도간 코사인 평균      (a 의 재현성)
      b_a = a 의 합의패턴과 다른 화합물 합의패턴의 코사인 평균 (a 의 변별)
    독립 단위는 화합물이므로 n 은 6, 7, 4 이고, d_a = w_a - b_a 에 대한
    대응 검정을 한다. n 이 작으므로 정규성 검정으로 t 와 비모수를 고르는 대신
    부호뒤집기 전수열거(2^n 가지)로 정확 p 를 낸다. n=6 이면 2^6 = 64 이므로
    단측 최소 p 는 1/64 = 0.0156 이고, n=4 면 1/16 = 0.0625 다. 즉 d3 의 두
    셀은 어떤 검정을 써도 0.05 밑으로 갈 수 없다. 그 사실을 숨기지 않는다.

    종전 p 는 p_pooled_pairs_superseded 로 남긴다. 지웠다고 적는 대신 왜
    쓰지 않는지 보이기 위한 것이다.
    """
    print('=' * 78)
    print('[게이트] 변별력 — 분석단위 = 화합물. 부호뒤집기 전수열거 정확검정')
    print('=' * 78)
    rows = []
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            D, dy = load(tag, mode)
            cen = centred(raw_patterns(D, dy))
            P = consensus(cen)
            names = sorted(cen)

            # 화합물 단위 요약
            w, b, per = [], [], []
            for a in names:
                ks = sorted(cen[a])
                wa = [float(np.dot(cen[a][i], cen[a][j]))
                      for i, j in itertools.combinations(ks, 2)]
                ba = [float(np.dot(P[a], P[o])) for o in names if o != a]
                if not wa or not ba:
                    continue
                w.append(float(np.mean(wa))); b.append(float(np.mean(ba)))
                # 2026-11-03 (T1). 평균 대비는 통과하지만 쌍별로는 어떤가.
                # 자기일치(농도간 평균 코사인)보다 더 가까운 타 화합물이
                # 있는지 센다. 있으면 그 셀에서 그 화합물은 '자기와 가장
                # 비슷하다'고 말할 수 없다. 초록과 본문이 평균 표현을 쓰는
                # 근거이자, 그 표현의 한계를 수치로 적기 위한 것이다.
                _worst = max(ba)
                _who = [o for o in names
                        if o != a and float(np.dot(P[a], P[o])) == _worst]
                per.append(dict(compound=a, n_concentrations=len(ks),
                                within=float(np.mean(wa)),
                                between=float(np.mean(ba)),
                                diff=float(np.mean(wa) - np.mean(ba)),
                                nearest_other=_who[0] if _who else None,
                                nearest_other_cos=float(_worst),
                                beats_nearest_other=bool(
                                    np.mean(wa) > _worst)))
            w = np.array(w); b = np.array(b); d = w - b
            n = len(d)

            # 부호뒤집기 전수열거. 귀무: d 의 부호가 교환가능.
            obs = float(d.mean())
            cnt = 0
            for signs in itertools.product([1, -1], repeat=n):
                if float((np.array(signs) * d).mean()) >= obs - 1e-15:
                    cnt += 1
            p_exact = cnt / (2 ** n)

            # 교차확인 두 가지와 효과크기
            try:
                p_wil = float(stats.wilcoxon(w, b, alternative='greater',
                                             mode='exact').pvalue)
            except Exception:
                p_wil = float('nan')
            sd = float(d.std(ddof=1)) if n > 1 else float('nan')
            dz = obs / sd if sd and sd == sd and sd > 0 else float('nan')
            # n 이 작아 t 분포 신뢰구간을 쓴다. 정확검정과 같은 가정을 쓰지
            # 않는다는 점을 SI 에 적는다.
            if n > 1 and sd > 0:
                tcrit = float(stats.t.ppf(0.975, n - 1))
                ci = [obs - tcrit * sd / np.sqrt(n), obs + tcrit * sd / np.sqrt(n)]
            else:
                ci = [float('nan'), float('nan')]

            # 종전(쌍 독립 가정) p. 쓰지 않지만 기록으로 남긴다.
            within_all, between_all = [], []
            for a in names:
                ks = sorted(cen[a])
                within_all += [float(np.dot(cen[a][i], cen[a][j]))
                               for i, j in itertools.combinations(ks, 2)]
            between_all = [float(np.dot(P[a], P[o]))
                           for a, o in itertools.combinations(names, 2)]
            p_old = float(stats.mannwhitneyu(within_all, between_all,
                                             alternative='greater').pvalue)

            _nviol = sum(1 for r in per if not r['beats_nearest_other'])
            print(f'  {tag} {mode:<7} n={n}  물질내 {w.mean():+.3f}  '
                  f'물질간 {b.mean():+.3f}  차 {obs:+.3f}  '
                  f'정확p={p_exact:.4f}  (최소가능 {1/2**n:.4f})  '
                  f'd_z={dz:+.2f}  쌍별위반 {_nviol}/{n}  '
                  f'[구판 쌍독립 p={p_old:.3g}]')
            rows.append(dict(dataset=tag, mode=mode,
                             unit='compound', n_compounds=n,
                             n_pairwise_violations=_nviol,
                             pairwise_violators=sorted(
                                 r['compound'] for r in per
                                 if not r['beats_nearest_other']),
                             within=float(w.mean()), between=float(b.mean()),
                             diff=obs, p=p_exact,
                             p_exact_signflip=p_exact,
                             p_min_attainable=1 / 2 ** n,
                             p_wilcoxon_exact=p_wil,
                             cohens_dz=dz, ci95_diff=ci,
                             n_within_pairs=len(within_all),
                             n_between_pairs=len(between_all),
                             p_pooled_pairs_superseded=p_old,
                             per_compound=per))
    return pd.DataFrame(rows)


def dmmp_question():
    """헤드라인 직접 검정: DMMP 의 실측 응답은 G 계열을 닮았는가 VX 를 닮았는가."""
    print()
    print('=' * 78)
    print('[헤드라인] DMMP 의 실측 29염료 패턴 (d1, 중심화 후)')
    print('=' * 78)
    rng = np.random.default_rng(SEED)
    rows = []
    for mode in ['uv', 'normal']:
        D, dy = load('d1', mode)
        cen = centred(raw_patterns(D, dy))
        P = consensus(cen)
        sims = {a: float(np.dot(P['DMMP'], P[a])) for a in P if a != 'DMMP'}
        g = float(np.mean([sims[a] for a in ['GB', 'GD', 'GF']]))
        diff = g - sims['VX']
        # 농도를 재표집한 부트스트랩 신뢰구간
        boots = []
        cs = sorted(set.intersection(*[set(cen[a]) for a in cen]))
        for _ in range(N_BOOT):
            pick = rng.choice(cs, size=len(cs), replace=True)
            Q = {a: unit(np.mean([cen[a][c] for c in pick], axis=0)) for a in cen}
            s = {a: float(np.dot(Q['DMMP'], Q[a])) for a in Q if a != 'DMMP'}
            boots.append(np.mean([s[a] for a in ['GB', 'GD', 'GF']]) - s['VX'])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        print(f'  [{mode}]  ' + '  '.join(f'{a} {v:+.3f}' for a, v in
                                          sorted(sims.items(), key=lambda t: -t[1])))
        print(f'          G계열(GB/GD/GF) 평균 {g:+.3f}   VX {sims["VX"]:+.3f}   '
              f'차 {diff:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  -> '
              f'{"G 계열" if diff > 0 else "VX"} 를 더 닮음')
        # 2026-11-03. g_series_mean 은 GB/GD/GF 3종 평균이다(위 273행).
        # 본문은 이것을 'the mean over the four G agents' 로 적고 있었고,
        # 숫자 자체는 산출물과 일치하므로 검증기 505건이 전부 통과했다.
        # 분모가 틀린 경우를 잡으려면 두 평균을 모두 예치해야 한다.
        g4 = float(np.mean([sims[a] for a in ['GA', 'GB', 'GD', 'GF']]))
        rows.append(dict(mode=mode, **{f'cos_{a}': v for a, v in sims.items()},
                         g_series_mean=g, g_series_n=3,
                         g4_series_mean=g4, g4_series_n=4,
                         vx=sims['VX'], g_minus_vx=diff,
                         ci_lo=float(lo), ci_hi=float(hi)))
    return pd.DataFrame(rows)


def space_test():
    """전제 검정: 어느 서술자 공간의 거리가 실측 패턴 비유사도를 예측하는가."""
    from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM,
                               OMEGA_DERIVED, PROTECTED, pick)
    print()
    print('=' * 78)
    print('[전제] 서술자 공간 거리 대 실측 ΔE 패턴 비유사도 (1 - 코사인)')
    print('=' * 78)
    df = pd.read_csv('descriptor_cache_xtb.csv')
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == 28 and c not in META_NUM]
    SPACES = {'structural': [c for c in num if c not in ELECTRONIC_28],
              'electronic': [c for c in num if c in ELECTRONIC_DIST],
              'merged': [c for c in num if c != OMEGA_DERIVED]}
    idx = {n: i for i, n in enumerate(df.compound_name)}
    DIST = {}
    for sp, cols in SPACES.items():
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Xz = (X - X.mean(0)) / X.std(0, ddof=1)
        DIST[sp] = np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1))

    rng = np.random.default_rng(SEED)
    rows = []
    for tag in ['d1', 'd2']:
        for mode in ['uv', 'normal']:
            D, dy = load(tag, mode)
            cen = centred(raw_patterns(D, dy))
            P = consensus(cen)
            ags = [a for a in sorted(P) if a in MAP]
            pairs = list(itertools.combinations(ags, 2))
            dissim = np.array([1 - float(np.dot(P[a], P[b])) for a, b in pairs])
            print(f'\n  --- {tag} {mode}  ({len(ags)}종, {len(pairs)}쌍: '
                  f'{", ".join(ags)})')
            for sp in SPACES:
                d = np.array([DIST[sp][idx[MAP[a]], idx[MAP[b]]] for a, b in pairs])
                r, p = stats.spearmanr(d, dissim)
                bs = [stats.spearmanr(d[k], dissim[k]).statistic
                      for k in (rng.integers(0, len(d), len(d))
                                for _ in range(N_BOOT))]
                bs = [b for b in bs if b == b]
                lo, hi = np.percentile(bs, [2.5, 97.5])
                print(f'      {sp:<11} Spearman {r:+.3f}  p={p:.3f}  '
                      f'95% CI [{lo:+.3f}, {hi:+.3f}]')
                rows.append(dict(dataset=tag, mode=mode, space=sp, n_pairs=len(pairs),
                                 spearman=float(r), p=float(p),
                                 ci_lo=float(lo), ci_hi=float(hi)))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    dis = discrimination()
    dm = dmmp_question()
    sp = space_test()
    json.dump(dict(seed=SEED, bootstrap=N_BOOT,
                   discrimination=dis.to_dict('records'),
                   dmmp=dm.to_dict('records'),
                   space_test=sp.to_dict('records')),
              open('outputs/deltaE_premise.json', 'w'), indent=1, ensure_ascii=False)
    # per_compound 는 중첩 레코드이므로 CSV 에서 뺀다. JSON 에 그대로 있다.
    dis.drop(columns=['per_compound']).to_csv(
        'outputs/deltaE_discrimination.csv', index=False)
    sp.to_csv('outputs/deltaE_space_test.csv', index=False)
    print('\n[done] outputs/deltaE_premise.json, deltaE_discrimination.csv, '
          'deltaE_space_test.csv')
