#!/usr/bin/env python3
"""통계 보완 (2026-10-14). GATE 8 4차 리뷰어 2 지적 R2-1 / R2-2 / R2-3 대응.

R2-1  18검정을 하나의 family 로 묶은 것이 자의적이다. 세 공간은 중첩이고
      (병합 = 구조 U 전자) 통계량 6종은 같은 거리행렬의 함수이므로 검정들이
      독립이 아니다. BH 는 양의 의존에서 보수적이므로 q=0.159 가 과대일 수 있고,
      그렇다면 '아무 것도 살아남지 않는다' 는 결론이 보정 방식의 선택에 달린다.
      논문의 논지가 '선택이 답을 정한다' 인데 그 논지를 받치는 보정이 선택이면
      자기 적용 실패다. 그래서
        (a) family 를 공간별 6검정으로 좁혔을 때
        (b) 통계량 1종(mean_diff)을 사전지정했을 때
        (c) 귀무분포들의 실제 상관에서 실질 검정 수 meff 를 추정했을 때
      의 보정값을 모두 낸다. meff 는 Li & Ji (2005) 고유값 방법.

R2-2  구조 서술자 38개는 상호 상관이 높다(PC1 36.5%). 상관된 풀에서 3개를
      뽑으면 실질 자유도가 3보다 작으므로 p<0.05 도달률 39.8% 가 부풀려졌을 수
      있다. ZCA 백색화로 좌표를 무상관화한 뒤 같은 통제를 돌려 대조한다.

R2-3  500회는 백분위를 0.2 까지만 분해한다. 2000회로 올린다.

출력: outputs/stat_supplement.json
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, ELECTRONIC_DIST, META_NUM, OMEGA_DERIVED, PROTECTED, pick

SEED = 0xC0DE
N_RANDSUB = 2000
STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']


def null_matrix(D, kl, sim_classes):
    """통계량 6종의 귀무분포를 (n_perm, 6) 으로 돌려준다. 관측치 인덱스도."""
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    pos = {g: j for j, g in enumerate(i_pool)}
    true = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    npool = len(i_pool)
    out = np.empty((len(ALL), 6))
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None], A3, np.inf).min(2)
        di = np.where(~m[None], A3, np.inf).min(2)
        out[a0:a0 + b, 0] = (ds - di).mean(0)
        out[a0:a0 + b, 1] = np.median(ds - di, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            out[a0:a0 + b, 2] = np.log(ds / di).mean(0)
        out[a0:a0 + b, 3] = ds.mean(0)
        out[a0:a0 + b, 4] = -(ds < di).sum(0)
        out[a0:a0 + b, 5] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    return out, ALL.index(true)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate(
        (p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def meff_li_ji(R):
    """Li & Ji (2005) 고유값 기반 실질 검정 수."""
    lam = np.linalg.eigvalsh(R)
    lam = np.clip(lam, 0, None)
    return float(sum((l >= 1) + (l - np.floor(l)) for l in lam))


def spaces_for(df, num):
    return {'structural': [c for c in num if c not in ELECTRONIC_28],
            'electronic': [c for c in num if c in ELECTRONIC_DIST],
            'merged': [c for c in num if c != OMEGA_DERIVED]}


def zdist(Z):
    return np.sqrt(((Z[:, None] - Z[None]) ** 2).sum(-1))


def family_sensitivity(path, sim_classes, label):
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    kl = df.compound_class.tolist()
    nulls, obs_p, names = [], [], []
    for sp, cols in spaces_for(df, num).items():
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        N, ti = null_matrix(zdist(Z), kl, sim_classes)
        for j, st in enumerate(STATS):
            col = N[:, j]
            good = np.isfinite(col)
            obs_p.append(float((col[good] <= col[ti]).sum() / good.sum()))
            nulls.append(np.where(good, col, np.nan))
            names.append(f'{sp}:{st}')
    P = np.array(obs_p)

    # (a) 전체 18 / (b) 공간별 6 / (c) 통계량 사전지정 3
    q_all = bh(P)
    q_by_space = np.concatenate([bh(P[i:i + 6]) for i in (0, 6, 12)])
    idx_md = [names.index(f'{s}:mean_diff')
              for s in ('structural', 'electronic', 'merged')]
    q_pre = bh(P[idx_md])

    # (c) 의존을 고려한 실질 검정 수
    M = np.vstack(nulls)
    R = pd.DataFrame(M.T).corr().to_numpy()
    R = np.nan_to_num(R, nan=0.0)
    np.fill_diagonal(R, 1.0)
    me = meff_li_ji(R)

    res = dict(
        label=label, n_tests=len(P), min_p=float(P.min()),
        bh_all_min_q=float(q_all.min()),
        bh_all_n_sig=int((q_all < 0.05).sum()),
        bh_by_space_min_q=float(q_by_space.min()),
        bh_by_space_n_sig=int((q_by_space < 0.05).sum()),
        prespecified_min_p=float(P[idx_md].min()),
        prespecified_bh_min_q=float(q_pre.min()),
        prespecified_n_sig=int((q_pre < 0.05).sum()),
        meff=me,
        meff_adjusted_min_p=float(min(1.0, P.min() * me)),
        bonferroni_min_p=float(min(1.0, P.min() * len(P))))
    print(f'\n### {label}')
    print(f'  최소 명목 p            {res["min_p"]:.4f}')
    print(f'  (a) 18검정 BH          최소 q {res["bh_all_min_q"]:.3f}  '
          f'유의 {res["bh_all_n_sig"]}건')
    print(f'  (b) 공간별 6검정 BH    최소 q {res["bh_by_space_min_q"]:.3f}  '
          f'유의 {res["bh_by_space_n_sig"]}건')
    print(f'  (c) 통계량 사전지정    최소 p {res["prespecified_min_p"]:.4f}  '
          f'BH 최소 q {res["prespecified_bh_min_q"]:.3f}  '
          f'유의 {res["prespecified_n_sig"]}건')
    print(f'  실질 검정 수 meff      {me:.2f} / 18   '
          f'-> meff 보정 최소 p {res["meff_adjusted_min_p"]:.3f}  '
          f'(Bonferroni {res["bonferroni_min_p"]:.3f})')
    return res


def whitened_control(path, sim_classes, label):
    """R2-2. 백색화 전후로 무작위 3좌표 부분공간의 p 분포를 비교한다."""
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    kl = df.compound_class.tolist()
    keep, _ = pick(df, [c for c in num if c not in ELECTRONIC_28], PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)

    # ZCA 백색화. 공분산이 단위행렬이 되어 좌표가 무상관이 된다.
    C = np.cov(Z, rowvar=False)
    lam, V = np.linalg.eigh(C)
    keep_ax = lam > 1e-8
    W = V[:, keep_ax] @ np.diag(1 / np.sqrt(lam[keep_ax])) @ V[:, keep_ax].T
    Zw = Z @ W
    Zw = (Zw - Zw.mean(0)) / Zw.std(0, ddof=1)

    out = {}
    for tag, Zx in [('raw', Z), ('whitened', Zw)]:
        rng = np.random.default_rng(SEED)
        ps = []
        for _ in range(N_RANDSUB):
            c = rng.choice(Zx.shape[1], 3, replace=False)
            N, ti = null_matrix(zdist(Zx[:, c]), kl, sim_classes)
            col = N[:, 0]
            ps.append(float((col <= col[ti]).sum() / len(col)))
        ps = np.array(ps)
        qs = np.percentile(ps, [5, 25, 50, 75, 95])
        out[tag] = dict(n_draws=N_RANDSUB, frac_below_05=float((ps < 0.05).mean()),
                        quantiles=dict(zip(['q05', 'q25', 'q50', 'q75', 'q95'],
                                           [float(v) for v in qs])),
                        mc_se_pct=float(100 * np.sqrt(
                            (ps < 0.05).mean() * (1 - (ps < 0.05).mean())
                            / N_RANDSUB)))
        print(f'  {label} {tag:<9} p<0.05 {100*out[tag]["frac_below_05"]:.1f}% '
              f'(MC SE {out[tag]["mc_se_pct"]:.2f}%p)  중위 {qs[2]:.3f}  '
              f'5-95% [{qs[0]:.3f}, {qs[4]:.3f}]')

    # 관측값 백분위를 2000회 분포에 다시 매긴다
    rng = np.random.default_rng(SEED)
    ps = []
    for _ in range(N_RANDSUB):
        c = rng.choice(Z.shape[1], 3, replace=False)
        N, ti = null_matrix(zdist(Z[:, c]), kl, sim_classes)
        col = N[:, 0]
        ps.append(float((col <= col[ti]).sum() / len(col)))
    ps = np.array(ps)
    return out, ps


def threshold_free_match(path, sim_classes, label):
    """유효매칭 2조건 중 '지정집합 밖보다 가깝다' 는 임계값을 쓰지 않는다.
    그 조건만 충족하는 작용제 수를 센다. 임계값을 무한대로 둔 매칭 수와 같다."""
    df = pd.read_csv(path)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    keep, _ = pick(df, [c for c in num if c != OMEGA_DERIVED], PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = zdist(Z)
    np.fill_diagonal(D, np.inf)
    kl = df.compound_class.tolist()
    ia = [i for i, k in enumerate(kl) if k == 'real_agent']
    isim = [i for i, k in enumerate(kl) if k in sim_classes]
    iout = [i for i, k in enumerate(kl)
            if k not in sim_classes and k != 'real_agent']
    k = int(sum(D[a, isim].min() <= D[a, iout].min() for a in ia))
    print(f'  {label:<26} 임계값 무관 매칭 {k}/11  '
          f'(귀무기대 "밖이 더 가까움" {len(ia)*len(iout)/(len(isim)+len(iout)):.2f}/11)')
    return dict(n_match_unbounded=k, n_agents=len(ia),
                null_expectation_outside_closer=float(
                    len(ia) * len(iout) / (len(isim) + len(iout))))


if __name__ == '__main__':
    res = {'seed': SEED, 'n_random_subsets': N_RANDSUB}
    print('=' * 78)
    print('[R2-1] family 정의 민감도와 의존을 고려한 실질 검정 수')
    print('=' * 78)
    res['family'] = {}
    for lab, path, sc in [
            ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
            ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
            ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
             ['simulant', 'literature_simulant'])]:
        res['family'][lab] = family_sensitivity(path, sc, lab)

    print()
    print('=' * 78)
    print('[임계값 무관성] 유효매칭의 두 번째 조건은 임계값을 쓰지 않는다')
    print('=' * 78)
    res['threshold_free'] = {
        'panel35_incumbent3': threshold_free_match(
            'descriptor_cache_xtb_35.csv', ['simulant'], 'panel35 incumbent3'),
        'panel28_incumbent3': threshold_free_match(
            'descriptor_cache_xtb.csv', ['simulant'], 'panel28 incumbent3')}

    print()
    print('=' * 78)
    print(f'[R2-2/R2-3] 백색화 대조, 무작위 3좌표 {N_RANDSUB}회')
    print('=' * 78)
    w28, ps28 = whitened_control('descriptor_cache_xtb.csv', ['simulant'],
                                 'panel28')
    res['whitening'] = {'panel28': w28}
    for lab, v in [('struct_full', 0.0279), ('struct_pc3', 0.0044),
                   ('electronic', 0.2750)]:
        res.setdefault('percentiles_2000', {})[lab] = float(
            100 * (ps28 <= v).mean())
    print('\n  2000회 분포에서의 관측 백분위:')
    for k, v in res['percentiles_2000'].items():
        print(f'    {k:<14} {v:.1f}%')

    with open('outputs/dim_percentiles.json', 'w') as fh:
        json.dump(dict(struct_full=res['percentiles_2000']['struct_full'],
                       struct_pc3=res['percentiles_2000']['struct_pc3'],
                       electronic=res['percentiles_2000']['electronic'],
                       n_subsets=N_RANDSUB), fh, indent=1)

    with open('outputs/stat_supplement.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('\n[stat_supplement] outputs/stat_supplement.json 기록')
