#!/usr/bin/env python3
"""Step 3 롤백 (방향 C) 분석 — 2026-10-13.

GATE 8 3차에서 4인 중 3인이 치명으로 지목한 두 지점을 분석 수준에서 해소하고,
주장을 메타 결과로 내린다. 수행 항목은 심사 지적 번호와 1:1 대응한다.

  C1/M3 전자공간 재정의.  omega = mu^2/2eta 는 HOMO/LUMO 의 결정론적 함수이므로
                          거리공간에서 제외한다. 전 패널에 결측 없는 전자 좌표는
                          E_HOMO, E_LUMO, dipole 세 개뿐이다 (P 특이 서술자는
                          비인 화합물 22/28 결측).
  C2      차원 비대칭 통제. 구조공간 51차원 대 전자공간 3차원을 그대로 비교하면
                          결과가 화학인지 기하인지 분리되지 않는다. (a) 구조공간을
                          주성분 3개로 절단, (b) 구조 서술자 무작위 3종 조합의
                          p 분포를 귀무 기준선으로 제시한다.
  S1/M1/P2 패널 두 구성(28, 35)을 동등하게 보고한다. 어느 쪽도 주 분석이 아니다.
  C3/S6   ΔE 공간 검정을 Mantel 순열로 교체한다. 쌍 15개는 독립 관측이 아니다.
  C4/M5b  중심화 후 물질간 코사인의 영가설을 0 이 아니라 -1/(n-1) 로 둔다.
  P5/M5a  DMMP 의 유사 작용제 전체 순위를 보고한다. GA 가 d1 uv 1순위다.

출력: outputs/rollback_c.json
"""
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats

from spaces import ELECTRONIC_28, META_NUM, PROTECTED, pick
from perm_statistics import STATS

SEED = 0xC0DE
N_MANTEL = 10000
N_RANDSUB = 500

# C1. omega 는 거리공간에서 제외한다. 표에는 유도량으로 남긴다.
OMEGA = 'electrophilicity_index_eV'
ELEC3 = [c for c in ELECTRONIC_28 if c != OMEGA]


# ----------------------------------------------------------------- 공통 유틸
def zdist(df, cols):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    return D, len(keep)


def exact_perm(D, kl, sim_classes):
    """원고와 같은 귀무. 한 라벨링에서 선택된 k 종이 모사체, 풀의 나머지가 간섭물."""
    D = D.copy()
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    null = {s: np.empty(len(ALL)) for s in STATS}
    npool = len(i_pool)
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
        null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
        null['mean_logratio'][a0:a0 + b] = np.log(ds / di).mean(0)
        null['mean_dsim'][a0:a0 + b] = ds.mean(0)
        null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
        null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    ti = ALL.index(true_combo)
    return {s: dict(observed=float(null[s][ti]),
                    p=float((null[s] <= null[s][ti]).sum() / len(ALL)))
            for s in STATS}, len(ALL)


# ------------------------------------------------- C1 omega 함수종속 정량화
def omega_dependence():
    print('=' * 78)
    print('[C1/M3] omega 의 함수 종속 — 전자공간의 실제 자유도')
    print('=' * 78)
    out = {}
    for tag, path in [('panel28', 'descriptor_cache_xtb.csv'),
                      ('panel35', 'descriptor_cache_xtb_35.csv')]:
        d = pd.read_csv(path)
        H, L, w = (d.E_HOMO_eV.values, d.E_LUMO_eV.values, d[OMEGA].values)
        mu, eta = (H + L) / 2, (L - H) / 2
        dev = float(np.abs(w - mu * mu / (2 * eta)).max())
        E4 = d[ELECTRONIC_28]
        Z4 = ((E4 - E4.mean()) / E4.std(ddof=1)).values
        sv = np.linalg.svd(Z4, compute_uv=False)
        cum3 = float((sv[:3] ** 2).sum() / (sv ** 2).sum())
        r_lumo = float(np.corrcoef(L, w)[0, 1])
        D4, _ = zdist(d, ELECTRONIC_28)
        D3, _ = zdist(d, ELEC3)
        iu = np.triu_indices(len(d), 1)
        rho = float(stats.spearmanr(D4[iu], D3[iu]).statistic)
        print(f'  {tag}: |omega - mu^2/2eta| 최대 {dev:.2e}   '
              f'r(E_LUMO, omega) = {r_lumo:+.3f}   '
              f'주성분3 누적분산 {100*cum3:.1f}%   '
              f'4좌표 대 3좌표 거리 Spearman {rho:.4f}')
        out[tag] = dict(max_deviation_from_definition=dev,
                        r_lumo_omega=r_lumo, pc3_cumulative_variance=cum3,
                        dist_rank_correlation_4d_vs_3d=rho,
                        singular_values=[float(s) for s in sv])
    print('  판정: omega 는 독립 좌표가 아니다. 거리공간에서 제외하고 표에만 남긴다.')
    return out


# --------------------------------------- 패널 2종 x 공간 3종 (omega 제외판)
def panels():
    print()
    print('=' * 78)
    print('[S1/M1/P2] 패널 두 구성을 동등하게. 전자공간은 omega 제외 3좌표')
    print('=' * 78)
    CASES = [
        ('panel28_designated3', 'descriptor_cache_xtb.csv', ['simulant']),
        ('panel35_designated3', 'descriptor_cache_xtb_35.csv', ['simulant']),
        ('panel35_designated_plus_literature8', 'descriptor_cache_xtb_35.csv',
         ['simulant', 'literature_simulant']),
    ]
    out = {}
    for label, path, sim_cls in CASES:
        df = pd.read_csv(path)
        n = len(df)
        num = [c for c in df.select_dtypes('number').columns
               if df[c].notna().sum() == n and c not in META_NUM]
        kl = df.compound_class.tolist()
        if not any(k in sim_cls for k in kl):
            print(f'  [건너뜀] {label}: 해당 class 없음 '
                  f'(존재: {sorted(set(kl))})')
            continue
        SP = {'structural': [c for c in num if c not in ELECTRONIC_28],
              'electronic3': [c for c in num if c in ELEC3],
              'merged': [c for c in num if c != OMEGA]}
        print(f'\n  ### {label}  ({n}종, 모사체 {sum(k in sim_cls for k in kl)})')
        out[label] = {}
        for sp, cols in SP.items():
            D, nk = zdist(df, cols)
            row, nperm = exact_perm(D, kl, sim_cls)
            out[label][sp] = dict(n_descriptors=nk, n_perm=nperm, statistics=row)
            ps = ', '.join(f'{s}={row[s]["p"]:.3f}' for s in STATS)
            print(f'    {sp:<12} ({nk:>2}차원, 귀무 {nperm:,})  {ps}')
    return out


# ------------------------------------------------------ C2 차원 비대칭 통제
def dimensionality_control(elec_p):
    global ELEC_P
    ELEC_P = elec_p
    print()
    print('=' * 78)
    print('[C2] 차원 비대칭 통제 — 구조공간을 3차원으로 내리면 어떻게 되는가')
    print('=' * 78)
    rng = np.random.default_rng(SEED)
    out = {}
    for tag, path, sim_cls in [
            ('panel28', 'descriptor_cache_xtb.csv', ['simulant']),
            ('panel35', 'descriptor_cache_xtb_35.csv', ['simulant'])]:
        df = pd.read_csv(path)
        n = len(df)
        num = [c for c in df.select_dtypes('number').columns
               if df[c].notna().sum() == n and c not in META_NUM]
        scols = [c for c in num if c not in ELECTRONIC_28]
        kl = df.compound_class.tolist()

        keep, _ = pick(df, scols, PROTECTED)
        X = df[keep].to_numpy(float)
        Z = (X - X.mean(0)) / X.std(0, ddof=1)
        full, _ = exact_perm(np.sqrt(((Z[:, None] - Z[None]) ** 2).sum(-1)),
                             kl, sim_cls)

        # (a) 주성분 3개로 절단
        U, S, _ = np.linalg.svd(Z, full_matrices=False)
        P3 = U[:, :3] * S[:3]
        pc3, _ = exact_perm(np.sqrt(((P3[:, None] - P3[None]) ** 2).sum(-1)),
                            kl, sim_cls)

        # (b) 무작위 3좌표 부분공간 통제는 stat_supplement.py 가 2000회로,
        # 백색화 전후 둘 다 수행한다 (2026-10-14, 리뷰어 2 R2-2/R2-3).
        # 같은 양을 두 스크립트가 다른 난수로 산출하면 본문이 어느 쪽을 쓰는지
        # 모호해지므로 여기서는 계산하지 않는다.
        print(f'\n  ### {tag}  (구조 {len(keep)}차원)')
        print(f'    구조 전차원        mean_diff p = {full["mean_diff"]["p"]:.4f}')
        print(f'    구조 주성분 3개    mean_diff p = {pc3["mean_diff"]["p"]:.4f}')
        print('    구조 무작위 3좌표   stat_supplement.py 참조 (여기서 미수행)')
        # 2026-10-24. 종전에는 미수행 경로가 random3_p_quantiles 를 전부 NaN,
        # random3_frac_below_05 를 0.000 으로, n_random_subsets 를 500 으로
        # 적어 기탁 JSON 에 내보냈다. 돌리지 않은 500회를 돌린 것처럼 읽히고,
        # 본문이 쓰는 값(2000회, 36.3%)과도 다르다. 키 자체를 내보내지 않고
        # 어디서 산출되는지만 남긴다.
        out[tag] = dict(n_structural_dims=len(keep),
                        full_p=full['mean_diff']['p'],
                        pc3_p=pc3['mean_diff']['p'],
                        random3_subspace_control='stat_supplement.py')
    return out


# ------------------------------------------------- C3 Mantel, C4 중심화 바닥
def deltaE_corrected():
    import deltaE_premise as m
    print()
    print('=' * 78)
    print('[C3/S6] ΔE 공간 검정을 Mantel 순열로. [C4/M5b] 중심화 바닥 영가설')
    print('=' * 78)
    df = pd.read_csv('descriptor_cache_xtb.csv')
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == 28 and c not in META_NUM]
    # 2026-10-24. 종전에는 구조·전자3 두 공간만 돌았다. 본문의 1차 공간인
    # 병합이 빠져 있었고, 빠진 셀 중 하나(d2 uv merged)가 계수 최대다.
    # 사후에 family 를 좁힌 꼴이고 좁힌 방향이 결과를 완화하는 쪽이었다.
    # 패널·공간·통계량을 사전 고정하라고 결론내는 논문에서 자기모순이다.
    _struct = [c for c in num if c not in ELECTRONIC_28]
    _elec3 = [c for c in num if c in ELEC3]
    SPACES = {'structural': _struct, 'electronic3': _elec3,
              'merged': _struct + _elec3}
    idx = {nm: i for i, nm in enumerate(df.compound_name)}
    DIST = {sp: zdist(df, cols)[0] for sp, cols in SPACES.items()}
    rng = np.random.default_rng(SEED)
    mantel, floor = [], []

    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                D, dy = m.load(tag, mode)
            except Exception as e:
                print(f'  [수집 불가: {tag} {mode} {e}]')
                continue
            P = m.consensus(m.centred(m.raw_patterns(D, dy)))
            ags = [a for a in sorted(P) if a in m.MAP]
            if len(ags) < 4:
                continue
            k = len(ags)
            pr = list(itertools.combinations(range(k), 2))
            dis = np.array([1 - float(np.dot(P[ags[i]], P[ags[j]]))
                            for i, j in pr])

            # C4. 중심화는 물질간 코사인 평균을 -1/(n-1) 에 고정한다.
            allc = [float(np.dot(P[a], P[b]))
                    for a, b in itertools.combinations(sorted(P), 2)]
            nall = len(P)
            floor.append(dict(dataset=tag, mode=mode, n_compounds=nall,
                              observed_mean_cosine=float(np.mean(allc)),
                              algebraic_floor=float(-1 / (nall - 1)),
                              gap=float(np.mean(allc) + 1 / (nall - 1))))

            for sp in SPACES:
                M = np.array([[DIST[sp][idx[m.MAP[a]], idx[m.MAP[b]]]
                               for b in ags] for a in ags])
                d = np.array([M[i, j] for i, j in pr])
                r = float(stats.spearmanr(d, dis).statistic)
                # Mantel: 화합물 라벨을 순열한다. 쌍을 순열하지 않는다.
                cnt = 0
                for _ in range(N_MANTEL):
                    pm = rng.permutation(k)
                    dp = np.array([M[pm[i], pm[j]] for i, j in pr])
                    rp = stats.spearmanr(dp, dis).statistic
                    if rp == rp and abs(rp) >= abs(r):
                        cnt += 1
                pm_ = (cnt + 1) / (N_MANTEL + 1)
                # 종전(쌍 독립 가정) p 를 나란히 둔다
                p_naive = float(stats.spearmanr(d, dis).pvalue)
                mantel.append(dict(dataset=tag, mode=mode, space=sp,
                                   n_compounds=k, n_pairs=len(pr),
                                   spearman=r, p_mantel=float(pm_),
                                   p_pairs_as_independent=p_naive))
                print(f'  {tag} {mode:<7} {sp:<12} rho {r:+.3f}   '
                      f'Mantel p = {pm_:.4f}   '
                      f'(쌍독립가정 p = {p_naive:.4f})')
    # 12셀은 하나의 family 다. 순열검정 18건을 보정하면서 이쪽을 보정하지 않으면
    # 논문이 스스로 비판하는 오류를 범하게 된다. 2026-10-14 추가,
    # 2026-10-24 에 병합공간 4셀을 되돌려 8셀 -> 12셀.
    pv = np.array([x['p_mantel'] for x in mantel])
    k = len(pv)
    o = np.sort(pv)
    bh = np.minimum.accumulate((o * k / np.arange(1, k + 1))[::-1])[::-1]
    mantel_bh = dict(n_tests=int(k), min_q=float(bh[0]),
                     n_sig_q=int((bh < 0.05).sum()))
    print(f'  [BH] Mantel {k}셀 최소 q = {bh[0]:.3f}, q<0.05 {int((bh<0.05).sum())}건')
    print()
    print('  중심화 바닥 대조 (관측 평균 코사인 대 -1/(n-1)):')
    for f in floor:
        print(f'    {f["dataset"]} {f["mode"]:<7} n={f["n_compounds"]}  '
              f'관측 {f["observed_mean_cosine"]:+.3f}  '
              f'바닥 {f["algebraic_floor"]:+.3f}  '
              f'여유 {f["gap"]:+.3f}')
    # 2026-10-24. 중복 채널 축약 기록을 기탁 산출물에 싣는다. SI 가 유효
    # 채널 수를 인용하므로 그 근거가 저장소 안에 있어야 한다.
    return dict(mantel=mantel, centring_floor=floor, mantel_bh=mantel_bh,
                channels=m.CHANNEL_LOG)


# ------------------------------------------------------------ P5 GA 전체 순위
def dmmp_full_ranking():
    import deltaE_premise as m
    print()
    print('=' * 78)
    print('[P5/M5a] DMMP 의 유사 작용제 전체 순위 — GA 를 포함해 전량 보고')
    print('=' * 78)
    out = []
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                D, dy = m.load(tag, mode)
            except Exception:
                continue
            P = m.consensus(m.centred(m.raw_patterns(D, dy)))
            if 'DMMP' not in P:
                continue
            r = sorted(((float(np.dot(P['DMMP'], P[k])), k)
                        for k in P if k != 'DMMP'), reverse=True)
            out.append(dict(dataset=tag, mode=mode,
                            ranking=[dict(agent=k, cosine=v) for v, k in r]))
            print(f'  {tag} {mode:<7} ' +
                  ', '.join(f'{k} {v:+.3f}' for v, k in r))
    return out


if __name__ == '__main__':
    res = dict(seed=SEED, n_mantel=N_MANTEL, n_random_subsets=N_RANDSUB,
               omega_excluded_from_distance=OMEGA,
               electronic_coordinates=ELEC3)
    res['omega_dependence'] = omega_dependence()
    res['panels'] = panels()
    res['dimensionality_control'] = dimensionality_control(
        res['panels']['panel28_designated3']['electronic3']['statistics']['mean_diff']['p'])
    res['deltaE'] = deltaE_corrected()
    res['dmmp_ranking'] = dmmp_full_ranking()
    with open('outputs/rollback_c.json', 'w') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print('\n[rollback_c] outputs/rollback_c.json 기록')
