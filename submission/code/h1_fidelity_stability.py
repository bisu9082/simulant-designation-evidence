#!/usr/bin/env python3
"""충실도 축이 전자구조 방법오차를 통과하는가 (2026-11-18 신설).

왜 필요한가. h1 의 fidelity 는 GFN2-xTB 세 좌표에서 최근접 진본 작용제까지의
표준화 거리다. 이 논문은 같은 세 좌표를 B3LYP 로 다시 만들면 순열 결과가
p=0.275 에서 0.0368 로 움직인다고 보고한다. 그런데 서술자 거리를 유일하게
건설적으로 쓰는 이 축에는 그 방법오차가 통과되지 않았다.

모형. 좌표별 표준화가 아핀 변환(shift·scale)을 정확히 흡수하므로, 방법
차이의 계통 성분은 이 축에 작용하지 않는다. 먼저 그것을 결정론적으로
확인하고(보정 전후 축이 같은지), 남는 비아핀 잔차만 전파한다.

  1. 보정집합 28종. GFN2 는 outputs/xtb_fast_vs_slow.csv 의 fast 열
     (h1 축이 쓰는 값), B3LYP 는 outputs/dft_dipole.csv.
  2. 좌표별 OLS  B3LYP = a + b * GFN2. 잔차 SD 를 기록한다.
  3. 보정축을 237 종 전체에 적용한다. 잔차는 B3LYP 척도에 있으므로
     반드시 보정축에 더한다. 원 GFN2 축에 더하면 LUMO 기울기 때문에
     잔차의 상대 크기가 과소평가된다.
  4. 28종의 잔차 3-벡터를 화합물 단위로 복원추출(좌표 간 상관 보존)해
     전체에 더하고, 표준화·최근접작용제거리를 다시 계산한다. B=4000.

[C-3] 새 양자화학 계산은 하지 않는다. 이미 예치된 28종 B3LYP 와 GFN2 를
쓰고, 그 둘의 차이만 전파한다. 207 종에 B3LYP 가 없다는 사실이 이 설계의
전제이고, 결과는 그 전제 안에서만 읽어야 한다.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

EL = ['E_HOMO_eV', 'E_LUMO_eV', 'dipole_moment_D']
DFT = ['homo_eV', 'lumo_eV', 'dipole_D']
SEED, B = 0xC0DE, 4000


def nearest_agent(Z, is_agent):
    mu, sd = Z.mean(0), Z.std(0, ddof=1)
    Zz = (Z - mu) / sd
    A = Zz[is_agent]
    return np.array([np.sqrt(((A - z) ** 2).sum(1)).min() for z in Zz])


def main():
    pool = pd.read_csv('outputs/h1_pool.csv')
    isag = (pool.compound_class == 'real_agent').to_numpy()
    cand = pd.read_csv('outputs/h1_candidates.csv')
    _fid = dict(zip(cand.label, cand.fidelity_dist))
    has = np.array([l in _fid for l in pool.label])
    X = pool[EL].to_numpy(float)

    fid0 = nearest_agent(X, isag)
    dev = np.max(np.abs(fid0[has]
                        - np.array([_fid[l] for l in pool.label[has]])))
    print(f'[check] 예치 fidelity 재현 최대편차 {dev:.3e}')
    assert dev < 1e-9

    fast = pd.read_csv('outputs/xtb_fast_vs_slow.csv')
    dft = pd.read_csv('outputs/dft_dipole.csv')
    m = fast.merge(dft, left_on='compound_name', right_on='name',
                   suffixes=('_g', '_d'))
    print(f'[calib] 보정집합 {len(m)}종')
    a, b, res = [], [], []
    for cg, cd in zip(EL, DFT):
        g, d = m[cg].to_numpy(float), m[cd].to_numpy(float)
        sl, ic = np.polyfit(g, d, 1)
        r = d - (ic + sl * g)
        a.append(ic)
        b.append(sl)
        res.append(r)
        print(f'   {cg:18s} a {ic:+.4f}  b {sl:+.4f}  '
              f'R2 {1 - r.var() / d.var():.3f}  잔차SD {r.std(ddof=2):.4f}')
    a, b = np.array(a), np.array(b)
    R = np.column_stack(res)

    Xc = a + b * X
    fid_c = nearest_agent(Xc, isag)
    aff = float(np.max(np.abs(fid_c - fid0)))
    rho0 = stats.spearmanr(fid0[has], pool.y_tox.to_numpy(float)[has])[0]
    rho_c = stats.spearmanr(fid_c[has], pool.y_tox.to_numpy(float)[has])[0]
    print(f'[affine] 보정 전후 축 최대차 {aff:.3e}, Spearman {rho0:+.6f} -> '
          f'{rho_c:+.6f}  (표준화가 아핀을 흡수한다)')

    rng = np.random.default_rng(SEED)
    y = pool.y_tox.to_numpy(float)
    cut0 = float(np.percentile(fid_c[has], 25))
    labs = pool.label[has].to_numpy()
    sl_names = json.load(open('outputs/h1_summary.json'))[
        'safety_tiers']['relaxed_10']['labels']
    dcm_i = int(np.where(labs == 'dichloromethane')[0][0])
    # 2026-11-18 (GATE 8 3차). 현행 경계 사분면에도 같은 방법오차를
    # 통과시킨다. 종전에는 완화 경계의 개수(여덟)에만 잔차를 전파하고,
    # 현행 경계의 개수(0)는 'outright' 라는 부사까지 붙여 확정값으로
    # 적었다. 두 값이 같은 축·같은 사분위 경계에서 나오므로 비대칭이다.
    _H = json.load(open('outputs/h1_summary.json'))
    bar_inc = float(_H['tox_cut'])
    alert = (pd.read_csv('outputs/h1_candidates.csv')
             .set_index('label').alert_flag.to_dict())
    noalert = np.array([not bool(alert.get(l, False)) for l in labs])
    safe_inc = (y[has] >= bar_inc) & noalert

    draws = np.empty((B, int(has.sum())))
    rhos = np.empty(B)
    top1 = []
    for k in range(B):
        idx = rng.integers(0, len(R), size=len(X))
        Xb = Xc + R[idx]
        fb = nearest_agent(Xb, isag)[has]
        draws[k] = fb
        rhos[k] = stats.spearmanr(fb, y[has])[0]
        top1.append(labs[int(np.argmin(fb))])

    within = draws.var(0).mean()
    between = fid_c[has].var()
    lam = 1.0 - within / between
    q = np.percentile(draws, 25, axis=1)
    keep = (draws <= q[:, None]).mean(0)
    better = (draws < draws[:, dcm_i][:, None]).mean(0)
    occ = (draws <= q[:, None]) & safe_inc[None, :]
    occ_n = occ.sum(1)
    _ex_dcm = safe_inc & (labs != 'dichloromethane')
    occ_nodcm = ((draws <= q[:, None]) & _ex_dcm[None, :]).sum(1)
    _entry = {l: float(occ[:, labs == l][:, 0].mean())
              for l in labs[safe_inc]}
    quad_inc = {
        'bar_log10': bar_inc,
        'published_count': int(((fid_c[has] <= cut0) & safe_inc).sum()),
        'quartile_cut_published': cut0,
        'quartile_cut_under_residual': {
            'median': float(np.median(q)),
            'ci95': [float(np.percentile(q, 2.5)),
                     float(np.percentile(q, 97.5))]},
        'frac_draws_nonempty': float((occ_n > 0).mean()),
        'mean_occupancy': float(occ_n.mean()),
        'frac_draws_nonempty_excluding_incumbent':
            float((occ_nodcm > 0).mean()),
        'entry_rate_by_compound': dict(sorted(_entry.items(),
                                              key=lambda kv: -kv[1])),
    }
    kept = {l: float(keep[labs == l][0]) for l in sl_names if (labs == l).any()}
    t1 = pd.Series(top1).value_counts(normalize=True)

    out = {
        'generated_by': 'h1_fidelity_stability.py', 'seed': SEED, 'n_draws': B,
        'note': ('전자구조 방법오차(GFN2->B3LYP)의 비아핀 잔차를 이 축에 '
                 '전파한 결과. 계통(아핀) 성분은 좌표별 표준화가 흡수한다.'),
        'calibration': {'n': int(len(m)), 'intercept': a.tolist(),
                        'slope': b.tolist(),
                        'residual_sd': R.std(0, ddof=2).tolist(),
                        'coords': EL},
        'affine_invariance': {'max_axis_shift': aff,
                              'spearman_before': float(rho0),
                              'spearman_after': float(rho_c)},
        'reliability_lambda': float(lam),
        'spearman_under_residual': {
            'median': float(np.median(rhos)),
            'ci95': [float(np.percentile(rhos, 2.5)),
                     float(np.percentile(rhos, 97.5))],
            'frac_abs_gt_published_ci': float((np.abs(rhos) > 0.199).mean())},
        # 2026-11-18. 상한을 하드코딩하고 있었다. 예치물에서 읽는다.
        'attenuation_corrected': {
            'rho': float(rho_c / np.sqrt(lam)),
            'ci_upper': float(_H['spearman_ci'][1] / np.sqrt(lam))},
        'quadrant_at_incumbent_bar_under_residual': quad_inc,
        'shortlist_quartile_retention': kept,
        'shortlist_mean_retention': float(np.mean(list(kept.values()))),
        'shortlist_p_more_faithful_than_dcm': {
            l: float(better[labs == l][0]) for l in kept},
        'n_candidates_beating_dcm_at_90pct': int((better >= 0.90).sum()),
        'most_frequent_top1': t1.head(5).round(4).to_dict(),
    }
    json.dump(out, open('outputs/h1_fidelity_stability.json', 'w'),
              ensure_ascii=False, indent=1)
    print(f'[lambda] 축 신뢰도 {lam:.4f}')
    print(f'[rho]    중앙 {np.median(rhos):+.4f}  '
          f'95% [{np.percentile(rhos, 2.5):+.3f}, '
          f'{np.percentile(rhos, 97.5):+.3f}]  '
          f"|rho|>0.199 비율 {out['spearman_under_residual']['frac_abs_gt_published_ci']:.4f}")
    print(f'[명단]   사분위 유지율 평균 {out["shortlist_mean_retention"]:.3f}  '
          f'DCM 보다 충실 90% 신뢰 후보 '
          f'{out["n_candidates_beating_dcm_at_90pct"]}종')
    for l, v in kept.items():
        print(f'   {l[:42]:44s} 유지 {v:.3f}  DCM우위 '
              f'{out["shortlist_p_more_faithful_than_dcm"][l]:.3f}')
    print('-> outputs/h1_fidelity_stability.json')


if __name__ == '__main__':
    main()
