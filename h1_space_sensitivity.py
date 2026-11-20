#!/usr/bin/env python3
"""충실도 축을 구조공간에서 다시 그린다 (2026-11-18 신설).

왜 필요한가. 이 논문은 공간 선택을 선언해야 할 일곱 자유도 가운데 하나로
들고, 서술자 공간이 답을 바꾼다는 것을 세 군데에서 보인다. 그런데 서술자
거리를 유일하게 *건설적으로* 쓰는 축 — h1 의 fidelity, 즉 전자공간에서
최근접 진본 작용제까지의 표준화 거리 — 에는 그 민감도가 없다. 보고된
민감도는 배제 8건(rho 0.062 -> 0.030)뿐이고 공간 복제본이 없다.

여기서는 같은 풀, 같은 규칙으로 구조공간 fidelity 를 만들어 나란히 놓는다.
서술자는 build_descriptors.descriptors() 가 만들고, 좌표 선택은 본문과 같은
spaces.pick() 이 한다. 새 측정도 새 양자화학 계산도 쓰지 않는다.

[C-3] 이 스크립트는 h1_decoupling.py 의 예치물(outputs/h1_candidates.csv)과
28종 캐시에서 입력을 읽는다. fidelity 재구성이 예치값과 일치하는지 먼저
확인하고, 일치할 때만 구조공간 값을 보고한다.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

from build_descriptors import descriptors
from spaces import ELECTRONIC_28, META_NUM, PROTECTED, pick

EL = ['E_HOMO_eV', 'E_LUMO_eV', 'dipole_moment_D']


def nearest_agent(Z, is_agent):
    """표준화 좌표에서 최근접 진본 작용제까지의 유클리드 거리."""
    mu, sd = Z.mean(0), Z.std(0, ddof=1)
    Zz = (Z - mu) / sd
    A = Zz[is_agent]
    return np.array([np.sqrt(((A - z) ** 2).sum(1)).min() for z in Zz])


def main():
    pool = pd.read_csv('outputs/h1_pool.csv')
    cand = pd.read_csv('outputs/h1_candidates.csv')
    # 기준 작용제 집합은 compound_class == 'real_agent' 다 (11종).
    # is_cwa 는 skeleton 기반으로 후보에서 '빼는' 규칙이고(21종) 거리의
    # 기준이 아니다. 두 규칙을 섞으면 축이 달라진다.
    pool['is_agent'] = (pool.compound_class == 'real_agent').to_numpy()
    _fid = dict(zip(cand.label, cand.fidelity_dist))
    pool['fidelity_deposit'] = [_fid.get(l, np.nan) for l in pool.label]
    isag = pool.is_agent.to_numpy()
    print(f'[pool] 전체 {len(pool)}종 (작용제 {int(isag.sum())}, '
          f'후보 {int(pool.fidelity_deposit.notna().sum())})')

    # (1) 전자공간 재구성 — 예치값과 맞는지 먼저 본다
    fid_el = nearest_agent(pool[EL].to_numpy(float), isag)
    _has = pool.fidelity_deposit.notna().to_numpy()
    dev = np.nanmax(np.abs(fid_el[_has] - pool.fidelity_deposit[_has]))
    print(f'[check] 전자공간 fidelity 재구성 최대편차 {dev:.3e}')
    if dev > 1e-9:
        raise SystemExit('[중단] 예치값과 재구성이 다르다. 입력을 맞출 것')

    # (2) 구조공간 — 본문과 같은 생성기·같은 필터
    rows = [descriptors(r.label, r.smiles) for r in pool.itertuples()]
    D = pd.DataFrame(rows)
    n = len(D)
    num = [c for c in D.select_dtypes('number').columns
           if D[c].notna().sum() == n and c not in META_NUM
           and c not in ELECTRONIC_28]
    keep, _ = pick(D, num, PROTECTED)
    print(f'[struct] 투입 {len(num)}열 -> 유지 {len(keep)}열')
    fid_st = nearest_agent(D[keep].to_numpy(float), isag)

    out = {'generated_by': 'h1_space_sensitivity.py',
           'note': ('충실도 축의 공간 선택 민감도. 전자공간은 예치값 재현, '
                    '구조공간은 같은 풀·같은 필터로 새로 그린 것이다.'),
           'n_candidates': int(_has.sum()),
           'n_agents': int(isag.sum()),
           'structural_cols_in': len(num), 'structural_cols_kept': len(keep),
           'electronic_reconstruction_max_dev': float(dev)}

    y = pool.y_tox.to_numpy(float)
    cm = _has                       # 예치된 후보 207종만 보고 대상

    # (2b) 구조공간 상관이 어느 블록에 기대는가 (2026-11-18 추가).
    #
    # 왜. meth:h1stab 이 구조공간의 음계수를 승격하지 않는 사유로 '그것이
    # 의존하는 형태 블록은 meth:shape 가 대수적 중복이라고 보인 바로 그
    # 블록' 이라고 적는다. 그것은 수치 주장인데 수치가 예치되어 있지
    # 않았다. 여기서 블록을 빼 가며 잰다.
    _SHAPE = json.load(open('outputs/shape_reduction.json'))
    _derived = [c for c in _SHAPE['shape_derived_dropped'] if c in keep]
    _prim = [c for c in _SHAPE['shape_primitive'] if c in keep]
    _kap = [c for c in keep if c.lower().startswith('kappa')]
    _blocks = {
        'full': [],
        'shape_derived_removed': _derived,
        'shape_block_removed': _derived + _prim,
        'kappa_removed': _kap,
        'shape_and_kappa_removed': _derived + _prim + _kap,
    }
    _ok_t = ~np.isnan(y[cm])
    block_out = {}
    for _tag, _drop in _blocks.items():
        _cols = [c for c in keep if c not in set(_drop)]
        _f = nearest_agent(D[_cols].to_numpy(float), isag)[cm]
        _r, _p = stats.spearmanr(_f[_ok_t], y[cm][_ok_t])
        block_out[_tag] = {'n_cols': len(_cols), 'dropped': sorted(_drop),
                           'spearman_rho': float(_r), 'spearman_p': float(_p)}
        print(f'[block] {_tag:26s} {len(_cols):3d}열  rho {_r:+.4f}  p {_p:.4f}')
    out['structural_block_dependence'] = block_out
    for tag, fid in (('electronic', fid_el), ('structural', fid_st)):
        f = fid[cm]
        t = y[cm]
        ok = ~np.isnan(t)
        rho, p = stats.spearmanr(f[ok], t[ok])
        cut = float(np.percentile(f, 25))
        # 안전 경계는 본문과 같다: 현행 모사체 dichloromethane 의
        # log10(AEGL-3 10분). 백분위가 아니다.
        _tx = pd.read_csv('outputs/toxicity_aegl.csv')
        _d = _tx.loc[_tx.compound_name == 'dichloromethane',
                     'AEGL3_10min_mgm3']
        tox_cut = float(np.log10(float(_d.iloc[0])))
        quad = pool.label[cm][(f <= cut) & (t >= tox_cut)].tolist()
        order = np.argsort(f)
        out[tag] = {
            'spearman_rho': float(rho), 'spearman_p': float(p),
            'n': int(ok.sum()),
            'quartile_cut': cut, 'tox_cut': tox_cut,
            'most_faithful': pool.label[cm].iloc[int(order[0])],
            'most_faithful_dist': float(f[order[0]]),
            'dichloromethane_dist': float(
                f[pool.label[cm].to_numpy() == 'dichloromethane'][0]),
            'faithful_and_safe': quad,
            'top8': pool.label[cm].iloc[order[:8]].tolist(),
        }
        print(f'[{tag:>10}] rho {rho:+.4f} (p={p:.4f})  사분위 경계 {cut:.3f}  '
              f'최충실 {out[tag]["most_faithful"]} '
              f'({out[tag]["most_faithful_dist"]:.3f})  '
              f'사분면 {len(quad)}종')

    # 두 축의 순위 일치
    rho_ax, p_ax = stats.spearmanr(fid_el[cm], fid_st[cm])
    out['axis_agreement'] = {'spearman_rho': float(rho_ax),
                             'spearman_p': float(p_ax)}
    print(f'[agree] 두 fidelity 축의 순위 일치 rho {rho_ax:+.4f}')

    # 논문이 완화 경계에서 이름을 적은 여덟 종(h1_summary 의 relaxed_10)이
    # 구조공간에서 어디에 있는가. 상위 8종이 아니라 그 명단이다.
    ranks = stats.rankdata(fid_st[cm]) / cm.sum() * 100
    lab = pool.label[cm].to_numpy()
    _sl = json.load(open('outputs/h1_summary.json'))[
        'safety_tiers']['relaxed_10']['labels']
    out['shortlist'] = _sl
    out['shortlist_structural_percentile'] = {
        l: float(ranks[lab == l][0]) for l in _sl if (lab == l).any()}
    _v = list(out['shortlist_structural_percentile'].values())
    out['shortlist_structural_percentile_range'] = [min(_v), max(_v)]
    print(f'[cross] 논문 명단 {len(_sl)}종의 구조공간 백분위 '
          f'{min(_v):.1f}~{max(_v):.1f}: '
          + ', '.join(f'{k} {v:.1f}' for k, v in
                      out['shortlist_structural_percentile'].items()))

    json.dump(out, open('outputs/h1_space_sensitivity.json', 'w'),
              ensure_ascii=False, indent=1)
    print('-> outputs/h1_space_sensitivity.json')


if __name__ == '__main__':
    main()
