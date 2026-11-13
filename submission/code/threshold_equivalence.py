#!/usr/bin/env python3
"""매칭 역치가 공간 사이에서 등가인지 검사한다.

왜 필요한가. 본문은 세 공간의 매칭 수를 나란히 놓고 전자공간을
'the exception and the informative one' 이라 부른다. 그 비교가 성립하려면
'정규화거리 <= 2' 라는 역치가 세 공간에서 같은 엄격도여야 한다.

median NN 으로 나누는 것은 국소 척도를 맞추지만 역치가 거리분포에서
차지하는 자리까지 맞추지는 못한다. 그래서 역치가 전 쌍거리의 몇 퍼센트를
포함하는지를 공간마다 센다. 그 값이 다르면 매칭 수의 공간 간 대비는
화학이 아니라 차원의 산물을 일부 담는다.

함께, 공간의 명목 좌표 수가 실제 자유도와 얼마나 다른지 적는다.
p > n 이므로 좌표 수는 선형독립 방향 수를 넘을 수 있다.

GFN2 LUMO 의 방법오차도 같은 파일에 담는다. 본문·SI 는 GFN2 와 B3LYP
LUMO 의 전역 순위상관 0.898 을 그 축의 안심 근거로 든다. 그런데 최근접
판정은 국소 성질이고, 편차의 폭이 좌표 자체의 화학적 폭보다 크면 전역
단조성은 국소 순서를 보증하지 못한다. 두 폭을 나란히 적는다.

[C-3] 이 파일은 예치 캐시에서만 읽는다. 새 양자계산을 하지 않는다.
"""
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, ELECTRONIC_DIST, META_NUM, OMEGA_DERIVED
from spaces import PROTECTED, pick

PANELS = [('panel28_incumbent3', 'descriptor_cache_xtb.csv'),
          ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv')]
THRESHOLD = 2.0


def spaces_of(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return {'structural': [c for c in num if c not in ELECTRONIC_28],
            'electronic3': [c for c in num if c in ELECTRONIC_DIST],
            'merged': [c for c in num if c != OMEGA_DERIVED]}


def one_space(df, cols):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    iu = np.triu_indices(len(df), 1)
    pair = D[iu]
    cut = THRESHOLD * nn_med
    # 유효차원: 고윳값 참여비 exp(H) 로 센다. 좌표 수가 아니라 방향 수다.
    C = np.cov(Z, rowvar=False)
    ev = np.clip(np.linalg.eigvalsh(C), 0, None)
    p = ev / ev.sum() if ev.sum() > 0 else ev
    p = p[p > 0]
    eff = float(np.exp(-(p * np.log(p)).sum()))
    return {
        'n_coords': len(keep),
        'matrix_rank': int(np.linalg.matrix_rank(Z - Z.mean(0))),
        'effective_dim_participation': round(eff, 3),
        'nn_median': round(nn_med, 4),
        'threshold_distance': round(cut, 4),
        'pairs_total': int(pair.size),
        'pairs_within_threshold': int((pair <= cut).sum()),
        'pct_pairs_within_threshold': round(
            100.0 * float((pair <= cut).mean()), 2),
        'relative_contrast': round(
            float((pair.max() - pair.min()) / pair.min()), 2),
    }


def lumo_method_error():
    """GFN2 LUMO 의 방법오차를 그 좌표의 화학적 폭에 대해 잰다."""
    df = pd.read_csv('outputs/dft_vs_gfn2.csv')
    g = df['lumo_eV'].to_numpy(float)          # B3LYP/def2-SVPD
    x = df['E_LUMO_eV'].to_numpy(float)        # GFN2-xTB (예치 캐시)
    ok = np.isfinite(g) & np.isfinite(x)
    g, x, nm = g[ok], x[ok], df['compound_name'].to_numpy()[ok]
    d = x - g
    rg = float(g.max() - g.min())
    rd = float(d.max() - d.min())
    rs = float(pd.Series(g).rank().corr(pd.Series(x).rank(), method='spearman'))
    shift = (pd.Series(x).rank() - pd.Series(g).rank()).to_numpy()
    j = int(np.argmax(np.abs(shift)))
    # 2026-11-13 (GATE 8 10차). 원시 eV 비는 이 공간에 작용하는 양이 아니다.
    # 이 논문의 모든 거리는 패널 내 z표준화 좌표 위에 있으므로, 축 하나에
    # 걸리는 평행이동과 축척은 거리에 영향이 없다. 두 LUMO 열을 각각
    # 표준화한 뒤의 화합물별 변위가 실제로 거리를 움직이는 양이다.
    zg = (g - g.mean()) / g.std(ddof=1)
    zx = (x - x.mean()) / x.std(ddof=1)
    disp = np.abs(zx - zg)
    return {
        'n_compounds': int(ok.sum()),
        'standardised_displacement_rmsd': round(
            float(np.sqrt((disp ** 2).mean())), 4),
        'standardised_displacement_max': round(float(disp.max()), 4),
        'standardised_displacement_note': (
            'In coordinate standard deviations, the unit the distances are '
            'actually computed in. Compare the median nearest-neighbour '
            'separation of the electronic space in the same units.'),
        'b3lyp_lumo_range_eV': round(rg, 4),
        'gfn2_minus_b3lyp_min_eV': round(float(d.min()), 4),
        'gfn2_minus_b3lyp_max_eV': round(float(d.max()), 4),
        'deviation_span_eV': round(rd, 4),
        'deviation_span_over_chemical_range': round(rd / rg, 3),
        'spearman_lumo': round(rs, 4),
        'largest_rank_shift': int(abs(shift[j])),
        'largest_rank_shift_compound': str(nm[j]),
    }


def main():
    out = {'generated_by': 'threshold_equivalence.py',
           'note': ('매칭 역치가 공간 사이에서 등가인지, 그리고 GFN2 LUMO 의 '
                    '방법오차가 그 좌표의 화학적 폭에 비해 얼마나 큰지.'),
           'threshold_rule': 'normalised distance <= 2 x median nearest neighbour',
           'panels': {}}
    for label, path in PANELS:
        df = pd.read_csv(path)
        out['panels'][label] = {k: one_space(df, v)
                                for k, v in spaces_of(df).items()}
    out['gfn2_lumo_method_error'] = lumo_method_error()

    with open('outputs/threshold_equivalence.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)

    for label, sp in out['panels'].items():
        print(f'  {label}')
        for k, v in sp.items():
            print(f"    {k:<12} 좌표 {v['n_coords']:>2} 랭크 {v['matrix_rank']:>2} "
                  f"유효 {v['effective_dim_participation']:>5} | "
                  f"역치 {v['threshold_distance']:>7} = 전 쌍거리의 "
                  f"{v['pct_pairs_within_threshold']:>5}% | "
                  f"상대대비 {v['relative_contrast']}")
    e = out['gfn2_lumo_method_error']
    print(f"\n  GFN2 LUMO 방법오차: {e['gfn2_minus_b3lyp_min_eV']:+.2f} ~ "
          f"{e['gfn2_minus_b3lyp_max_eV']:+.2f} eV (폭 {e['deviation_span_eV']:.2f}) "
          f"대 B3LYP LUMO 전 범위 {e['b3lyp_lumo_range_eV']:.2f} eV "
          f"-> {e['deviation_span_over_chemical_range']:.2f}배")
    print(f"  전역 순위상관 {e['spearman_lumo']:.3f}, 최대 순위이동 "
          f"{e['largest_rank_shift']}위 ({e['largest_rank_shift_compound']})")
    print('  -> outputs/threshold_equivalence.json')


if __name__ == '__main__':
    main()
