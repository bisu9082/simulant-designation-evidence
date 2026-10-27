#!/usr/bin/env python3
"""전자공간을 B3LYP 궤도에너지로 갈아끼운 대조 (2026-10-25 신설).

왜. 심사 시늄레이션에서 "GFN2 의 omega 는 DFT 로 다시 재면서, 결론을 나르는
3좌표 전자공간(E_HOMO, E_LUMO, 쌍극자)은 GFN2 그대로 두었다"는 지적이 나왔다.
dft_vs_gfn2.csv 에 28종 B3LYP 궤도에너지가 이미 있으므로 돌려 볼 수 있다.

설계. 먼저 GFN2 로 정본값(유효매칭, 전수치환 p)이 재현되는지 확인한다. 그
관문을 통과하지 못하면 B3LYP 결과도 읽을 수 없다. 통과한 뒤 두 궤도에너지만
B3LYP 로 바꾸고 같은 계산을 반복한다.

한계를 먼저 적는다. 쌍극자는 여전히 GFN2 다. B3LYP 기하는 RDKit ETKDG +
MMFF94 이고 GFN2 기하가 아니다. 따라서 차이에는 방법과 기하가 함께 들어
있으며, 이 대조는 "전자공간 결론이 방법 선택에 얼마나 둔감한가"의 상한이
아니라 한 점 추정이다.

출력: outputs/dft_electronic_space.json
"""
import json
from itertools import combinations

import numpy as np
import pandas as pd

THRESHOLD = 2.0
PANEL = 'descriptor_cache_xtb.csv'
COLS = ['dipole_moment_D', 'E_LUMO_eV', 'E_HOMO_eV']


def dmat(df, cols):
    X = df[cols].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D


def matches(df, D, sim_cls=('simulant',)):
    """본문 2.5절의 매칭 2조건. 정규화거리 <= 2 AND 모든 비모사체보다 가까움."""
    nn_med = float(np.median(D.min(1)))
    kl = df.compound_class.tolist()
    nm = df.compound_name.tolist()
    i_ag = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sm = [i for i, k in enumerate(kl) if k in sim_cls]
    i_ot = [i for i, k in enumerate(kl)
            if k != 'real_agent' and k not in sim_cls]
    rows = []
    for i in i_ag:
        js = min(i_sm, key=lambda j: D[i, j])
        d_sim, d_oth = D[i, js], min(D[i, j] for j in i_ot)
        dn = d_sim / nn_med
        rows.append(dict(agent=nm[i], nearest=nm[js], d_norm=round(float(dn), 3),
                         match=bool(dn <= THRESHOLD and d_sim < d_oth)))
    return nn_med, rows


def perm_p(df, D):
    """C(17,3)=680 전수 치환. 통계량은 본문 헤드라인인 평균차."""
    kl = df.compound_class.to_numpy()
    i_real = np.flatnonzero(kl == 'real_agent')
    i_pool = np.flatnonzero(kl != 'real_agent')
    i_sim = set(np.flatnonzero(kl == 'simulant').tolist())
    pool_sim = tuple(sorted(list(i_pool).index(j) for j in i_sim))
    Dsub = D[np.ix_(i_real, i_pool)]
    allc = list(combinations(range(len(i_pool)), 3))
    assert len(allc) == 680
    diff = np.empty(680)
    ti = None
    for t, combo in enumerate(allc):
        m = np.zeros(len(i_pool), bool)
        m[list(combo)] = True
        diff[t] = (Dsub[:, m].min(1) - Dsub[:, ~m].min(1)).mean()
        if combo == pool_sim:
            ti = t
    obs = diff[ti]
    return float(obs), float((diff <= obs).mean())


def run(df, label):
    D = dmat(df, COLS)
    nn, rows = matches(df, D)
    obs, p = perm_p(df, D)
    n_match = sum(r['match'] for r in rows)
    print(f'  {label:10s} 전형최근접 {nn:.3f}  유효매칭 {n_match}/11  '
          f'관측 {obs:+.4f}  일측 p {p:.4f}')
    return dict(label=label, nn_median=nn, n_match=int(n_match),
                observed_diff=obs, p_one_sided=p, agents=rows)


def run_merged(df, label):
    """병합공간을 같은 교체로 다시 잰다 (2026-10-27 추가).

    전자공간만 보면 교체가 그 공간 안의 이야기로 끝난다. 그런데 병합공간은
    구조좌표와 전자좌표를 함께 쓰고 본문 결론을 나르는 p 를 낸다. 전자좌표
    세 개를 갈아끼우면 상관필터가 살리는 좌표 집합도 달라질 수 있으므로
    좌표 수도 함께 낸다.
    """
    from match_by_panel import merged_space, zdist
    D, nk = zdist(df, merged_space(df))
    np.fill_diagonal(D, np.inf)
    nn, rows = matches(df, D)
    obs, p = perm_p(df, D)
    n_match = sum(r['match'] for r in rows)
    print(f'  {label:14s} 좌표 {nk}개  전형최근접 {nn:.3f}  '
          f'유효매칭 {n_match}/11  관측 {obs:+.4f}  일측 p {p:.4f}')
    return dict(label=label, n_coords=int(nk), nn_median=nn,
                n_match=int(n_match), observed_diff=obs, p_one_sided=p,
                agents=rows)


def main():
    df = pd.read_csv(PANEL)
    dft = pd.read_csv('outputs/dft_vs_gfn2.csv').set_index('name')
    if not bool(dft.converged.all()):
        raise SystemExit('[중단] DFT 미수렴 항목이 있다. 대조를 만들 수 없다.')

    print('[1] 기준 재현 (GFN2 그대로)')
    base = run(df.copy(), 'GFN2')

    print('[2] 궤도에너지만 B3LYP/def2-SVPD 로 교체 (쌍극자는 GFN2 유지)')
    sw = df.copy()
    sw['E_HOMO_eV'] = sw.compound_name.map(dft.homo_eV).to_numpy(float)
    sw['E_LUMO_eV'] = sw.compound_name.map(dft.lumo_eV).to_numpy(float)
    assert sw[['E_HOMO_eV', 'E_LUMO_eV']].notna().all().all()
    swap = run(sw, 'B3LYP-2of3')

    print('[3] 세 좌표 전부 B3LYP (쌍극자도 교체)')
    dip = pd.read_csv('outputs/dft_dipole.csv').set_index('name')
    if not bool(dip.converged.all()):
        raise SystemExit('[중단] 쌍극자 계산에 미수렴 항목이 있다.')
    # 재현 관문. 쌍극자를 뽑은 SCF 가 궤도에너지 정본을 되살리지 못하면
    # 같은 계산이 아니다. 그 경우 쌍극자도 읽을 수 없다.
    for c in ('homo_eV', 'lumo_eV'):
        a = dft[c].reindex(dip.index).to_numpy(float)
        b = dip[c].to_numpy(float)
        bad = np.flatnonzero(np.abs(a - b) > 1e-6)
        if len(bad):
            raise SystemExit(
                f'[중단] 쌍극자 SCF 가 {c} 정본을 재현하지 못한다: '
                f'{[dip.index[i] for i in bad]}')
    print(f'  재현 관문 통과: 궤도에너지 {len(dip)}종 x 2 열 일치 (1e-6 eV)')
    fu = sw.copy()
    fu['dipole_moment_D'] = fu.compound_name.map(dip.dipole_D).to_numpy(float)
    assert fu['dipole_moment_D'].notna().all()
    full = run(fu, 'B3LYP-3of3')

    print('[4] 병합공간을 같은 교체로 재측정')
    m_gfn2 = run_merged(df.copy(), 'merged-GFN2')
    m_b3lyp = run_merged(fu.copy(), 'merged-B3LYP')

    out = dict(
        panel=PANEL, cols=COLS, threshold=THRESHOLD,
        gfn2=base, b3lyp=swap, b3lyp_full=full,
        merged_gfn2=m_gfn2, merged_b3lyp=m_b3lyp,
        caveat=('B3LYP 기하는 RDKit ETKDG+MMFF94 이고 GFN2 기하가 아니다. '
                '차이에 방법과 기하가 함께 들어 있다. b3lyp 배치는 쌍극자를 '
                'GFN2 로 남긴 2좌표 교체이고, b3lyp_full 은 세 좌표 전부 '
                'B3LYP 이다.'))
    json.dump(out, open('outputs/dft_electronic_space.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n[done] outputs/dft_electronic_space.json')


if __name__ == '__main__':
    main()
