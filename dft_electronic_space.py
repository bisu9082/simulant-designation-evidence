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
    swap = run(sw, 'B3LYP')

    out = dict(
        panel=PANEL, cols=COLS, threshold=THRESHOLD,
        gfn2=base, b3lyp=swap,
        caveat=('쌍극자는 GFN2 유지. B3LYP 기하는 RDKit ETKDG+MMFF94 이고 '
                'GFN2 기하가 아니다. 차이에 방법과 기하가 함께 들어 있다.'))
    json.dump(out, open('outputs/dft_electronic_space.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n[done] outputs/dft_electronic_space.json')


if __name__ == '__main__':
    main()
