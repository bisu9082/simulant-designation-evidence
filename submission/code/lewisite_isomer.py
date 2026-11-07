#!/usr/bin/env python3
"""루이사이트 1 의 기하이성질체 민감도 (2026-11-04 신설, 작업 E2).

왜. 원고는 루이사이트 1 을 E(trans) 이성질체로 계산한다
(input_smiles = Cl/C=C/[As](Cl)Cl, RDKit 이 STEREOE 로 읽는다). 그런데
본문이 그 선택을 밝히지 않았다. 둘째 런의 재료는 순도 67% 이고 잔부가
무엇인지 원논문이 적지 않으므로, 이성질체 선택이 결과를 움직이는지 재서
적는다. 주장하지 않고 측정한다.

변종: Z(cis) 이성질체 Cl/C=C\\[As](Cl)Cl 하나만 바꾸고 나머지 패널은
그대로 둔다. 35종 패널의 병합공간에서 루이사이트 행이 어떻게 움직이는지
본다.

산출: outputs/lewisite_isomer.json
"""
import json

import numpy as np
import pandas as pd

import build_descriptors as B
from match_by_panel import merged_space
from spaces import ELECTRONIC_28, META_NUM, PROTECTED, pick

PANEL = 'descriptor_cache_xtb_35.csv'
TARGET = 'lewisite 1 (L)'
Z_SMILES = r'Cl/C=C\[As](Cl)Cl'


def structural_keep(df):
    """구조공간. 공간 정의를 복제하지 않고 정본(match_by_panel.merged_space)
    에서 전자 좌표만 뺀다. 전자 서술자는 xTB 바이너리가 있어야 다시 만들 수
    있는데 이 환경에 없으므로, 전부 재생성되는 구조공간만 본다."""
    cols = [c for c in merged_space(df) if c not in ELECTRONIC_28]
    keep, _ = pick(df, cols, PROTECTED)
    return keep


def lewisite_row(df):
    keep = structural_keep(df)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    names = df.compound_name.tolist()
    kl = df.compound_class.tolist()
    i = names.index(TARGET)
    sim = [j for j, k in enumerate(kl) if k in ('simulant', 'literature_simulant')]
    other = [j for j, k in enumerate(kl) if k != 'real_agent' and j not in sim]
    js = min(sim, key=lambda j: D[i, j])
    jo = min(other, key=lambda j: D[i, j]) if other else None
    d_sim = float(D[i, js])
    return {'space': 'structural', 'n_coords': len(keep), 'nn_median': nn_med,
            'nearest_simulant': names[js], 'd_simulant': d_sim,
            'd_norm': d_sim / nn_med,
            'nearest_nonsimulant': names[jo] if jo is not None else None,
            'd_nonsimulant': float(D[i, jo]) if jo is not None else None,
            'is_match': bool(d_sim / nn_med <= 2.0
                             and (jo is None or D[i, js] <= D[i, jo]))}


def main():
    df = pd.read_csv(PANEL)
    base = lewisite_row(df)
    row = df[df.compound_name == TARGET].iloc[0]
    print('=' * 72)
    print('[E2] 루이사이트 1 기하이성질체 민감도 (구조공간)')
    print('=' * 72)
    print(f"  정본 input_smiles : {row.input_smiles}")
    z = B.descriptors(TARGET, Z_SMILES, row.compound_class, row.agent_series)
    alt = df.copy()
    miss = [c for c in alt.columns if c not in z]
    for c in alt.columns:
        if c in z:
            alt.loc[alt.compound_name == TARGET, c] = z[c]
    _el = [c for c in miss if c in ELECTRONIC_28]
    print(f"  Z 이성질체 재계산. 재생성되지 않은 열 {len(miss)}개 중 "
          f"전자 좌표 {len(_el)}개 — 이 때문에 병합공간은 보지 않는다")
    zz = lewisite_row(alt)
    print()
    for k in ['n_coords', 'nn_median', 'nearest_simulant', 'd_simulant',
              'd_norm', 'nearest_nonsimulant', 'is_match']:
        a, b = base[k], zz[k]
        same = '' if a == b else '   <- 바뀜'
        fa = f'{a:.4f}' if isinstance(a, float) else str(a)
        fb = f'{b:.4f}' if isinstance(b, float) else str(b)
        print(f'  {k:22s} E(정본) {fa:22s} Z {fb:22s}{same}')
    out = {'generated_by': 'lewisite_isomer.py', 'panel': PANEL,
           'scope': ('structural space only; the electronic descriptors need '
                     'the xTB binary, which was not available, so they were '
                     'not regenerated for the Z isomer and the merged space '
                     'is not reported'),
           'canonical_input_smiles': row.input_smiles,
           'canonical_stereo': 'E (trans)',
           'alternative_input_smiles': Z_SMILES,
           'columns_not_regenerated': miss,
           'E_trans': base, 'Z_cis': zz}
    with open('outputs/lewisite_isomer.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('\n-> outputs/lewisite_isomer.json')


if __name__ == '__main__':
    main()
