#!/usr/bin/env python3
"""본문의 쌍별 거리 주장을 analysis_main.py 와 똑같은 공간에서 재계산한다 (2026-10-11).

에디터 지적 두 건을 대조하기 위해 만들었다.
  (1) "dichloromethane ... is beaten by ethyl acetate for phosgene in every
      space examined" — 전자 공간에서 참이 아닐 수 있다.
  (2) VX-DPGME 거리 3.13 이 구조 공간 값인지 병합 공간 값인지.
analysis_main.py 의 열 선택·표준화·상관 필터를 그대로 재현하므로,
여기서 나오는 값은 본문 표와 같은 정의의 값이다.
"""
import numpy as np
import pandas as pd

from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM,
                           PROTECTED, pick)

df = pd.read_csv('descriptor_cache_xtb.csv')
assert len(df) == 28, len(df)
num = [c for c in df.select_dtypes('number').columns
       if df[c].notna().sum() == 28 and c not in META_NUM]
SPACES = {
    'structural (A_rdkit)': [c for c in num if c not in ELECTRONIC_28],
    'electronic (B_electronic)': [c for c in num if c in ELECTRONIC_DIST],
    'merged (C_merged)': num,
}
names = df.compound_name.tolist()
idx = {n: i for i, n in enumerate(names)}


def dist(cols):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Xz = (X - X.mean(0)) / X.std(0, ddof=1)
    return np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1)), len(keep)


PAIRS = [('phosgene (CG)', 'dichloromethane'),
         ('phosgene (CG)', 'ethyl acetate'),
         ('VX', 'DPGME'),
         ('VX', 'DMMP')]

rows = []
for sp, cols in SPACES.items():
    D, k = dist(cols)
    print(f'\n=== {sp}  ({k}개 열) ===')
    for a, b in PAIRS:
        v = D[idx[a], idx[b]]
        print(f'  {a:<16} - {b:<16} {v:.4f}')
        rows.append(dict(space=sp, a=a, b=b, d=round(float(v), 4)))
    cg_dcm = D[idx['phosgene (CG)'], idx['dichloromethane']]
    cg_ea = D[idx['phosgene (CG)'], idx['ethyl acetate']]
    verdict = ('ethyl acetate 가 더 가깝다' if cg_ea < cg_dcm
               else 'dichloromethane 이 더 가깝다')
    print(f'  -> 포스겐 기준: {verdict}  (DCM {cg_dcm:.4f} vs EtOAc {cg_ea:.4f})')

pd.DataFrame(rows).to_csv('outputs/pair_claim_check.csv', index=False)
print('\n[done] outputs/pair_claim_check.csv')
