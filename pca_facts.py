#!/usr/bin/env python3
"""공간별 주성분 분산 (2026-10-24 신설).

왜 필요한가. 본문 4.1절의 36.5 / 70.4 / 45.9 퍼센트는 파이프라인이 만들지 않는
파일(outputs/raw_results.csv, make_step4_outputs.py 는 run_analysis.sh 에 없다)에서
왔고, 45.9 는 현재 어떤 패널·공간과도 맞지 않았다. 검증기가 PCA 분산을 사실로
들고 있지 않아 240건 대조가 이것을 덮지 못했다 (2026-10-24 페르소나 심사 적발).

출력: outputs/pca_facts.json
"""
import json

import pandas as pd
from sklearn.decomposition import PCA

from spaces import ELECTRONIC_28, ELECTRONIC_DIST, META_NUM, P_LOCAL, PROTECTED, pick

PANELS = {'panel28': 'descriptor_cache_xtb.csv',
          'panel35': 'descriptor_cache_xtb_35.csv'}


def build(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM and c not in P_LOCAL]
    struct = [c for c in num if c not in ELECTRONIC_28]
    A, _ = pick(df, struct, PROTECTED)
    B = [c for c in ELECTRONIC_DIST if c in num]
    C, _ = pick(df, struct + B, PROTECTED)
    return {'A_rdkit': A, 'B_electronic': B, 'C_merged': C}


def main():
    out = {}
    for tag, path in PANELS.items():
        df = pd.read_csv(path)
        out[tag] = {}
        for name, cols in build(df).items():
            X = df[cols].to_numpy(float)
            Z = (X - X.mean(0)) / X.std(0, ddof=1)
            v = PCA().fit(Z).explained_variance_ratio_
            out[tag][name] = {'k': len(cols),
                              'pc1': round(float(v[0]), 5),
                              'pc1_pc2': round(float(v[0] + v[1]), 5)}
            print(f'  {tag} {name:14s} k={len(cols):2d} '
                  f'PC1={v[0]*100:.1f}% PC1+PC2={(v[0]+v[1])*100:.1f}%')
    with open('outputs/pca_facts.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print('[pca] outputs/pca_facts.json')


if __name__ == '__main__':
    main()
