#!/usr/bin/env python3
"""최근접 모사체 거리를 '그 공간의 전형적 최근접거리'로 정규화.
   거리 12에서 이긴 최근접이웃은 매칭이 아니다 — 이걸 수치로 만든다."""
import json
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
S = json.load(open('outputs/experiment_summary.json'))
import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument('--data', default='descriptor_cache_xtb.csv')
_ap.add_argument('--table', default='outputs/distance_table.csv')
_ap.add_argument('--suffix', default='')
_A = _ap.parse_args()
df = pd.read_csv(_A.data)
T = pd.read_csv(_A.table)
rows = []
for cond in S['conditions']:
    keep = cond['kept']
    X = df[keep].to_numpy(float)
    Xz = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Xz[:, None, :]-Xz[None, :, :])**2).sum(-1)); np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))            # 패널 전체의 최근접거리 중앙값
    t = T[T.condition == cond['condition']].copy()
    t['nn_median_space'] = nn_med
    t['d_sim_norm'] = t.d_simulant / nn_med
    # 본문 2.4절의 매칭 정의는 두 조건의 논리곱이다: 정규화 거리가 전형
    # 최근접거리의 2배 이내 AND 모든 간섭물보다 가까움. 종전에는 이 열이
    # 앞쪽 조건만 담아 배포 CSV 가 본문 정의와 어긋났다 (2026-10-12 교정).
    t['within_2x'] = t.d_sim_norm <= 2.0
    t['is_match'] = t.within_2x & ~t.interferent_closer
    rows.append(t)
    print(f"[{cond['condition']}] 전형 최근접거리(중앙값) = {nn_med:.3f}")
R = pd.concat(rows, ignore_index=True)
R.to_csv(f'outputs/distance_table_normalized{_A.suffix}.csv', index=False)
print()
for c in ['A_rdkit','B_electronic','C_merged']:
    t = R[R.condition == c].sort_values('d_sim_norm')
    n_match = int(t.within_2x.sum())
    n_win = int(t.is_match.sum())
    print(f'--- {c} ---  모사체가 전형거리 2배 이내: {n_match}/11  |  '
          f'그중 간섭물까지 이긴 경우: {n_win}/11')
    print(t[['real_agent','nearest_simulant','d_sim_norm','within_2x','interferent_closer','is_match']]
          .to_string(index=False, float_format=lambda v: f'{v:.2f}'))
    print()
