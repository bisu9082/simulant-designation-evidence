#!/usr/bin/env python3
"""
STEP 4 / simulant — GFN2-xTB 전자 서술자 단독 분석
입력: xtb_log_electronic.csv  (xtb_run.log 에 실제 출력된 값만 전사)
출력: electronic_space_summary.json, electronic_distance_table.csv, Fig_electronic.png

주의: q_P / fukui_f_plus_P 는 인 함유 6종에만 존재한다. 전역 PCA 에는
      28종 모두가 가진 4개 열(HOMO, LUMO, omega, dipole)만 사용한다.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

SEED = 42
N_BOOT = 10000
rng = np.random.default_rng(SEED)

# omega 는 HOMO/LUMO 의 결정론적 함수이므로 거리/PCA 에서 제외한다 (2026-10-13, M3/C1).
# 표 보고용으로는 아래 P 국소 표에 남는다.
GLOBAL_COLS = ['E_HOMO_eV', 'E_LUMO_eV', 'dipole_moment_D']
OMEGA_DERIVED = 'electrophilicity_index_eV'
P_COLS = ['q_P', 'fukui_f_plus_P']

# 2026-10-12. 이 스크립트만 xtb_log_electronic.csv 라는 구판 파일을 읽고
# 있었다. 그 결과 같은 양이 두 산출물에서 갈렸다 (간섭물 우세 8 대 7,
# 쌍체 t p=0.0034 대 0.0354). 정본 서술자 표로 통일한다.
import argparse as _ap_mod
_ap = _ap_mod.ArgumentParser()
_ap.add_argument('--data', default='descriptor_cache_xtb.csv')
_ap.add_argument('--suffix', default='')
_A = _ap.parse_args()
df = pd.read_csv(_A.data)
df = df.rename(columns={'compound_name': 'name', 'compound_class': 'klass'})
df['gap_eV'] = df['E_LUMO_eV'] - df['E_HOMO_eV']
print(f'[load] {len(df)} compounds  |  '
      f"real={sum(df.klass=='real_agent')} sim={sum(df.klass=='simulant')} "
      f"intf={sum(df.klass=='interferent')}")
print(f'[load] P-bearing (q_P present): {df.q_P.notna().sum()}  -> '
      f"{', '.join(df.loc[df.q_P.notna(),'name'])}")
assert len(df) >= 28 and df.q_P.notna().sum() >= 6

# ---------------------------------------------------------------- 1. 전역 PCA
X = df[GLOBAL_COLS].to_numpy(float)
Xz = (X - X.mean(0)) / X.std(0, ddof=1)
pca = PCA(n_components=len(GLOBAL_COLS), random_state=SEED).fit(Xz)
Y = pca.transform(Xz)
evr = pca.explained_variance_ratio_
print(f'\n[PCA] explained variance ratio = '
      f"{', '.join(f'{v:.3f}' for v in evr)}  (PC1+PC2 = {evr[:2].sum():.3f})")
for i in range(2):
    load = sorted(zip(GLOBAL_COLS, pca.components_[i]), key=lambda t: -abs(t[1]))
    print(f'  PC{i+1} loadings: ' + ' | '.join(f'{n} {v:+.3f}' for n, v in load))

# ------------------------------------------- 2. 표준화 유클리드 거리 (전자공간)
D = np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1))
np.fill_diagonal(D, np.inf)
name = df['name'].tolist()
klass = df['klass'].tolist()
i_sim = [i for i, k in enumerate(klass) if k == 'simulant']
i_int = [i for i, k in enumerate(klass) if k == 'interferent']
i_real = [i for i, k in enumerate(klass) if k == 'real_agent']

rows = []
for i in i_real:
    js = i_sim[int(np.argmin(D[i, i_sim]))]
    ji = i_int[int(np.argmin(D[i, i_int]))]
    rows.append(dict(real_agent=name[i],
                     nearest_simulant=name[js], d_simulant=D[i, js],
                     nearest_interferent=name[ji], d_interferent=D[i, ji],
                     interferent_closer=bool(D[i, ji] < D[i, js])))
T = pd.DataFrame(rows).sort_values('d_simulant').reset_index(drop=True)
T.to_csv(f'outputs/electronic_distance_table{_A.suffix}.csv', index=False)
print(f'\n[dist] 간섭물이 더 가까운 진본 작용제: '
      f'{int(T.interferent_closer.sum())} / {len(T)}')
print(T[['real_agent', 'nearest_simulant', 'd_simulant',
         'nearest_interferent', 'd_interferent', 'interferent_closer']]
      .to_string(index=False, float_format=lambda v: f'{v:.3f}'))

# ------------------------------------------------- 3. 쌍체 검정 + 효과크기 + CI
a, b = T.d_simulant.to_numpy(), T.d_interferent.to_numpy()
d_pair = a - b                                      # >0 이면 간섭물이 더 가깝다
sw = stats.shapiro(d_pair)
if sw.pvalue > 0.05:
    test, res = 't-test (paired)', stats.ttest_rel(a, b)
else:
    test, res = 'Wilcoxon signed-rank', stats.wilcoxon(a, b)
cohen_dz = d_pair.mean() / d_pair.std(ddof=1)
boot = np.array([rng.choice(d_pair, d_pair.size, replace=True).mean()
                 for _ in range(N_BOOT)])
ci = np.percentile(boot, [2.5, 97.5])
print(f'\n[stat] Shapiro-Wilk W={sw.statistic:.4f} p={sw.pvalue:.4f} -> {test}')
print(f'[stat] mean(d_sim - d_intf) = {d_pair.mean():+.4f}  '
      f"95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]  (bootstrap n={N_BOOT})")
print(f'[stat] {test}: stat={res.statistic:.4f} p={res.pvalue:.4f}  '
      f"Cohen's dz={cohen_dz:+.4f}")

# ------------------------------------------------ 4. 인 중심 국소 반응성 (6종)
P = df[df.q_P.notna()][['name', 'klass'] + P_COLS + ['electrophilicity_index_eV']]
P = P.sort_values('fukui_f_plus_P', ascending=False).reset_index(drop=True)
print('\n[P-local] 인 원자 국소 친전자성 (f+(P) 내림차순)')
print(P.to_string(index=False, float_format=lambda v: f'{v:.5f}'))
gb = float(P.loc[P.name == 'sarin (GB)', 'fukui_f_plus_P'].iloc[0])
dm = float(P.loc[P.name == 'DMMP', 'fukui_f_plus_P'].iloc[0])
vx = float(P.loc[P.name == 'VX', 'fukui_f_plus_P'].iloc[0])
print(f'\n[P-local] DMMP / GB = {dm/gb:.3f}   VX / GB = {vx/gb:.3f}   '
      f'DMMP / VX = {dm/vx:.3f}')
g = P[P.name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])]['fukui_f_plus_P']
print(f'[P-local] G-series(GB,GD,GF) f+(P) = {g.mean():.5f} +/- {g.std(ddof=1):.5f}'
      f'  -> DMMP z = {(dm - g.mean())/g.std(ddof=1):+.2f}')

# --------------------------------------------------------------- 5. 요약 JSON
out = dict(
    n_compounds=int(len(df)), n_phosphorus=int(df.q_P.notna().sum()),
    descriptors_global=GLOBAL_COLS, seed=SEED, bootstrap_samples=N_BOOT,
    pca=dict(explained_variance_ratio=[float(v) for v in evr],
             pc1_pc2_cumulative=float(evr[:2].sum()),
             loadings={f'PC{i+1}': {c: float(v) for c, v in
                                    zip(GLOBAL_COLS, pca.components_[i])}
                       for i in range(2)}),
    distance_test=dict(test_used=test, statistic=float(res.statistic),
                       p_value=float(res.pvalue), cohens_dz=float(cohen_dz),
                       ci_95=[float(ci[0]), float(ci[1])],
                       shapiro_W=float(sw.statistic), shapiro_p=float(sw.pvalue),
                       mean_diff_sim_minus_intf=float(d_pair.mean()),
                       n_interferent_closer=int(T.interferent_closer.sum()),
                       n_real_agents=int(len(T))),
    p_local=dict(f_plus_P={r['name']: float(r['fukui_f_plus_P'])
                           for _, r in P.iterrows()},
                 q_P={r['name']: float(r['q_P']) for _, r in P.iterrows()},
                 DMMP_over_GB=float(dm / gb), VX_over_GB=float(vx / gb),
                 G_series_mean=float(g.mean()), G_series_sd=float(g.std(ddof=1)),
                 DMMP_z_vs_G_series=float((dm - g.mean()) / g.std(ddof=1))),
)
with open(f'outputs/electronic_space_summary{_A.suffix}.json', 'w') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
np.save('outputs/_pca_scores.npy', Y)
df.to_csv('outputs/_electronic_with_gap.csv', index=False)
print('\n[done] outputs/electronic_space_summary.json, electronic_distance_table.csv')
