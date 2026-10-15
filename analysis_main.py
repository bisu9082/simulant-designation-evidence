#!/usr/bin/env python3
"""
simulant / STEP 4 — 병합 서술자 공간 분석 (RDKit + GFN2-xTB)

조건 3개 ablation:
  A  rdkit      : RDKit 서술자만
  B  electronic : GFN2 전자 서술자만 (28종 전부가 가진 4개)
  C  merged     : A + B

각 조건에서 분산필터 → 상관필터(|r|>0.95, PROTECTED 우선) → z표준화 → PCA + 거리.
H1 검정: 진본 작용제별 최근접 모사체 vs 최근접 간섭물 (쌍체, bootstrap CI, 효과크기).

클래스 라벨은 descriptor_cache_xtb.csv 의 compound_class 열을 정본으로 사용한다.
"""
import argparse, json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

SEED, N_BOOT, R_CUT = 42, 10000, 0.95
rng = np.random.default_rng(SEED)

# 공간 정의와 pick() 은 spaces.py 가 정본이다 (2026-10-13 분리).
# 이 파일은 스크립트이므로 임포트하면 argparse 와 분석 전체가 돌아간다.
# 다른 스크립트는 반드시 spaces 에서 가져올 것.
from spaces import (ELECTRONIC_28, ELECTRONIC_DIST, META_NUM,  # noqa: F401
                    OMEGA_DERIVED, P_LOCAL, PROTECTED, R_CUT, pick)


def analyse(df, cols, tag, protected=PROTECTED):
    keep, dropped = pick(df, cols, protected)
    X = df[keep].to_numpy(float)
    Xz = (X - X.mean(0)) / X.std(0, ddof=1)
    pca = PCA(n_components=min(len(keep), len(df)), random_state=SEED).fit(Xz)
    evr = pca.explained_variance_ratio_
    load1 = sorted(zip(keep, pca.components_[0]), key=lambda t: -abs(t[1]))[:6]

    D = np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    name, klass = df.compound_name.tolist(), df.compound_class.tolist()
    i_s = [i for i, k in enumerate(klass) if k in SIM_CLASSES]
    i_i = [i for i, k in enumerate(klass) if k == 'interferent']
    i_r = [i for i, k in enumerate(klass) if k == 'real_agent']

    rows = []
    for i in i_r:
        js, ji = i_s[int(np.argmin(D[i, i_s]))], i_i[int(np.argmin(D[i, i_i]))]
        rows.append(dict(condition=tag, real_agent=name[i],
                         nearest_simulant=name[js], d_simulant=D[i, js],
                         nearest_interferent=name[ji], d_interferent=D[i, ji],
                         interferent_closer=bool(D[i, ji] < D[i, js])))
    T = pd.DataFrame(rows)
    a, b = T.d_simulant.to_numpy(), T.d_interferent.to_numpy()
    dp = a - b
    sw = stats.shapiro(dp)
    if sw.pvalue > 0.05:
        test, res = 't-test (paired)', stats.ttest_rel(a, b)
    else:
        test, res = 'Wilcoxon signed-rank', stats.wilcoxon(a, b)
    boot = np.array([rng.choice(dp, dp.size, replace=True).mean()
                     for _ in range(N_BOOT)])
    ci = np.percentile(boot, [2.5, 97.5])

    print(f'\n=== [{tag}] 서술자 {len(cols)} -> 보존 {len(keep)} '
          f'(상관 제거 {len(dropped)}) ===')
    print(f'  PC1 {evr[0]*100:5.1f}%  PC2 {evr[1]*100:5.1f}%  '
          f'누적 {evr[:2].sum()*100:5.1f}%')
    print('  PC1 적재 상위: ' + ' | '.join(f'{n} {v:+.3f}' for n, v in load1))
    print(f'  간섭물이 더 가까움: {int(T.interferent_closer.sum())} / {len(T)}')
    print(f'  Shapiro p={sw.pvalue:.4f} -> {test}: stat={res.statistic:.4f} '
          f'p={res.pvalue:.4f}  dz={dp.mean()/dp.std(ddof=1):+.4f}')
    print(f'  mean(d_sim - d_intf)={dp.mean():+.4f}  '
          f'95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]')

    return dict(condition=tag, n_in=len(cols), n_kept=len(keep),
                kept=keep, dropped=[list(x) for x in dropped],
                evr=[float(v) for v in evr[:5]],
                pc1_pc2=float(evr[:2].sum()),
                pc1_loadings={n: float(v) for n, v in load1},
                n_interferent_closer=int(T.interferent_closer.sum()),
                n_real=len(T), test_used=test,
                statistic=float(res.statistic), p_value=float(res.pvalue),
                cohens_dz=float(dp.mean() / dp.std(ddof=1)),
                mean_diff=float(dp.mean()), ci_95=[float(ci[0]), float(ci[1])],
                shapiro_p=float(sw.pvalue)), T, Xz, keep


ap = argparse.ArgumentParser()
ap.add_argument('--data', default='descriptor_cache_xtb.csv')
ap.add_argument('--simulant-classes', default='simulant',
                help="최근접 '모사체' 로 볼 compound_class (쉼표 구분). "
                     "확장 패널에서 literature_simulant 를 포함하려면 지정.")
ap.add_argument('--suffix', default='', help='산출 파일 접미사')
ap.add_argument('--drop-pathological', action='store_true',
                help='소분자에서 음수가 되는 kappa2/kappa3 등 정의붕괴 서술자 제외')
A = ap.parse_args()
SIM_CLASSES = tuple(s.strip() for s in A.simulant_classes.split(','))

df = pd.read_csv(A.data)
N_PANEL = len(df)
assert N_PANEL >= 28, N_PANEL
# 2026-10-12: 패널이 28종으로 고정돼 있었다. 확장 패널(35종)을 같은 코드로
# 돌리기 위해 크기를 데이터에서 읽는다.
print(f'[load] {len(df)}종  ' +
      '  '.join(f'{k}={v}' for k, v in df.compound_class.value_counts().items()))
print('[load] 모사체 정본:', ', '.join(df[df.compound_class == 'simulant'].compound_name))

num = [c for c in df.select_dtypes('number').columns
       if df[c].notna().sum() == N_PANEL and c not in META_NUM]
patho = [c for c in num if (df[c] < 0).any() and c.startswith('kappa')]
print(f'[filter] {N_PANEL}종 전부 보유 수치열 {len(num)}개')
if patho:
    print(f'[filter] 소분자에서 음수가 되는 형상지수: {patho} '
          f'(min {", ".join(f"{c}={df[c].min():.1f}" for c in patho)})')
if A.drop_pathological and patho:
    num = [c for c in num if c not in patho]
    print(f'[filter] --drop-pathological 적용 -> {len(num)}개')

rdkit_cols = [c for c in num if c not in ELECTRONIC_28]
elec_cols = [c for c in num if c in ELECTRONIC_DIST]
# 병합공간에서도 omega 를 뺀다. 거리에 쓰이는 모든 공간에서 제외한다.
merged_cols = [c for c in num if c != OMEGA_DERIVED]
print(f'[filter] RDKit {len(rdkit_cols)} / 전자 {len(elec_cols)} / '
      f'병합 {len(merged_cols)}')

res, tabs = [], []
for tag, cols in [('A_rdkit', rdkit_cols), ('B_electronic', elec_cols),
                  ('C_merged', merged_cols)]:
    r, T, _, _ = analyse(df, cols, tag)
    res.append(r); tabs.append(T)

TT = pd.concat(tabs, ignore_index=True)
suffix = (A.suffix or '') + ('_nopatho' if A.drop_pathological else '')
TT.to_csv(f'outputs/distance_table{suffix}.csv', index=False)

# ---------------------------------------------------- 인 중심 국소 반응성 (6종)
P = df[df.q_P.notna()][['compound_name', 'compound_class', 'agent_series'] + P_LOCAL]
P = P.sort_values('fukui_f_plus_P', ascending=False)
print('\n=== [P-local] 인 함유 6종 (f+(P) 내림차순) ===')
print(P.to_string(index=False, float_format=lambda v: f'{v:.5f}'))
g = P[P.agent_series == 'G'].fukui_f_plus_P
dm = float(P.loc[P.compound_name == 'DMMP', 'fukui_f_plus_P'].iloc[0])
gb = float(P.loc[P.compound_name == 'sarin (GB)', 'fukui_f_plus_P'].iloc[0])
vx = float(P.loc[P.compound_name == 'VX', 'fukui_f_plus_P'].iloc[0])
gG = P[P.agent_series == 'G']
print(f'\n  G계열(GB,GD,GA,GF) f+(P) = {g.mean():.5f} +/- {g.std(ddof=1):.5f}')
gn = P[P.compound_name.isin(['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)'])].fukui_f_plus_P
print(f'  GB/GD/GF만    f+(P) = {gn.mean():.5f} +/- {gn.std(ddof=1):.5f}  '
      f'-> DMMP z = {(dm - gn.mean())/gn.std(ddof=1):+.2f}')
print(f'  DMMP/GB = {dm/gb:.3f}   VX/GB = {vx/gb:.3f}   DMMP/VX = {dm/vx:.3f}')
print('  wiberg_PO:  ' + ', '.join(
    f'{r.compound_name}={r.wiberg_PO:.4f}' for r in P.itertuples()))

out = dict(seed=SEED, bootstrap_samples=N_BOOT, corr_cutoff=R_CUT,
           drop_pathological=bool(A.drop_pathological),
           pathological_detected=patho,
           class_counts={k: int(v) for k, v in df.compound_class.value_counts().items()},
           simulants=list(df[df.compound_class == 'simulant'].compound_name),
           conditions=res,
           p_local=dict(
               table={r.compound_name: dict(q_P=r.q_P, q_phosphoryl_O=r.q_phosphoryl_O,
                                            wiberg_PO=r.wiberg_PO,
                                            fukui_f_plus_P=r.fukui_f_plus_P)
                      for r in P.itertuples()},
               G3_mean=float(gn.mean()), G3_sd=float(gn.std(ddof=1)),
               DMMP_z_vs_G3=float((dm - gn.mean()) / gn.std(ddof=1)),
               DMMP_over_GB=float(dm / gb), VX_over_GB=float(vx / gb),
               DMMP_over_VX=float(dm / vx)))
with open(f'outputs/experiment_summary{suffix}.json', 'w') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print(f'\n[done] outputs/experiment_summary{suffix}.json, distance_table{suffix}.csv')
