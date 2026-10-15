#!/usr/bin/env python3
"""Step 4 마감 산출물: raw_results.csv (long) + results_table.tex"""
import json
import numpy as np
import pandas as pd

S = json.load(open('outputs/experiment_summary.json'))
SAN = json.load(open('outputs/sanity_checks.json'))
H1 = json.load(open('outputs/h1_summary.json'))
QS = json.load(open('outputs/qsar_tox_summary.json'))
NORM = pd.read_csv('outputs/distance_table_normalized.csv')
TOX = pd.read_csv('outputs/toxicity_aegl.csv')
CMP = pd.read_csv('outputs/xtb_fast_vs_slow.csv')
C = pd.read_csv('descriptor_cache_xtb.csv')

rows = []


def add(block, cond, item, metric, value, unit='', note=''):
    rows.append(dict(block=block, condition=cond, item=item, metric=metric,
                     value=value, unit=unit, note=note))


# --- A. 서술자 공간 ablation
for c in S['conditions']:
    k = c['condition']
    add('descriptor_space', k, '-', 'n_descriptors_kept', c['n_kept'])
    add('descriptor_space', k, '-', 'PC1_variance_explained', round(c['evr'][0], 4))
    add('descriptor_space', k, '-', 'PC1_PC2_cumulative', round(c['pc1_pc2'], 4))
    for n, v in c['pc1_loadings'].items():
        add('descriptor_space', k, n, 'PC1_loading', round(v, 4))
    add('descriptor_space', k, '-', 'n_interferent_closer', c['n_interferent_closer'],
        'of 11')

# --- B. 최근접 이웃 거리 (정규화 포함)
for r in NORM.itertuples():
    add('nearest_neighbour', r.condition, r.real_agent, 'd_nearest_simulant',
        round(r.d_simulant, 4), 'z-dist', r.nearest_simulant)
    add('nearest_neighbour', r.condition, r.real_agent, 'd_nearest_interferent',
        round(r.d_interferent, 4), 'z-dist', r.nearest_interferent)
    add('nearest_neighbour', r.condition, r.real_agent, 'd_simulant_normalised',
        round(r.d_sim_norm, 4), 'x typical NN',
        'valid match' if (r.is_match and not r.interferent_closer) else '')

# --- C. Sanity
m = SAN['sanity_monotonicity']
for k, v in zip(m['k_grid'], m['mean_distance']):
    add('sanity_monotonicity', 'C_merged', f'k={k}', 'mean_pairwise_distance',
        round(v, 4))
add('sanity_monotonicity', 'C_merged', '-', 'passed', int(m['passed']))
for k, v in SAN['sanity_baseline'].items():
    add('sanity_baseline', k, '-', 'null_mean_count', round(v['null_mean_count'], 3),
        'of 11', 'C(17,3)=680 전수')
    add('sanity_baseline', k, '-', 'observed_count', v['observed_count'], 'of 11')
    add('sanity_baseline', k, '-', 'observed_mean_diff', round(v['observed_diff'], 4))
    add('sanity_baseline', k, '-', 'exact_p', round(v['p_exact'], 4))

# --- D. 프로토콜 일관성 (빠른 vs 정밀 xtb)
for col in ['E_HOMO_eV', 'E_LUMO_eV', 'electrophilicity_index_eV', 'dipole_moment_D']:
    s = CMP[[col, col + '_slow']].dropna()
    add('xtb_protocol', 'fast_vs_precise', col, 'pearson_r',
        round(float(np.corrcoef(s[col], s[col + '_slow'])[0, 1]), 4), '', f'n={len(s)}')
    add('xtb_protocol', 'fast_vs_precise', col, 'mean_abs_diff',
        round(float(np.abs(s[col] - s[col + '_slow']).mean()), 4))

# --- E. 인 중심 반응성
for n, v in S['p_local']['table'].items():
    for k2, v2 in v.items():
        add('p_local', '-', n, k2, round(v2, 5))

# --- F. 독성 실측
for r in TOX.itertuples():
    add('toxicity_measured', '-', r.compound_name, 'AEGL3_10min',
        (round(r.AEGL3_10min_mgm3, 4) if r.AEGL3_10min_mgm3 == r.AEGL3_10min_mgm3
         else ''), 'mg/m3', r.note if isinstance(r.note, str) else '')

# --- G. QSAR
for c in QS['cv']:
    add('qsar_tox', c['tag'], '-', 'R2', round(c['r2_mean'], 4), '',
        f"sd {c['r2_sd']:.4f}, seeds {c['n_seeds']}")
    add('qsar_tox', c['tag'], '-', 'RMSE_log10', round(c['rmse_mean'], 4), '',
        f"95% CI [{c['rmse_ci95'][0]:.3f}, {c['rmse_ci95'][1]:.3f}]")

# --- H. H1
for k in ['spearman_r', 'spearman_p', 'pearson_r', 'pearson_p', 'r_squared',
          'n_candidates']:
    add('h1_decoupling', '-', '-', k, round(H1[k], 5) if isinstance(H1[k], float)
        else H1[k], '',
        f"95% CI [{H1['spearman_ci'][0]:.3f}, {H1['spearman_ci'][1]:.3f}]"
        if k == 'spearman_r' else '')
for k, v in H1['quadrants'].items():
    add('rq2_quadrant', '-', k, 'n_compounds', v)

R = pd.DataFrame(rows)
R.to_csv('outputs/raw_results.csv', index=False)
print(f'[raw] outputs/raw_results.csv  {len(R)}행')

# ================================================================ LaTeX
N = {c['condition']: c for c in S['conditions']}
B = SAN['sanity_baseline']
tex = []
tex.append(r"""% simulant / Step 4 — 결과 표 (Step 5 본문에 직접 삽입)
\begin{table}[t]
\centering
\caption{Descriptor-space ablation. Nearest-neighbour comparison for the eleven
authentic agents, with the exact permutation test over all $\binom{17}{3}=680$ ways
of designating three of the seventeen non-agents as simulants.}
\label{tab:ablation}
\begin{tabular}{lccccc}
\hline
Descriptor space & Kept & PC1 identity & PC1+PC2 & Interferent closer & Exact $p$ \\
\hline""")
LBL = {'A_rdkit': 'Structural / physicochemical', 'B_electronic': 'GFN2 electronic',
       'C_merged': 'Merged'}
PC1 = {'A_rdkit': 'molecular size', 'B_electronic': 'global electrophilicity',
       'C_merged': 'molecular size'}
for k in ['A_rdkit', 'B_electronic', 'C_merged']:
    c = N[k]
    tex.append(f"{LBL[k]} & {c['n_kept']} & {PC1[k]} & {c['pc1_pc2']*100:.1f}\\% & "
               f"{c['n_interferent_closer']} / 11 & {B[k]['p_exact']:.3f} \\\\")
tex.append(r"""\hline
\multicolumn{6}{l}{\footnotesize Null expectation for ``interferent closer'' is
%.2f / 11, set by the 14:3 group-size asymmetry alone.}\\
\end{tabular}
\end{table}
""" % B['A_rdkit']['null_mean_count'])

M = NORM[NORM.condition == 'C_merged'].sort_values('d_sim_norm')
tex.append(r"""\begin{table}[t]
\centering
\caption{Nearest designated simulant for each authentic agent, expressed as a
multiple of the typical nearest-neighbour distance of the same descriptor space
(merged space shown). A match requires both a normalised distance $\le 2$ and a
distance smaller than that to every interferent.}
\label{tab:match}
\begin{tabular}{llcc}
\hline
Authentic agent & Nearest simulant & Normalised distance & Valid match \\
\hline""")
SH = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
      'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD', 'lewisite 1 (L)': 'L',
      'hydrogen cyanide (AC)': 'AC', 'cyanogen chloride (CK)': 'CK',
      'phosgene (CG)': 'CG', 'chloropicrin (PS)': 'PS', 'VX': 'VX'}
for r in M.itertuples():
    ok = r.is_match and not r.interferent_closer
    tex.append(f"{SH.get(r.real_agent, r.real_agent)} & "
               f"{r.nearest_simulant.replace('dichloromethane','DCM')} & "
               f"{r.d_sim_norm:.2f} & {'yes' if ok else 'no'} \\\\")
tex.append(r"""\hline
\end{tabular}
\end{table}
""")

P = S['p_local']['table']
tex.append(r"""\begin{table}[t]
\centering
\caption{GFN2-xTB local reactivity at phosphorus for the six phosphorus-bearing
compounds. $\omega$ is a global index and is reported for completeness only.}
\label{tab:plocal}
\begin{tabular}{lccccc}
\hline
Compound & $q_{\mathrm{P}}$ (e) & $q_{\mathrm{P=O}}$ (e) & $W_{\mathrm{P=O}}$ &
$f^{+}(\mathrm{P})$ & $\omega$ (eV) \\
\hline""")
OM = dict(zip(C.compound_name, C.electrophilicity_index_eV))
for n in ['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)', 'DMMP', 'tabun (GA)', 'VX']:
    v = P[n]
    nm = SH.get(n, n)
    tex.append(f"{'\\textbf{' + nm + '}' if n == 'DMMP' else nm} & "
               f"{v['q_P']:.5f} & {v['q_phosphoryl_O']:.5f} & "
               f"{v['wiberg_PO']:.4f} & {v['fukui_f_plus_P']:.5f} & {OM[n]:.2f} \\\\")
g3 = [P[x]['fukui_f_plus_P'] for x in ['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)']]
tex.append(r"""\hline
\multicolumn{6}{l}{\footnotesize GB/GD/GF mean $f^{+}(\mathrm{P})$ = %.4f $\pm$ %.4f;
DMMP lies %.2f SD below (%.0f\%% of sarin) and %.2f$\times$ above VX.}\\
\end{tabular}
\end{table}
""" % (np.mean(g3), np.std(g3, ddof=1), abs(S['p_local']['DMMP_z_vs_G3']),
       100 * S['p_local']['DMMP_over_GB'], S['p_local']['DMMP_over_VX']))

tex.append(r"""\begin{table}[t]
\centering
\caption{Decoupling test. Fidelity is the standardised electronic-space distance to
the nearest authentic agent; safety is the measured EPA AEGL-3 (10~min).}
\label{tab:h1}
\begin{tabular}{lc}
\hline
Quantity & Value \\
\hline""")
_D = pd.read_csv('outputs/h1_candidates.csv')
_cl = _D[~_D.alert_flag]
_fid, _bar = H1['fidelity_cut'], H1['tox_cut']
_n_str = int(((_cl.fidelity_dist <= _fid) & (_cl.y_tox >= _bar)).sum())
_rel = _cl[(_cl.fidelity_dist <= _fid) & (_cl.y_tox >= _bar - 1.0)]
_dcmf = float(_D.loc[_D.label == 'dichloromethane', 'fidelity_dist'].iloc[0])
_top = _rel.nsmallest(1, 'fidelity_dist')
tex.append(f"Candidates with both axes measured & {H1['n_candidates']} \\\\")
tex.append(f"Spearman $r$ & ${H1['spearman_r']:+.3f}$ "
           f"(95\\% CI ${H1['spearman_ci'][0]:+.3f}$ to ${H1['spearman_ci'][1]:+.3f}$) \\\\")
tex.append(f"Pearson $r$ & ${H1['pearson_r']:+.3f}$ "
           f"(95\\% CI ${H1['pearson_ci'][0]:+.3f}$ to ${H1['pearson_ci'][1]:+.3f}$), "
           f"$p = {H1['pearson_p']:.2f}$ \\\\")
tex.append(f"Variance of toxicity explained by fidelity & "
           f"{100*H1['r_squared']:.1f}\\% \\\\")
tex.append(r"\hline")
tex.append(f"Safety bar: dichloromethane, {H1['dcm_aegl3_mgm3']:,.0f}~mg~m$^{{-3}}$ "
           f"& $\\log_{{10}} = {_bar:.2f}$ \\\\")
tex.append(f"Faithful and safe at that bar & {_n_str} \\\\")
tex.append(f"Faithful and safe one order of magnitude below it & {len(_rel)} \\\\")
tex.append(f"\\quad clearing the structural-alert gate & "
           f"{int((~_rel.alert_flag).sum())} \\\\")
tex.append(f"\\quad more faithful than dichloromethane itself & "
           f"{int((_rel.fidelity_dist < _dcmf).sum())} \\\\")
tex.append(f"Most faithful of them & {_top.label.iloc[0]} "
           f"({_top.fidelity_dist.iloc[0]:.2f} vs {_dcmf:.2f}) \\\\")
tex.append(r"""\hline
\end{tabular}
\end{table}
""")
open('outputs/results_table.tex', 'w').write('\n'.join(tex))
print('[tex] outputs/results_table.tex  표 4종')
