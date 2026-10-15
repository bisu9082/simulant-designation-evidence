#!/usr/bin/env python3
"""본문 수치를 산출물과 대조한다 (2026-10-11 전면 재작성, 방향을 뒤집었다).

종전 설계의 결함. 핵심 수치 목록을 손으로 적어 두고 그 문자열이 main.tex 에
'존재하는지'만 확인했다. 존재 확인은 모순을 잡지 못한다. 그래서
  - 본문이 같은 양을 208 과 218 로 동시에 적어도 통과했고,
  - 같은 환산값을 41,683 과 41,685 로 적어도 통과했으며,
  - 19,531 과 19,530 도 통과했다.
GATE 8 3차에서 전부 지적됐다.

재설계. 산출물(JSON/CSV)에서 사실표를 만든 뒤, 본문에서 그 양을 가리키는
문맥을 모두 찾아 값을 추출하고 사실표와 비교한다. 두 가지를 본다.
  (A) 문맥 대조: 한 양을 가리키는 모든 등장이 산출물 값과 같은가.
      같은 양을 두 값으로 적으면 반드시 걸린다.
  (B) 근사 불일치: 본문의 수치 중 사실표의 어떤 값과 0.5% 안쪽으로
      가까우면서 같지는 않은 것. 반올림 분기(84.93 대 84.933)를 잡는다.
"""
import json
import re
import sys

import pandas as pd

TEX = 'paper/main.tex'


def load_facts():
    """산출물에서만 읽는다. 여기에 손으로 적은 수치는 없다."""
    f = {}
    q = json.load(open('outputs/qsar_tox_summary.json'))
    f['n_train'] = q['n_train']
    f['n_heldout'] = q['n_heldout_cwa']
    f['n_features'] = q['n_features']
    f['ad_threshold'] = q['ad_threshold']
    for row in q['cv']:
        if row['tag'].startswith('RandomForest'):
            f['rf_r2'] = row['r2_mean']
            f['rf_r2_sd'] = row['r2_sd']
            f['rf_rmse'] = row['rmse_mean']
            f['rf_rmse_lo'], f['rf_rmse_hi'] = row['rmse_ci95']
        if row['tag'].startswith('baseline-1'):
            f['b1_rmse'] = row['rmse_mean']
            f['b1_lo'], f['b1_hi'] = row['rmse_ci95']
        if row['tag'].startswith('baseline-2'):
            f['b2_rmse'] = row['rmse_mean']

    h = json.load(open('outputs/h1_summary.json'))
    f['n_candidates'] = h['n_candidates']
    f['spearman_r'] = h['spearman_r']
    f['spearman_lo'], f['spearman_hi'] = h['spearman_ci']
    f['pearson_r'] = h['pearson_r']
    f['pearson_lo'], f['pearson_hi'] = h['pearson_ci']
    f['dcm_aegl3'] = h['dcm_aegl3_mgm3']
    t = h['safety_tiers']
    f['n_relaxed_10'] = t['relaxed_10']['n']
    f['n_relaxed_100'] = t['relaxed_100']['n']
    f['n_incumbent'] = t['incumbent']['n']
    f['relaxed_bar'] = t['relaxed_10']['bar_mgm3']

    m = json.load(open('outputs/minority_regression.json'))
    f['n_ache'] = m['n_ache']
    f['n_ache_op'] = m['n_by_group']['organophosphorus_ester']
    f['n_ache_carb'] = m['n_by_group']['carbamate_ester']
    f['low10_err'] = m['low10_err']
    f['low10_n'] = m['low10_n']
    f['high10_err'] = m['high10_err']
    f['phorate_fold'] = m['most_toxic_ache']['fold_too_safe']
    f['malathion_fold'] = m['least_toxic_ache']['fold_too_toxic']

    v = json.load(open('outputs/aegl_parse_validation.json'))
    f['source_with_aegl3'] = v['source_with_aegl3']
    f['trainset_checked'] = v['trainset_checked']
    f['panel_checked'] = v['panel_checked']
    f['comma_rows'] = v['comma_rows_source']

    T = pd.read_csv('outputs/toxicity_aegl.csv')
    f['panel_with_aegl3'] = int(T.AEGL3_10min_mgm3.notna().sum())
    f['panel_from_footnote'] = int((T.AEGL3_provenance == 'footnote').sum())
    _H = T[T.AEGL3_10min_mgm3.notna()]
    _dv = float(_H[_H.compound_name == 'dichloromethane'].AEGL3_10min_mgm3.iloc[0])
    f['dcm_panel_rank'] = int((_H.AEGL3_10min_mgm3 > _dv).sum()) + 1
    f['_dcm_is_safest'] = f['dcm_panel_rank'] == 1
    f['_panel_safer_than_dcm'] = sorted(
        _H[_H.AEGL3_10min_mgm3 > _dv].compound_name.tolist())
    f['panel_without'] = int(T.AEGL3_10min_mgm3.isna().sum())

    E = pd.read_csv('outputs/qsar_heldout.csv')
    f['heldout_n'] = int(len(E))
    f['heldout_in_domain'] = int(E.in_domain.sum())
    f['heldout_worst_in_domain'] = float(E[E.in_domain].error_log10.max())
    f['heldout_err_min'] = float(E.error_log10.min())
    f['heldout_err_max'] = float(E.error_log10.max())

    C = pd.read_csv('outputs/h1_candidates.csv')
    _d = C[C.label == 'dichloromethane']
    _t = C[C.label == 'Trichloroethylene']
    if len(_d) and len(_t):
        f['tce_ratio'] = float(_d.fidelity_dist.iloc[0] / _t.fidelity_dist.iloc[0])
        f['dcm_above'] = int((C.y_tox > float(_d.y_tox.iloc[0])).sum())
    f['panel_without'] = int(T.AEGL3_10min_mgm3.isna().sum())
    for nm, key in [('dichloromethane', 'dcm'), ('methanol', 'methanol'),
                    ('n-hexane', 'hexane'), ('chloroform', 'chloroform'),
                    ('toluene', 'toluene'), ('acetone', 'acetone'),
                    ('benzene', 'benzene'), ('p-xylene', 'pxylene'),
                    ('cyclosarin (GF)', 'gf'), ('VX', 'vx')]:
        r = T[T.compound_name == nm]
        if len(r) and r.AEGL3_10min_mgm3.notna().all():
            f[f'aegl3_{key}'] = float(r.AEGL3_10min_mgm3.iloc[0])

    P = pd.read_csv('outputs/toxicity_predicted.csv')
    d = P[P.compound_name == 'DMMP']
    if len(d):
        f['dmmp_ad'] = float(d.ad_dist.iloc[0])

    s = json.load(open('outputs/permutation_sidedness.json'))
    f['p_struct'] = s['A_rdkit']['one_sided']
    f['p_elec'] = s['B_electronic']['one_sided']
    f['p_merged'] = s['C_merged']['one_sided']

    # 롤백 C 의 헤드라인 수치. 2026-10-13 신설. 이 값들이 본문의 주장을 떠받치므로
    # 전량 등재한다. 종전에는 헤드라인이 검증기에 없어서 전파 실패가 통과했다.
    rb = json.load(open('outputs/rollback_c.json'))
    dc = rb['dimensionality_control']
    f['pc3_p28'] = dc['panel28']['pc3_p']
    f['omega_dev'] = rb['omega_dependence']['panel28'][
        'max_deviation_from_definition']
    f['r_lumo_omega'] = rb['omega_dependence']['panel28']['r_lumo_omega']
    f['omega_distrank'] = rb['omega_dependence']['panel28'][
        'dist_rank_correlation_4d_vs_3d']
    f['omega_pc3_var'] = 100.0 * rb['omega_dependence']['panel28'][
        'pc3_cumulative_variance']
    for m in rb['deltaE']['mantel']:
        if m['dataset'] == 'd2' and m['mode'] == 'uv' and \
                m['space'] == 'structural':
            f['mantel_d2uv_struct'] = m['p_mantel']
            f['mantel_d2uv_struct_rho'] = m['spearman']
            f['naive_d2uv_struct'] = m['p_pairs_as_independent']
        if m['dataset'] == 'd2' and m['mode'] == 'normal' and \
                m['space'] == 'electronic3':
            f['mantel_d2nor_elec'] = m['p_mantel']
            f['mantel_d2nor_elec_rho'] = m['spearman']
            f['naive_d2nor_elec'] = m['p_pairs_as_independent']
    for fl in rb['deltaE']['centring_floor']:
        if fl['dataset'] == 'd1' and fl['mode'] == 'uv':
            f['floor_d1uv_obs'] = fl['observed_mean_cosine']
            f['floor_d1uv_floor'] = fl['algebraic_floor']
    for rk in rb['dmmp_ranking']:
        if rk['dataset'] == 'd1' and rk['mode'] == 'uv':
            for e in rk['ranking']:
                if e['agent'] == 'GA':
                    f['dmmp_ga_uv'] = e['cosine']
                if e['agent'] == 'GB':
                    f['dmmp_gb_uv'] = e['cosine']

    f['p_struct_fig'] = f['p_struct']
    f['p_elec_fig'] = f['p_elec']
    f['p_merged_fig'] = f['p_merged']

    f['mantel_minq'] = rb['deltaE']['mantel_bh']['min_q']

    ss = json.load(open('outputs/stat_supplement.json'))
    w = ss['whitening']['panel28']
    f['rand3_frac28'] = 100.0 * w['raw']['frac_below_05']
    f['rand3_med28'] = w['raw']['quantiles']['q50']
    f['white_frac28'] = 100.0 * w['whitened']['frac_below_05']
    f['white_med28'] = w['whitened']['quantiles']['q50']
    s2 = json.load(open('outputs/stat_supplement2.json'))
    wy = s2['westfall_young']
    f['wy35'] = wy['panel35_incumbent3']['wy_adjusted_p']
    f['wy35lit'] = wy['panel35_incumbent_plus_lit8']['wy_adjusted_p']
    f['wy28'] = wy['panel28_incumbent3']['wy_adjusted_p']
    bt = s2['distance_bootstrap']
    f['boot_inc_med'] = bt['incumbent3']['median']
    f['boot_inc_lo'] = bt['incumbent3']['ci'][0]
    f['boot_inc_hi'] = bt['incumbent3']['ci'][1]
    f['boot_inc_zero'] = 100.0 * bt['incumbent3']['frac_zero']
    f['boot_lit_med'] = bt['incumbent_plus_lit8']['median']
    f['boot_diff_med'] = bt['difference']['median']
    f['boot_diff_pos'] = round(bt['difference']['frac_positive'] * bt['n_boot'])
    _r3 = bt['incumbent3']['retention']
    _r8 = bt['incumbent_plus_lit8']['retention']
    f['ret_gb3'] = 100.0 * _r3.get('sarin (GB)', 0)
    f['ret_gd3'] = 100.0 * _r3.get('soman (GD)', 0)
    f['ret_gf3'] = 100.0 * _r3.get('cyclosarin (GF)', 0)
    f['ret_ga8'] = 100.0 * _r8.get('tabun (GA)', 0)

    tf = ss['threshold_free']
    f['unbounded35'] = tf['panel35_incumbent3']['n_match_unbounded']
    f['null_exp35'] = tf['panel35_incumbent3']['null_expectation_outside_closer']
    f['null_exp28'] = tf['panel28_incumbent3']['null_expectation_outside_closer']

    fam = ss['family']['panel28_incumbent3']
    f['q_by_space28'] = fam['bh_by_space_min_q']
    f['q_prespec28'] = fam['prespecified_bh_min_q']
    f['meff28'] = fam['meff']
    f['meff_p28'] = fam['meff_adjusted_min_p']

    mb = json.load(open('outputs/match_by_panel.json'))
    f['match35_inc'] = mb['panel35_incumbent3']['n_match']
    f['match35_lit'] = mb['panel35_incumbent_plus_lit8']['n_match']
    f['match28_inc'] = mb['panel28_incumbent3']['n_match']
    _lit = {a['agent']: a for a in mb['panel35_incumbent_plus_lit8']['agents']}
    _inc = {a['agent']: a for a in mb['panel35_incumbent3']['agents']}
    f['dfp_gb'] = _lit['sarin (GB)']['d_norm']
    f['dfp_gd'] = _lit['soman (GD)']['d_norm']
    f['dfp_gf'] = _lit['cyclosarin (GF)']['d_norm']
    f['dmmp_gb_35'] = _inc['sarin (GB)']['d_norm']
    f['cees_hd'] = _lit['sulfur mustard (HD)']['d_norm']

    pc = json.load(open('outputs/dim_percentiles.json'))
    f['pct_struct'] = pc['struct_full']
    f['pct_elec'] = pc['electronic']
    f['pct_pc3'] = pc['struct_pc3']

    mp = json.load(open('outputs/multiplicity_c.json'))
    for lab, _m, _minp, _bonf, bhq, _k in mp:
        if lab == 'panel28_designated3':
            f['bh_q28'] = bhq
        if lab == 'panel35_designated3':
            f['bh_q35'] = bhq
        if lab == 'panel35_designated_plus_literature8':
            f['bh_q35lit'] = bhq
    return f


def nums_in(s):
    return [float(x.replace(',', '').replace('\\,', ''))
            for x in re.findall(r'(?<![\d.])(\d[\d,]*(?:\\,\d{3})*(?:\.\d+)?)', s)]


# (사실표 키, 본문에서 그 양을 가리키는 문맥 정규식, 허용 오차)
# 문맥 정규식의 첫 캡처 그룹이 검사할 수치다.
CHECKS = [
    ('n_train', [r'(\d[\d,]*) training (?:chemicals|compounds)',
                 r'leaving (\d[\d,]*) for\s*\n?training',
                 r'trained on (\d[\d,]*) AEGL'], 0),
    ('n_heldout', [r'(?:the|all) (twenty|seventeen|nineteen|eighteen) '
                   r'(?:military )?agents held out',
                   r'On the (twenty|seventeen) agents held out'], 0),
    ('n_candidates', [r'(?:the|Across the|across the) (\d[\d,]*) candidates',
                      r'excluding the authentic agents\s*\n?themselves \(\$n=(\d+)\$\)',
                      r'(\d[\d,]*)-compound\s*\n?AEGL pool',
                      r'all but \w+ of the (\d[\d,]*) compounds in the',
                      r'independent across the (\d[\d,]*) compounds',
                      r'for the (\d[\d,]*) candidates\s*\n?carrying both axes',
                      r'of the (\d[\d,]*), none combines'], 0),
    ('dcm_aegl3', [r'dichloromethane[,\s]+(\d[\d,\\]*)~mg~m',
                   r'at\s*\n?(\d[\d,\\]*)~mg~m\$\^\{-3\}\$ it sits at the',
                   r"dichloromethane's (\d[\d,\\]*)~mg~m"], 1.0),
    ('aegl3_chloroform', [r'chloroform \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_toluene', [r'toluene \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_methanol', [r'methanol \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_hexane', [r'n-hexane\s*\n?\((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_acetone', [r'acetone \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_benzene', [r'benzene \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_pxylene', [r'p-xylene \((\d[\d,\\]*)\)'], 1.0),
    ('rf_r2', [r'R\^\{2\}=([\d.]+) \\pm'], 0.001),
    ('rf_r2_sd', [r'R\^\{2\}=[\d.]+ \\pm ([\d.]+)'], 0.001),
    ('rf_rmse', [r'RMSE \$?([\d.]+)\$? in'], 0.001),
    ('b1_rmse', [r'training-mean baseline of\s*\n?([\d.]+) \('], 0.001),
    ('b2_rmse', [r'linear baseline of ([\d.]+)\.'], 0.001),
    ('dmmp_ad', [r"DMMP's domain distance, ([\d.]+),"], 0.005),
    ('n_ache', [r'(Thirteen|Eleven|Twelve) of the\s*\n?\d+ training compounds carry'], 0),
    ('n_relaxed_10', [r'(Eight|Seven|Nine) compounds then combine',
                      r'there are (eight|seven|nine) candidates'], 0),
    ('n_relaxed_100', [r'the count rises to (fifteen|seventeen|sixteen)'], 0),
    ('panel_with_aegl3', [r'(Nineteen|Thirteen|Fourteen|Eighteen) of the '
                          r'twenty-eight compounds have an AEGL-3',
                          r'the (nineteen|thirteen|fourteen) panel compounds',
                          r'(nineteen|thirteen|fourteen) compounds that have one'], 0),
    ('source_with_aegl3', [r'agree on all (\d[\d,]*) entries',
                           r'which leaves (\d[\d,]*) \(250 from a table cell'], 0),
    ('trainset_checked', [r'on all (\d[\d,]*) of\s*\n?those that resolved to a structure'], 0),
    ('comma_rows', [r'covers all (\d[\d,]*) values at or above'], 0),
    ('panel_from_footnote',
     [r'(Five|Nine|Four) of the nineteen values come from the',
      r'(Five|Nine|Four) of the nineteen panel values'], 0),
    ('panel_without', [r'\(b\) The (nine|fifteen|ten) compounds with no regulator-grade',
                       r'(Nine|Fifteen|Ten) compounds have no value'], 0),
    ('heldout_worst_in_domain',
     [r'worst of those fourteen being\s*\n?wrong by ([\d.]+) log units',
      r'is still wrong\s*\n?by ([\d.]+) log units'], 0.005),
    ('heldout_in_domain',
     [r'(Fourteen|Thirteen|Fifteen) of the twenty lie inside',
      r'(fourteen|thirteen|fifteen) of them pass the domain check'], 0),
    ('tce_ratio', [r'trichloroethylene[, ]+is ([\d.]+) times closer',
                   r'trichloroethylene ([\d.]+) times closer'], 0.05),
    ('dcm_above', [r'with (seven|eight|two) compounds\s*\n?above it',
                   r'all but (seven|eight|two) of the 214 compounds'], 0),
    ('panel_without', [r'inhalation value exists for any of these (\d+)'], 0),
    ('dcm_panel_rank', [r'designated simulant, ranks (\d+) of \d+'], 0),
    ('panel_with_aegl3', [r'designated simulant, ranks \d+ of (\d+)'], 0),
    ('heldout_in_domain', [r'(\d+) of the 20 agents pass the domain check'], 0),
    ('heldout_n', [r'and all (\d+) are wrong by'], 0),
    ('heldout_worst_in_domain',
     [r'worst inside the domain is wrong by ([\d.]+) log'], 0.005),
    ('spearman_r', [r'Spearman \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?',
                    r'Spearman \$r=([-+][\d.]+)\$'], 0.0006),
    ('spearman_lo', [r'Spearman[^\n]*?95\\% CI \$?([-+][\d.]+)\$? to'], 0.001),
    ('spearman_hi', [r'Spearman[^\n]*?95\\% CI [^\n]*?to\s*\n?\$?([-+][\d.]+)\$?'], 0.001),
    ('pearson_r', [r'Pearson \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?'], 0.0006),
    # 2026-10-13: 본문이 재집필되어 종전 정규식의 문맥이 사라졌다. kappa 민감도
    # 문장과 Table 1 의 '28종 패널' 블록을 동시에 대조한다. 표에는 패널 3구성이
    # 들어 있으므로 '28 cpd.' 행에 앵커를 걸어야 한다. 앵커 없이 쓰면 35종 행의
    # 값까지 물어와 거짓 불일치가 난다.
    ('p_struct', [r'condition \(structural \$p=([\d.]+)\$',
                  r'28 cpd\., incumbent 3 & structural\s+& 38 & \\textbf\{([\d.]+)\}'],
     0.0005),
    ('p_elec', [r'electronic \$p=([\d.]+)\$',
                r'28 cpd\.[^\n]*\n[^\n]*& electronic\s+&\s+3 & ([\d.]+) &'], 0.0005),
    ('p_merged', [r'merged \$p=([\d.]+)\$',
                  r'28 cpd\.[^\n]*\n[^\n]*\n[^\n]*& merged\s+& 41 & \\textbf\{([\d.]+)\}'],
     0.0005),
    ('pc3_p28', [r'from \$p=0.028\$ to\s*\n?\$p=([\d.]+)\$'], 0.0005),
    ('rand3_frac28', [r'and ([\d.]+)\\% of the subspaces reach'], 0.05),
    ('rand3_med28', [r'panel the median \$p\$ is ([\d.]+)'], 0.0005),
    ('r_lumo_omega', [r'r\(E_\{\\mathrm\{LUMO\}\}, \\omega\) =\s*\n?([-+][\d.]+)'],
     0.0005),
    ('omega_distrank', [r'\(Spearman ([\d.]+)\)'], 0.0005),
    ('omega_pc3_var', [r'account for ([\d.]+)\\% of its variance'], 0.05),
    ('bh_q28', [r'smallest adjusted value of ([\d.]+)',
                r'Benjamini-Hochberg \$q =\s*\n?([\d.]+)\$'], 0.0005),
    ('bh_q35', [r'Benjamini-Hochberg over 18\s+& ([\d.]+) &'], 0.0005),
    ('bh_q35lit', [r'Benjamini-Hochberg over 18\s+& [\d.]+ & ([\d.]+) &'], 0.0005),
    ('mantel_d2uv_struct', [r'\\rho = \+0.539\$, Mantel \$p = ([\d.]+)\$'], 0.0005),
    ('mantel_d2nor_elec', [r'\$\\rho = \+0\.571\$,\s*\n?Mantel \$p = ([\d.]+)\$'], 0.0005),
    ('mantel_d2uv_struct_rho', [r'ultraviolet channel \(\$\\rho = ([-+][\d.]+)\$'], 0.0015),
    ('floor_d1uv_obs', [r'The observed means clear that floor by 0.013 to\s*\n?0.066: \$([-+][\d.]+)\$'], 0.0015),
    ('dmmp_ga_uv', [r'tabun, at \$\+?([\d.]+)\$ against sarin'], 0.0015),
    ('dmmp_gb_uv', [r'against sarin\'s \$\+?([\d.]+)\$'], 0.0015),
    # 도판 (d) 에 박힌 p 도 같은 사실표와 대조한다.
    ('p_struct_fig', [r'A_rdkit[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_elec_fig', [r'B_electronic[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_merged_fig', [r'C_merged[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('pct_struct', [r'value of 0.028 sits at the ([\d.]+)(?:nd|th|st|rd) percentile'], 0.06),  # 2000회 분해능 0.05
    ('pct_elec', [r'electronic value of 0.275 at the ([\d.]+)th'], 0.06),  # 2000회 분해능 0.05
    ('pct_pc3', [r'value of 0\.0044 sits at the ([\d.]+)(?:st|nd|rd|th)\s*\n?percentile'], 0.06),  # 2000회 분해능 0.05
    ('mantel_minq', [r'smallest adjusted value is\s*\n?([\d.]+) and no cell',
                     r'smallest \$q\$ = ([\d.]+)\)'], 0.0005),
    ('white_frac28', [r'from 36.3\\% to ([\d.]+)\\%'], 0.05),
    ('white_med28', [r'from 0.094 to ([\d.]+)\.'], 0.0005),
    ('q_by_space28', [r'over 6, per space & [\d.]+ & [\d.]+ & ([\d.]+)'], 0.0005),
    ('q_prespec28', [r'would survive at\s*\n?\$q=([\d.]+)\$'], 0.0005),
    ('meff28', [r'estimates ([\d.]+) effective tests'], 0.005),
    ('meff_p28', [r'and gives\s*\n?([\d.]+), so the approximation'], 0.0005),
    ('match35_inc', [r'twenty-eight-compound panel and (\w+) on the primary one'], 0),
    ('match35_lit', [r'Admitting the\s*\n?five as simulants gives (\w+) matched agents'], 0),
    ('match28_inc', [r'agents of eleven, against\s*\n?(\w+) on the twenty-eight-compound panel'], 0),
    ('dfp_gb', [r'fluorophosphate is present, at ([\d.]+),'], 0.005),
    ('dfp_gd', [r'is present, at [\d.]+, ([\d.]+) and'], 0.005),
    ('dfp_gf', [r'is present, at [\d.]+, [\d.]+ and ([\d.]+) from'], 0.005),
    ('dmmp_gb_35', [r'1.14 typical nearest-neighbour separations from DMMP on the\s*\n?smaller panel and\s*\n?([\d.]+) on the larger'], 0.005),
    ('cees_hd', [r'sulfide at ([\d.]+) from sulfur mustard'], 0.005),
    # 9.625 의 소수 둘째자리 반올림 경계이므로 허용오차를 0.006 으로 둔다
    ('null_exp35', [r'expectation of ([\d.]+) set by the 21:3'], 0.006),
    ('null_exp28', [r'smaller panel against ([\d.]+)'], 0.005),
    ('wy35', [r'Westfall-Young min-\$P\$ over 18\s+& ([\d.]+) &'], 0.0005),
    ('wy28', [r'Westfall-Young min-\$P\$ over 18\s+& [\d.]+ & [\d.]+ & ([\d.]+)'], 0.0005),
    ('wy35lit', [r'Westfall-Young min-\$P\$ over 18\s+& [\d.]+ & ([\d.]+) &'], 0.0005),
    ('boot_inc_zero', [r'exactly zero\s*\n?in ([\d.]+)\\% of replicates'], 0.05),
    ('boot_diff_pos', [r'positive in (\d+) of 1000'], 0),
    ('ret_gb3', [r'as its nearest compound in ([\d.]+)\\% of replicates'], 0.05),
    ('ret_gd3', [r'replicates, soman in ([\d.]+)\\%'], 0.05),
    ('ret_gf3', [r'and\s*\n?cyclosarin in ([\d.]+)\\%'], 0.05),
    ('ret_ga8', [r'cyanophosphonate in ([\d.]+)\\%'], 0.05),
    ('phorate_fold', [r'is predicted (\d[\d,]*) times\s*\n?too safe'], 1),
    ('malathion_fold', [r'is predicted (\d[\d,]*) times too toxic'], 1),
]

WORDS = {'zero': 0, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6,
         'seven': 7, 'eight': 8, 'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12,
         'thirteen': 13, 'fourteen': 14, 'fifteen': 15, 'sixteen': 16,
         'seventeen': 17, 'eighteen': 18, 'nineteen': 19, 'twenty': 20}


def tonum(s):
    s = s.strip()
    if s.lower() in WORDS:
        return float(WORDS[s.lower()])
    # LaTeX 천단위 구분자(\,)와 쉼표를 모두 제거한다. 순서를 뒤집으면
    # '41\\,685' 가 '41\\685' 로 남아 파싱이 깨진다.
    return float(s.replace('\\,', '').replace(',', ''))


def main():
    f = load_facts()
    def strip_comments(s):
        # LaTeX 주석은 PDF 에 들어가지 않는다. 주석을 읽으면 '보관용' 문구가
        # 검증 통과의 그림자 출처가 되므로 대조 전에 제거한다. \% 는 보존.
        out = []
        for line in s.split('\n'):
            if line.lstrip().startswith('%'):
                continue
            out.append(re.sub(r'(?<!\\)%.*$', '', line))
        return '\n'.join(out)

    tex = strip_comments(open(TEX).read())
    body = tex[tex.index('\\begin{abstract}'):]
    # 2026-10-13. 종전에는 main.tex 만 읽었다. 그래서 커버레터와 하이라이트가
    # 본문과 다른 Spearman·후보수를 싣고도 PASS 가 나왔다. 투고 패키지에 들어가는
    # 글은 전부 같은 산출물과 대조한다.
    import os
    for extra in ('paper/cover_letter.tex', 'paper/highlights.tex',
                  'paper/supplementary.tex'):
        if os.path.exists(extra):
            body += '\n' + strip_comments(open(extra).read())
    # 도판 주석도 같은 규칙으로 검사한다. 본문만 보던 설계 때문에 Figure 2 의
    # 구판 문자열('the least acutely toxic of the thirteen')이 통과했다.
    import fig_annot
    annot = fig_annot.all_text()
    if not annot:
        print('[verify] 경고: outputs/figure_annotations.json 이 비어 있다. '
              '도판 스크립트를 먼저 실행해야 주석이 검사된다.')
    body = body + '\n' + annot
    fails, checked, unseen = [], 0, []

    # (D) 도판 신선도. 2026-10-13 신설. 도판은 outputs/ 에 생성되고 paper/ 로
    # 복사되는데 그 복사가 수작업이었다. 분석을 다시 돌린 뒤 복사를 잊으면
    # 원고에 옛 그림이 박히고, 수치 대조는 전부 통과한다. GATE 8 3차에서
    # 도판과 본문이 어긋난 경로가 이것이므로 시각을 직접 비교한다.
    FIGS = ['Fig_main.png', 'Fig_panel.png', 'Fig_toxicity.png',
            'Fig_qsar_failure.png', 'Fig_h1.png']
    for fig in FIGS:
        src, dst = f'outputs/{fig}', f'paper/{fig}'
        if not os.path.exists(src):
            unseen.append(f'{fig}: outputs/ 에 없음 (도판 스크립트 미실행)')
            continue
        if not os.path.exists(dst):
            fails.append(f'  도판 누락 {fig}: paper/ 에 없다')
            continue
        if os.path.getmtime(dst) < os.path.getmtime(src) - 1:
            fails.append(f'  도판 구판 {fig}: paper/ 사본이 outputs/ 보다 '
                         f'오래됐다 (sh build_paper.sh 로 복사할 것)')
        checked += 1
    # 분석 산출물보다 도판이 오래되면 그림이 낡은 수치를 그리고 있다.
    newest = max(os.path.getmtime(x) for x in
                 ['outputs/experiment_summary.json', 'outputs/rollback_c.json']
                 if os.path.exists(x))
    for fig in FIGS:
        src = f'outputs/{fig}'
        if os.path.exists(src) and os.path.getmtime(src) < newest - 1:
            fails.append(f'  도판 구판 {fig}: 분석 산출물보다 오래됐다 '
                         f'(sh build_paper.sh --figs 로 재생성할 것)')

    for key, pats, tol in CHECKS:
        if key not in f:
            unseen.append(f'{key}: 사실표에 없음')
            continue
        want = float(f[key])
        found = []
        for pat in pats:
            for m in re.finditer(pat, body):
                found.append((tonum(m.group(1)), m.group(0)[:62].replace('\n', ' ')))
        if not found:
            unseen.append(f'{key}: 본문에서 문맥을 찾지 못함 (값 {want:g})')
            continue
        for got, ctx in found:
            checked += 1
            if abs(got - want) > tol:
                fails.append(f'  불일치 {key}: 본문 {got:g} vs 산출물 {want:g}'
                             f'   [{ctx}]')

    # (C) 데이터가 반증하는 문구. 수치가 맞아도 주장이 틀릴 수 있다.
    # Figure 2 의 구판 주석은 숫자를 하나도 틀리지 않으면서
    # 'the least acutely toxic of the thirteen' 이라고 주장했다.
    contradictions = []
    if not f.get('_dcm_is_safest', True):
        # 부정·정정 문맥은 제외한다. '더 이상 최소독성이 아니다' 같은 서술은
        # 반증된 주장이 아니라 그 주장을 거둔 문장이다.
        NEG = re.compile(r'no longer|is not the|are not the|rather than|'
                         r'but is not|corrected|withdraw', re.I)
        for m in re.finditer(r'[^.\n]*least\s+acutely\s+toxic[^.]*\.', body):
            seg = m.group(0)
            if NEG.search(seg):
                continue
            if 'dichloromethane' in seg or 'simulant' in seg:
                contradictions.append(
                    '  반증된 주장: 본문/도판이 dichloromethane 을 '
                    f'최소독성이라 하지만 실제 순위는 {f["dcm_panel_rank"]}위이고 '
                    f'{", ".join(f["_panel_safer_than_dcm"])} 가 더 안전하다.'
                    f'  [{" ".join(seg.split())[:90]}]')

    # (B) 근사 불일치: 사실표 값에 0.5% 안쪽으로 가깝지만 같지 않은 수치
    # 근사 검사는 1000 이상의 값에만 적용한다. 반올림 분기가 실제로 문제가 된
    # 것은 천단위 구분자가 붙는 환산값이었고(41,683 대 41,685), 그 아래에서는
    # 서로 무관한 양끼리 우연히 0.5% 안쪽으로 가까워져 오경보만 낸다
    # ('270 chemical entries' 와 malathion 의 268.8배).
    targets = {k: float(v) for k, v in f.items()
               if isinstance(v, (int, float)) and abs(float(v)) >= 1000}
    near = []
    for m in re.finditer(r'(?<![\d.])(\d[\d,]*(?:\\,\d{3})*(?:\.\d+)?)', body):
        try:
            v = tonum(m.group(1))
        except ValueError:
            continue
        if v < 1000:
            continue
        raw = m.group(1)
        dec = len(raw.split('.')[1]) if '.' in raw else 0
        sig = len(re.sub(r'[^0-9]', '', raw).lstrip('0')) or 1
        for k, t in targets.items():
            if t == 0:
                continue
            # 본문이 쓴 정밀도로 반올림해 일치하면 정상 표기다.
            # 41,684.91 을 41,685 로 적는 것은 결함이 아니고,
            # 41,683 으로 적는 것은 결함이다.
            if round(t, dec) == round(v, dec):
                continue
            try:
                from decimal import Decimal
                if float(f'%.{sig}g' % t) == v:
                    continue
            except Exception:
                pass
            d = abs(v - t) / t
            if 0 < d <= 0.005:
                near.append(f'  근사 불일치 {k}: 본문 {v:,.4g} vs 산출물 {t:,.4g} '
                            f'(차이 {100*d:.3f}%)  [...{body[max(0,m.start()-40):m.end()+12]}]'
                            .replace('\n', ' '))
    near = sorted(set(near))

    print(f'[verify] 문맥 대조 {checked}건')
    if unseen:
        print(f'[verify] 문맥 미발견 {len(unseen)}건 (검사 누락 — 정규식 보강 필요)')
        for u in unseen:
            print('  ', u)
    if fails:
        print(f'[verify] 불일치 {len(fails)}건')
        print('\n'.join(fails))
    if near:
        print(f'[verify] 근사 불일치 {len(near)}건 (반올림 분기 의심)')
        print('\n'.join(near))
    if contradictions:
        print(f'[verify] 데이터가 반증하는 주장 {len(contradictions)}건')
        print('\n'.join(sorted(set(contradictions))))
    ok = not fails and not near and not contradictions
    print(f'[verify] 판정: {"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
