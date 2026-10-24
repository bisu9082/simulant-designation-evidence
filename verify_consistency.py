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

# 사실표가 읽는 산출물 전부 (2026-10-17 신설).
#
# 왜 필요한가. 저장소를 새로 클론해 검증기를 돌리자 outputs/mde.json 하나가
# 없다는 FileNotFoundError 로 즉시 죽었고, 그 뒤에 또 무엇이 없는지는 알 수
# 없었다. 실제로는 산출물 4건과 도판 1건이 빠져 있었다. 한 건씩 터지는 대신
# 빠진 것을 한 번에 전부 보고한다. 이 목록은 '저장소가 담아야 하는 산출물'의
# 계약서 구실도 한다.
REQUIRED_OUTPUTS = [
    'outputs/aegl_parse_validation.json',
    'outputs/blockbalance.json',
    'outputs/dim_percentiles.json',
    'outputs/distance_table_normalized.csv',
    'outputs/experiment_summary.json',
    'outputs/figure_annotations.json',
    'outputs/h1_candidates.csv',
    'outputs/h1_summary.json',
    'outputs/h1_screen_sensitivity.json',
    'outputs/match_by_panel.json',
    'outputs/mde.json',
    'outputs/minority_regression.json',
    'outputs/multiplicity_c.json',
    'outputs/permutation_sidedness.json',
    'outputs/pesticide_panel.json',
    'outputs/cite_facts.json',
    'outputs/pesticide_family.json',
    'outputs/qsar_heldout.csv',
    'outputs/qsar_tox_summary.json',
    'outputs/rollback_c.json',
    'outputs/stat_supplement.json',
    'outputs/stat_supplement2.json',
    'outputs/toxicity_aegl.csv',
    'outputs/toxicity_predicted.csv',
]

# 도판 스크립트가 추가로 읽는 것. 검증기 자체는 쓰지 않으나 이것이 없으면
# build_paper.sh --figs 가 깨지고, 그러면 원고에 구판 그림이 남는다.
REQUIRED_FOR_FIGURES = [
    'outputs/_electronic_with_gap.csv',
    'outputs/_pca_scores.npy',
    'outputs/_null_A_rdkit.npy',
    'outputs/_null_B_electronic.npy',
    'outputs/_null_C_merged.npy',
    'outputs/electronic_distance_table.csv',
    'outputs/electronic_space_summary.json',
    'outputs/h1_gate_tiers.json',
    'outputs/sanity_checks.json',
    'outputs/fig_label_counts.json',
    'outputs/nearest_agent_bootstrap.json',
    'outputs/null_match_expectation.json',
]

FIGURES = ['Fig_main.png', 'Fig_spaces.png', 'Fig_toxicity.png',
           'Fig_plocal.png', 'Fig_h1.png', 'Fig_mde.png']

# 도판 신선도 검사의 mtime 허용오차(초). 2026-10-22. 종전 1초 기준은 갓 클론한
# 저장소에서 거짓으로 걸린다. git clone 은 모든 파일의 mtime 을 체크아웃 시각으로
# 새로 쓰는데, 13 MB 짜리 outputs/ 를 풀다 보면 초 경계를 넘어가기 때문이다.
# 논문이 인용하는 v1.0 을 리뷰어가 클론해서 바로 이 검증기를 돌린다는 점에서
# 거짓 실패는 치명적이다. 실제 구판은 날 단위로 벌어지므로(그래픽 초록 적발
# 당시 11일) 120초면 둘을 가른다.
MTIME_TOL = 120


def preflight():
    """산출물·도판 존재를 먼저 전수 확인한다. 없으면 전부 나열하고 멈춘다."""
    import os
    missing = [p for p in REQUIRED_OUTPUTS if not os.path.exists(p)]
    figless = [p for p in REQUIRED_FOR_FIGURES if not os.path.exists(p)]
    nofig = [f for f in FIGURES if not os.path.exists(f'paper/{f}')]
    # 그래픽 초록은 build_paper.sh 바깥에 있어서 분석을 다시 돌려도 갱신되지
    # 않았다. 2026-10-22 투고 직전에 p=0.297, Spearman r=-0.022, n=214 를
    # 그대로 들고 있는 10-11 판이 발견됐다. 산출물보다 오래되면 멈춘다.
    if os.path.exists('outputs/graphical_abstract.png'):
        _ga = os.path.getmtime('outputs/graphical_abstract.png')
        _src = max(os.path.getmtime(p) for p in
                   ('outputs/sanity_checks.json', 'outputs/h1_summary.json',
                    'outputs/h1_candidates.csv', 'outputs/stat_supplement.json')
                   if os.path.exists(p))
        if _ga < _src - MTIME_TOL:
            nofig.append('graphical_abstract.png (산출물보다 오래됐다 — '
                         'python3 make_ga.py 로 재생성할 것)')
    # 2026-10-20 신설. figure_annotations.json 은 누적만 하므로, 폐기된 도판의
    # 주석이 남아 본문 검사를 대신 통과시킬 수 있다. 등재부에 현행 도판이 아닌
    # 항목이 하나라도 있으면 멈춘다.
    stale = []
    if os.path.exists('outputs/figure_annotations.json'):
        import fig_annot as _fa
        # graphical_abstract 는 원고에 실리는 도판이 아니라 투고 시스템에
        # 따로 올리는 첨부다. FIGURES 에는 없지만 등재부는 쓴다 (2026-10-22).
        live = {f[:-4] for f in FIGURES} | {'graphical_abstract'}
        stale = sorted(set(_fa.load()) - live)
    if missing or figless or nofig or stale:
        print('[verify] 사전 점검 실패 — 필요한 산출물이 없다.')
        for p in missing:
            print(f'   사실표 필수 누락: {p}')
        for p in figless:
            print(f'   도판 생성 필수 누락: {p}')
        for f in nofig:
            print(f'   원고 도판 누락: paper/{f}')
        for f in stale:
            print(f'   도판 주석 등재부에 폐기 항목: {f} '
                  f'(outputs/figure_annotations.json 을 지우고 도판 재생성)')
        print('[verify] sh run_analysis.sh && sh build_paper.sh --figs 를 먼저 '
              '실행하거나, 저장소에서 빠진 파일을 추가할 것.')
        return False
    return True


def load_facts():
    """산출물에서만 읽는다. 여기에 손으로 적은 수치는 없다."""
    f = {}
    f['ga_q28'] = json.load(open('outputs/stat_supplement.json'))['family'][
        'panel28_incumbent3']['bh_all_min_q']
    NE = json.load(open('outputs/null_match_expectation.json'))['sets']
    f['nullexp_k3'] = NE['incumbent3']['null_mean']
    f['nullexp_k8'] = NE['incumbent_plus_lit8']['null_mean']
    NA = json.load(open('outputs/nearest_agent_bootstrap.json'))['simulants']
    f['na_dmmp_top'] = NA['DMMP']['top_pct']
    f['na_dmmp_g'] = NA['DMMP']['g_series_pct']
    f['na_dcm_g'] = NA['dichloromethane']['g_series_pct']
    f['na_dpgme_top'] = NA['DPGME']['top_pct']
    f['na_dpgme_g'] = NA['DPGME']['g_series_pct']
    L = json.load(open('outputs/fig_label_counts.json'))
    f['fig3_from_panel'] = L['n_from_panel']
    f['fig3_other_labelled'] = L['n_other_labelled']
    f['fig3_unlabelled'] = L['n_unlabelled']
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
    f['n_more_faithful_10'] = t['relaxed_10']['n_more_faithful_than_dcm']
    f['relaxed_bar'] = t['relaxed_10']['bar_mgm3']
    f['fid_cut'] = h['fidelity_cut']
    f['n_screened_distinct'] = h['n_screened_distinct']

    hs = json.load(open('outputs/h1_screen_sensitivity.json'))
    f['sens_kept_n'] = hs['kept']['n_candidates']
    f['sens_kept_rho'] = hs['kept']['spearman_r']

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

    md = json.load(open('outputs/mde.json'))
    f['mde_delta35'] = md['panel35_incumbent3']['mde_delta']
    f['mde_gap_median'] = md['panel35_incumbent3']['gap_median']
    f['mde_gap_min'] = md['panel35_incumbent3']['gap_min']
    f['mde_p0_35'] = md['panel35_incumbent3']['p_at_zero']
    f['mde_p0_28'] = md['panel28_incumbent3']['p_at_zero']

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

    # 2026-10-17. 표 3 (tab:match) 은 이 산출물을 쓰면서도 검증기에 등록된 적이
    # 없었다. 144건 검사를 통과하는 동안 11행 중 9행이 구판 값이었다. 동료
    # 검토에서 지적돼 전수 대조로 돌린다. 같은 이유로 res:match 산문의
    # 공간별 수치도 전부 사실표에 올린다.
    _p28 = {a['agent']: a for a in mb['panel28_incumbent3']['agents']}
    _T3 = {'GB': 'sarin (GB)', 'GD': 'soman (GD)',
           'CK': 'cyanogen chloride (CK)', 'GF': 'cyclosarin (GF)',
           'AC': 'hydrogen cyanide (AC)', 'CG': 'phosgene (CG)',
           'HD': 'sulfur mustard (HD)', 'GA': 'tabun (GA)',
           'L': 'lewisite 1 (L)', 'PS': 'chloropicrin (PS)', 'VX': 'VX'}
    for _code, _nm in _T3.items():
        f[f't3_{_code}'] = _p28[_nm]['d_norm']

    # GFN2 가상궤도 에너지의 범위. 한계 절이 이 수치로 'Koopmans 로 읽지
    # 말라'를 건다 (2026-10-23 신설).
    _x35 = pd.read_csv('descriptor_cache_xtb_35.csv')
    f['lumo_min'] = float(_x35.E_LUMO_eV.min())
    f['lumo_max'] = float(_x35.E_LUMO_eV.max())
    f['n_lumo_below4'] = int((_x35.E_LUMO_eV < -4).sum())

    # 공간별 정규화 거리와 포스겐의 최근접 경쟁자. 산문이 세 공간을 섞어
    # 적었던 자리이므로 공간을 키에 박아 둔다.
    _dt = pd.read_csv('outputs/distance_table_normalized.csv')
    _SP = {'A_rdkit': 'struct', 'B_electronic': 'elec', 'C_merged': 'merged'}
    for _cond, _g in _dt.groupby('condition'):
        _tag = _SP[_cond]
        _nn = float(_g.nn_median_space.iloc[0])
        for _, _r in _g.iterrows():
            _code = re.search(r'\((\w+)\)$', _r.real_agent)
            _code = _code.group(1) if _code else 'VX'
            f[f'{_tag}_{_code}'] = float(_r.d_sim_norm)
            if _code == 'CG':
                f[f'{_tag}_CG_rival'] = float(_r.d_interferent) / _nn

    pc = json.load(open('outputs/dim_percentiles.json'))
    f['pct_struct'] = pc['struct_full']
    f['pct_elec'] = pc['electronic']
    f['pct_pc3'] = pc['struct_pc3']

    # 블록 가중 민감도 (2026-10-17 신설, SI res:blockbalance).
    bb = json.load(open('outputs/blockbalance.json'))
    _S6 = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
           'n_closer', 'n_valid']
    _BBTAG = {'panel35_incumbent3': '35',
              'panel35_incumbent_plus_lit8': '35lit',
              'panel28_incumbent3': '28'}
    for _lab, _tag in _BBTAG.items():
        _p = bb[_lab]['block_balanced']['p']
        f[f'bb_minp{_tag}'] = min(_p[s] for s in _S6)
        f[f'bb_match{_tag}'] = bb[_lab]['block_balanced']['n_valid']
        f[f'bb_cur_match{_tag}'] = bb[_lab]['equal_per_descriptor']['n_valid']
        f[f'bb_meandiff{_tag}'] = _p['mean_diff']
    f['bb_nstruct'] = bb['panel28_incumbent3']['n_struct']
    f['bb_nelec'] = bb['panel28_incumbent3']['n_elec']

    # 블록균등 병합공간 + 기존 구조·전자 = 18검정 가족의 BH 최소 q.
    def _bh_min(ps):
        o = sorted(ps)
        m = len(o)
        q = [o[i] * m / (i + 1) for i in range(m)]
        for i in range(m - 2, -1, -1):
            q[i] = min(q[i], q[i + 1])
        return q[0]
    _rc = json.load(open('outputs/rollback_c.json'))
    _RCMAP = {'panel28_designated3': ('panel28_incumbent3', '28'),
              'panel35_designated3': ('panel35_incumbent3', '35'),
              'panel35_designated_plus_literature8':
                  ('panel35_incumbent_plus_lit8', '35lit')}
    for _lab, (_bblab, _tag) in _RCMAP.items():
        _sp = _rc['panels'][_lab]
        _ps = [_sp[s]['statistics'][k]['p']
               for s in ('structural', 'electronic3') for k in _S6]
        _ps += [bb[_bblab]['block_balanced']['p'][k] for k in _S6]
        f[f'bbq_{_tag}'] = _bh_min(_ps)
        # 사전지정 통계량 3검정도 가중에 따라 달라질 수 있으므로 같이 둔다.
        f[f'bbpre_{_tag}'] = _bh_min(
            [_sp['structural']['statistics']['mean_diff']['p'],
             _sp['electronic3']['statistics']['mean_diff']['p'],
             bb[_bblab]['block_balanced']['p']['mean_diff']])

    # 네 번째 구성 — 농약 2종을 모사체로 분류 (2026-10-20 신설).
    # 본문 §2.1 이 '쟁점의 여지가 있다'고 적고 돌리지 않았던 라벨 결정을
    # 실제로 돌린 결과다. 이 논문에서 다중검정 보정을 통과하는 유일한 구성이라
    # 수치 하나가 틀리면 주장이 통째로 바뀐다. 표와 산문을 모두 등록한다.
    pp = json.load(open('outputs/pesticide_panel.json'))
    f['pp_min_p'] = pp['min_p']
    f['pp_bh_q'] = pp['bh_min_q']
    f['pp_n_perm'] = pp['spaces']['merged']['n_perm']
    f['pp_n_sim'] = pp['n_simulants']
    for _sp, _tag in [('structural', 'st'), ('electronic3', 'el'),
                      ('merged', 'mg')]:
        _d = pp['spaces'][_sp]
        for _k, _short in [('mean_diff', 'md'), ('median_diff', 'mdn'),
                           ('mean_logratio', 'lr'), ('mean_dsim', 'ds'),
                           ('n_closer', 'nc'), ('n_valid', 'nv')]:
            f[f'pp_{_tag}_{_short}'] = _d['p'][_k]
        f[f'pp_{_tag}_match'] = _d['n_valid']

    # 외부 문헌에서 옮겨 적은 수치 (2026-10-20). 계산값이 아니라 전사값이므로
    # 어떤 산출물과도 대조되지 않는 상태에 있었다. 표 3 사고와 같은 경로다.
    cf = json.load(open('outputs/cite_facts.json'))['papers']
    f['j18_identified'] = cf['james2018']['n_simulants_identified']
    f['j18_unsuitable'] = cf['james2018']['n_unsuitable']
    f['j21_identified'] = cf['james2021']['n_chemicals_identified']
    f['j21_discounted'] = cf['james2021']['n_discounted']
    f['j21_suitable'] = cf['james2021']['n_suitable']
    f['j21_possible'] = cf['james2021']['n_possible']

    # 네 번째 구성의 다중검정 체계 5종 (2026-10-20). 표 2 네 번째 열.
    pf = json.load(open('outputs/pesticide_family.json'))
    f['pf_min_p'] = pf['min_p']
    f['pf_bh18'] = pf['bh_all_min_q']
    f['pf_bh6'] = pf['bh_by_space_min_q']
    f['pf_wy'] = pf['westfall_young_adjusted_p']
    f['pf_meff'] = pf['meff']
    f['pf_bonf_meff'] = pf['meff_adjusted_min_p']
    f['pf_prespec'] = pf['prespecified_bh_min_q']

    # 서술자 감소 정산 (2026-10-20). 본문이 51 과 38 을 따로 적고 잇지 않았다.
    f['struct_removed'] = 13
    f['struct_candidates'] = 51

    mp = json.load(open('outputs/multiplicity_c.json'))
    for lab, _m, _minp, _bonf, bhq, _k in mp:
        if lab == 'panel28_designated3':
            f['bh_q28'] = bhq
            f['minp28'] = _minp
        if lab == 'panel35_designated3':
            f['bh_q35'] = bhq
            f['minp35'] = _minp
        if lab == 'panel35_designated_plus_literature8':
            f['bh_q35lit'] = bhq
            f['minp35lit'] = _minp
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
                      r'association across the (\d[\d,]*)\s*\n?compounds',
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
    ('n_incumbent', [r'Faithful and safe at that bar & (\d+) '], 0),
    ('lumo_min', [r'they run from \$([-+][\d.]+)\$~eV for chloropicrin'], 0.005),
    ('lumo_max', [r'to \$([-+][\d.]+)\$~eV for water'], 0.005),
    ('n_lumo_below4', [r'and (twenty-four|twenty-three|twenty-five) of the\s*\n?'
                       r'thirty-five fall below'], 0),
    ('n_more_faithful_10', [r'more faithful than dichloromethane itself & (\d+) ',
                            r'and all (eight|seven|nine) are more faithful than'], 0),
    ('fid_cut', [r'own fidelity of\s*\n?0\.78 falling outside the ([\d.]+) quartile',
                 r'methanol at 2\.00 against the ([\d.]+) boundary'], 0.005),
    ('n_screened_distinct', [r'(Eight|Seven|Nine) further records were set aside',
                             r'Retaining all (eight|seven|nine)\s*\n?records'], 0),
    ('sens_kept_n', [r'raises the candidate count to (\d+) and'], 0),
    ('sens_kept_rho', [r'moves the Spearman coefficient\s*\n?from \$\+0\.062\$ to '
                       r'\$([-+][\d.]+)\$'], 0.0005),
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
    ('nullexp_k3',
     [r'against (\d+\.\d+) and \d+\.\d+ agents expected by chance'], 0.05),
    ('nullexp_k8', [r'against \d+\.\d+ and (\d+\.\d+) agents expected'], 0.05),
    ('na_dmmp_top', [r'DMMP is nearest sarin in\s*\n?(\d+\.\d+)\\% of replicates'], 0.05),
    ('na_dmmp_g',
     [r'nearest one of the four G-series agents in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dcm_g',
     [r'Dichloromethane is nearest a G-series agent in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dpgme_top',
     [r'DPGME divides evenly, sulfur mustard in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dpgme_g', [r'and the\s*\n?G series in (\d+\.\d+)\\%'], 0.05),
    ('ga_q28', [r'corrected over the 18 tests: \$q\$ = (\d+\.\d+)'], 0.0005),
    ('fig3_from_panel', [r"name the (ten|nine|eleven) agents of this paper's panel"], 0),
    ('fig3_other_labelled',
     [r'the\s*\n?(three|two|four) nitrogen mustards'], 0),
    ('fig3_unlabelled',
     [r'the remaining (seven|six|eight) held-out agents are plotted'], 0),
    ('panel_without', [r'\(b\) The (nine|fifteen|ten)\s+compounds with no\s+regulator-grade',
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
    ('panel_without', [r'inhalation\s*\n?\s*value exists for any of these (\d+)'], 0),
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
    ('meff_p28', [r'and gives\s*\n?([\d.]+), conservative here by a factor'], 0.0005),
    ('match35_inc', [r'twenty-eight-compound panel and (\w+) on the primary one'], 0),
    ('match35_lit', [r'Admitting the\s*\n?five as simulants gives (\w+) matched agents'], 0),
    ('match28_inc', [r'agents of eleven, against\s*\n?(\w+) on the twenty-eight-compound panel'], 0),
    ('dfp_gb', [r'fluorophosphate is present, at ([\d.]+),'], 0.005),
    ('dfp_gd', [r'is present, at [\d.]+, ([\d.]+) and'], 0.005),
    ('dfp_gf', [r'is present, at [\d.]+, [\d.]+ and ([\d.]+) from'], 0.005),
    ('dmmp_gb_35', [r'sarin reads 1\.14 against DMMP on the\s*\n?smaller panel and (\d+\.\d+) on the larger'], 0.005),
    ('cees_hd', [r'sulfide at ([\d.]+) from sulfur mustard'], 0.005),
    # 최소검출효과 (2026-10-17). 본문 Limitations + SI + 도판 주석 3경로 대조.
    ('mde_delta35', [r'moved ([\d.]+) typical nearest-neighbour',
                     r'returns \$p<0\.05\$ at\s*\n?\$\\delta = ([\d.]+)\$',
                     r'"mde_delta35": "([\d.]+)"'], 0.005),
    ('mde_gap_median', [r'carry has a median of\s*\n?([\d.]+)',
                        r'median of ([\d.]+), more than twice',
                        r'"gap_median35": "([\d.]+)"'], 0.0005),
    ('mde_gap_min', [r'only at sarin, where it is\s*\n?([\d.]+)',
                     r'it is ([\d.]+), does it fall below',
                     r'"gap_min35": "([\d.]+)"'], 0.0005),
    ('mde_p0_35', [r'it returns ([\d.]+) on\s*\n?the primary panel',
                   r'"p_at_zero35": "([\d.]+)"'], 0.0005),
    ('mde_p0_28', [r'and ([\d.]+) on the twenty-eight-compound panel, reproducing',
                   r'"p_at_zero28": "([\d.]+)"'], 0.0005),
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
    # 표 3 (tab:match) 전수 대조. 2026-10-17 신설. 이 표는 그때까지 검증기에
    # 등록되지 않아 11행 중 9행이 구판 값인 채로 144건 검사를 통과했다.
    ('t3_GB', [r'GB & DMMP & ([\d.]+) &', r'Sarin \(([\d.]+)\), soman'], 0.0005),
    ('t3_GD', [r'GD & DMMP & ([\d.]+) &',
               r'Sarin \([\d.]+\), soman \(([\d.]+)\)'], 0.0005),
    ('t3_CK', [r'CK & DCM & ([\d.]+) &'], 0.0005),
    ('t3_GF', [r'GF & DMMP & ([\d.]+) &',
               r'cyclosarin \(([\d.]+)\) are matched'], 0.0005),
    ('t3_AC', [r'AC & DCM & ([\d.]+) &'], 0.0005),
    ('t3_CG', [r'CG & DCM & ([\d.]+) &',
               r'merged space \([\d.]+ against ([\d.]+)\)'], 0.0005),
    ('t3_HD', [r'HD & DCM & ([\d.]+) &'], 0.0005),
    ('t3_GA', [r'GA & DMMP & ([\d.]+) &', r'tabun at ([\d.]+) and VX',
               r'passes its value of ([\d.]+)'], 0.0005),
    ('t3_L', [r'L & DCM & ([\d.]+) &'], 0.0005),
    ('t3_PS', [r'PS & DMMP & ([\d.]+) &'], 0.0005),
    ('t3_VX', [r'VX & DPGME & ([\d.]+) &', r'VX at\s*\n?([\d.]+) are not'], 0.0005),
    # res:match 산문의 공간별 수치. 종전 산문은 세 공간의 원거리와 정규화
    # 거리를 섞어 적었고 어느 값도 산출물과 대조되지 않았다.
    ('struct_GD', [r'soman at ([\d.]+)\s*\n?and cyclosarin at'], 0.0005),
    ('struct_GF', [r'cyclosarin at ([\d.]+) falling outside'], 0.0005),
    ('struct_VX', [r'and ([\d.]+) in the structural\s*\n?space that is not a match'],
     0.0005),
    ('merged_VX', [r'at ([\d.]+) typical separations in the merged space'], 0.0005),
    ('elec_AC', [r'dichloromethane at ([\d.]+),\s*\n?because acetonitrile'], 0.0005),
    ('struct_CG', [r'structural space \([\d.]+ typical separations\s*\n?against ([\d.]+)\)'],
     0.005),
    ('minp35', [r'Smallest nominal \$p\$ \(18 tests\)\s+& ([\d.]+) &',
                r'eighteen-test family there is ([\d.]+)\s*\n?with equal weight per'],
     0.0005),
    ('minp35lit', [r'Smallest nominal \$p\$ \(18 tests\)\s+& [\d.]+ & ([\d.]+) &'],
     0.0005),
    ('minp28', [r'Smallest nominal \$p\$ \(18 tests\)\s+& [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('elec_CG', [r'electronic space \([\d.]+ against ([\d.]+)\)'], 0.005),
    ('struct_CG_rival', [r'structural space \(([\d.]+) typical separations'], 0.005),
    ('merged_CG_rival', [r'ethyl acetate in the merged space \(([\d.]+) against'], 0.005),
    ('elec_CG_rival', [r'chloroform in the electronic space \(([\d.]+) against'], 0.005),
    # 블록 가중 민감도 (2026-10-17). SI 표 + SI 산문 + 본문 한 문장 3경로.
    ('bb_nstruct', [r'Weighting the (\d+) structural', r'\$n_\{s\}=(\d+)\$',
                    r'puts (\d+) structural coordinates'], 0),
    ('bb_nelec', [r'and (\d+) electronic descriptors as two equal blocks',
                  r'\$n_\{e\}=(\d+)\$'], 0),
    ('bb_minp35', [r'lowers that value to ([\d.]+)',
                   r'Thirty-five, incumbent three & \d+/11 & [\d.]+ & ([\d.]+) &',
                   r'0\.151 with equal weight per descriptor and ([\d.]+) with'],
     0.0005),
    ('bb_minp35lit',
     [r'Thirty-five, plus literature eight & \d+/11 & [\d.]+ & ([\d.]+) &'], 0.0005),
    ('bb_minp28',
     [r'Twenty-eight, incumbent three & \d+/11 & [\d.]+ & ([\d.]+) &'], 0.0005),
    ('bbq_35', [r'Thirty-five, incumbent three & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_35lit',
     [r'Thirty-five, plus literature eight & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_28',
     [r'Twenty-eight, incumbent three & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bb_match35', [r'Thirty-five, incumbent three & (\d+)/11'], 0),
    ('bb_match35lit', [r'Thirty-five, plus literature eight & (\d+)/11'], 0),
    ('bb_match28', [r'Twenty-eight, incumbent three & (\d+)/11'], 0),
    ('bb_meandiff28', [r'there falls from 0\.025 to ([\d.]+)'], 0.0005),
    ('bbpre_28', [r'which stays at ([\d.]+) on the twenty-eight-compound panel'],
     0.0005),
    # 네 번째 구성 (2026-10-20). 표 1 네 번째 블록 + §3.1 산문 + Limitations.
    # 표 3 사건의 교훈대로 표를 먼저 등록한다.
    ('pp_st_md', [r'35 cpd\., incumbent 3 & structural  & 38 & \\textbf\{([\d.]+)\}'],
     0.0005),
    ('pp_st_lr', [r'35 cpd\., incumbent 3 & structural  & 38 & \\textbf\{[\d.]+\} & '
                  r'\\textbf\{\$<\$0\.001\} & ([\d.]+) &'], 0.0005),
    ('pp_st_ds', [r'35 cpd\., incumbent 3 & structural  & 38 & \\textbf\{[\d.]+\} & '
                  r'\\textbf\{\$<\$0\.001\} & [\d.]+ & ([\d.]+) &'], 0.0005),
    ('pp_el_md', [r'\\quad \+ literature 5 & electronic  &  3 & ([\d.]+) &'], 0.0005),
    ('pp_el_mdn', [r'\\quad \+ literature 5 & electronic  &  3 & [\d.]+ & ([\d.]+) &'],
     0.0005),
    ('pp_el_nv', [r'literature 5 & electronic  &  3 & '
                  r'(?:[\d.]+ & ){5}([\d.]+) &'], 0.0005),
    ('pp_mg_md', [r'\\quad \+ pesticides 2 & merged      & 41 & \\textbf\{([\d.]+)\}'],
     0.0005),
    ('pp_mg_lr', [r'\\quad \+ pesticides 2 & merged      & 41 & \\textbf\{[\d.]+\} & '
                  r'\\textbf\{\$<\$0\.001\} & \\textbf\{([\d.]+)\}'], 0.0005),
    ('pp_mg_ds', [r'\\quad \+ pesticides 2 & merged      & 41 & \\textbf\{[\d.]+\} & '
                  r'\\textbf\{\$<\$0\.001\} & \\textbf\{[\d.]+\} & ([\d.]+) &'], 0.0005),
    ('pp_n_perm', [r'\\binom\{24\}\{10\}=([\d\\,]+)'], 0),
    ('pp_min_p', [r'smallest nominal value of ([\d.]+) and a'], 0.00001),
    ('pf_meff', [r'\$m_\{\\mathrm\{eff\}\}\$ of 18 \(Li-Ji\)     & [\d.]+  & [\d.]+  & [\d.]+  & ([\d.]+)'], 0.005),
    ('pf_prespec', [r'Statistic specified in advance \(3 tests\) & [\d.]+ & [\d.]+ & '
                    r'\\textbf\{[\d.]+\} & \\textbf\{([\d.]+)\}'], 0.0005),
    ('struct_removed', [r'removes (\d+) of the 51'], 0),
    ('struct_candidates', [r'removes \d+ of the (\d+)',
                           r'descriptors \((\d+) after excluding run metadata\)'], 0),
    ('pp_bh_q', [r'Benjamini-Hochberg value of ([\d.]+) over the same',
                 r'produces an adjusted\s*\n?value of (\d+\.\d+)',
                 r'into one with an adjusted value of\s*\n?(\d+\.\d+)'], 0.00005),
    # 인용 문헌의 수치 (2026-10-20). outputs/cite_facts.json 이 단일 출처다.
    ('j18_identified', [r'identified\s*\n?(\d+) simulants'], 0),
    ('j18_unsuitable', [r'judged (\d+) unsuitable'], 0),
    ('j21_identified', [r'a later one screened (\d+) chemicals'], 0),
    ('j21_discounted', [r'discounted (\d+) on the risk'], 0),
    ('j21_suitable', [r'repository of (\d+) suitable'], 0),
    ('j21_possible', [r'suitable and (\d+) provisional'], 0),
]

WORDS = {'zero': 0, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6,
         'seven': 7, 'eight': 8, 'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12,
         'thirteen': 13, 'fourteen': 14, 'fifteen': 15, 'sixteen': 16,
         'seventeen': 17, 'eighteen': 18, 'nineteen': 19, 'twenty': 20,
         'twenty-three': 23, 'twenty-four': 24, 'twenty-five': 25}


def tonum(s):
    s = s.strip()
    if s.lower() in WORDS:
        return float(WORDS[s.lower()])
    # LaTeX 천단위 구분자(\,)와 쉼표를 모두 제거한다. 순서를 뒤집으면
    # '41\\,685' 가 '41\\685' 로 남아 파싱이 깨진다.
    return float(s.replace('\\,', '').replace(',', ''))


def main():
    if not preflight():
        return 1
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
    FIGS = ['Fig_main.png', 'Fig_spaces.png', 'Fig_toxicity.png',
            'Fig_plocal.png', 'Fig_h1.png']
    for fig in FIGS:
        src, dst = f'outputs/{fig}', f'paper/{fig}'
        if not os.path.exists(src):
            unseen.append(f'{fig}: outputs/ 에 없음 (도판 스크립트 미실행)')
            continue
        if not os.path.exists(dst):
            fails.append(f'  도판 누락 {fig}: paper/ 에 없다')
            continue
        if os.path.getmtime(dst) < os.path.getmtime(src) - MTIME_TOL:
            fails.append(f'  도판 구판 {fig}: paper/ 사본이 outputs/ 보다 '
                         f'오래됐다 (sh build_paper.sh 로 복사할 것)')
        checked += 1
    # 분석 산출물보다 도판이 오래되면 그림이 낡은 수치를 그리고 있다.
    newest = max(os.path.getmtime(x) for x in
                 ['outputs/experiment_summary.json', 'outputs/rollback_c.json']
                 if os.path.exists(x))
    for fig in FIGS:
        src = f'outputs/{fig}'
        if os.path.exists(src) and os.path.getmtime(src) < newest - MTIME_TOL:
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
