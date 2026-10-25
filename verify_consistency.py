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
import os
import re
import numpy as _np
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
    'outputs/pca_facts.json',
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
    f['xreader_pairs'] = v['cross_reader_pairs']
    f['xreader_mismatch'] = v['cross_reader_mismatch']
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
        if m['dataset'] == 'd2' and m['mode'] == 'uv' and m['space'] == 'merged':
            f['mantel_d2uv_merged'] = m['p_mantel']
            f['mantel_d2uv_merged_rho'] = m['spearman']
        if m['dataset'] == 'd2' and m['mode'] == 'normal' and \
                m['space'] == 'electronic3':
            f['mantel_d2nor_elec'] = m['p_mantel']
            f['mantel_d2nor_elec_rho'] = m['spearman']
            f['naive_d2nor_elec'] = m['p_pairs_as_independent']
    # 2026-10-24 신설. 중복 채널 축약은 본문이 '23 effective channels' 라고
    # 쓰는 새 수치다. 손으로 적은 숫자가 아니라 산출물에서 온 것임을 묶는다.
    # 화합물 단위 검정력. deltaE_power.py 와 같은 식을 쓴다.
    from scipy import stats as _st
    from scipy.optimize import brentq as _bq
    def _pw(r, n, a=0.05):
        z = _np.arctanh(r) * _np.sqrt(n - 3)
        zc = _st.norm.ppf(1 - a / 2)
        return _st.norm.cdf(z - zc) + _st.norm.cdf(-z - zc)
    # 2026-10-24 2차 교정. Fisher z 근사(0.924)는 실제로 돌린 정확 Mantel 의
    # 검정력을 과소평가한다. 시뮬레이션 값(mantel_power.py)이 정본이다.
    f['mde_rho_compound'] = float(_bq(lambda r: _pw(r, 6) - 0.8, .01, .999))
    import os as _os
    if _os.path.exists('outputs/mantel_power.json'):
        _mp = json.load(open('outputs/mantel_power.json'))
        f['mde_rho_mantel'] = _mp['mde_rho_mantel']
        f['mantel_power_perm'] = _mp['n_perm']
    # 2026-10-25. 전자공간을 B3LYP 궤도에너지로 갈아끼운 대조. 심사
    # 시늄레이션이 제기하고 재현으로 확인된 축이다. 기준 재현(GFN2)이 정본
    # 0.275 와 1/11 을 내는 것까지 같은 파일에 들어 있다.
    if _os.path.exists('outputs/dft_electronic_space.json'):
        _es = json.load(open('outputs/dft_electronic_space.json'))
        f['dft_es_gfn2_p'] = float(_es['gfn2']['p_one_sided'])
        f['dft_es_gfn2_match'] = int(_es['gfn2']['n_match'])
        f['dft_es_b3lyp_p'] = float(_es['b3lyp']['p_one_sided'])
        f['dft_es_b3lyp_match'] = int(_es['b3lyp']['n_match'])

    # 2026-10-25. 병리지수 제거 대조. 구판은 omega 가 아직 든 4좌표 전자공간에서
    # 돌아간 산출물이었고 본문은 그 값을 '모든 분석을 반복했다'며 인용했다.
    if _os.path.exists('outputs/permutation_nopatho.json'):
        _nop = json.load(open('outputs/permutation_nopatho.json'))
        if _nop['B_electronic']['n_kept'] != 3:
            raise SystemExit('[verify] 중단: permutation_nopatho 의 전자공간이 '
                             f"{_nop['B_electronic']['n_kept']}좌표다. omega 제거 "
                             '이후 재실행되지 않은 산출물이다.')
        f['nopatho_struct_p'] = float(_nop['A_rdkit']['p_one_sided'])
        f['nopatho_elec_p'] = float(_nop['B_electronic']['p_one_sided'])
        f['nopatho_merged_p'] = float(_nop['C_merged']['p_one_sided'])

    # 2026-10-25 신설. DFT 대조 수치가 SI 에 들어갔다. GFN2 의 omega 축에서
    # 끌어낸 순서가 방법 종속임을 보이는 값들이므로 산출물에 묶는다.
    if _os.path.exists('outputs/dft_vs_gfn2.csv'):
        import pandas as _pd
        _dv = _pd.read_csv('outputs/dft_vs_gfn2.csv')
        _ag = _dv[_dv.compound_class == 'real_agent']
        _it = _dv[_dv.compound_class == 'interferent']
        _gb = float(_ag[_ag.name.str.startswith('sarin')].omega_eV.iloc[0])
        _dm = float(_dv[_dv.name == 'DMMP'].omega_eV.iloc[0])
        f['dft_omega_min'] = float(_dv.omega_eV.min())
        f['dft_omega_max'] = float(_dv.omega_eV.max())
        f['gfn2_omega_min'] = float(_dv.gfn2_w.min())
        f['gfn2_omega_max'] = float(_dv.gfn2_w.max())
        f['dft_ps'] = float(_dv[_dv.name.str.startswith('chloropicrin')].omega_eV.iloc[0])
        f['dft_cg'] = float(_dv[_dv.name.str.startswith('phosgene')].omega_eV.iloc[0])
        f['dft_gb'] = _gb
        f['dft_dmmp'] = _dm
        f['dft_itf_above_gb'] = int((_it.omega_eV > _gb).sum())
        f['dft_n_dft'] = int(len(_dv))
        _lo = _ag.sort_values('omega_eV').head(3)
        f['dft_lo1'] = float(_lo.omega_eV.iloc[0])
        f['dft_lo2'] = float(_lo.omega_eV.iloc[1])
        f['dft_lo3'] = float(_lo.omega_eV.iloc[2])
        f['dft_rho_omega'] = float(_dv[['gfn2_w', 'omega_eV']]
                                   .corr(method='spearman').iloc[0, 1])
        f['dft_rho_lumo'] = float(_dv[['E_LUMO_eV', 'lumo_eV']]
                                  .corr(method='spearman').iloc[0, 1])

    # 2026-10-24. 다른 세션이 올린 AsCl3 민감도 표(SI)가 검증기에 등록돼 있지
    # 않았다. 그 표의 28종 행은 이번 정확종속 규칙으로 3/11 에서 4/11 로
    # 움직였고, 등록돼 있었다면 기계가 잡았을 자리다.
    if _os.path.exists('outputs/match_by_panel_ascl3.json'):
        _ac = json.load(open('outputs/match_by_panel_ascl3.json'))
        _AC = {'panel28_incumbent3': 'a28', 'panel35_incumbent3': 'a35',
               'panel35_incumbent_plus_lit8': 'a35l',
               'panel36_incumbent3': 'a36', 'panel36_incumbent_plus_lit': 'a36l'}
        for _k, _t in _AC.items():
            if _k in _ac:
                f[f'ascl3_{_t}_match'] = _ac[_k]['n_match']
                _L = [x for x in _ac[_k]['agents'] if 'lewisite' in x['agent']]
                if _L:
                    f[f'ascl3_{_t}_dist'] = _L[0]['d_norm']

    # 2026-10-24 신설. 하단꼬리와 중첩 조건부 귀무는 초록까지 올라간 새 수치다.
    _nm = json.load(open('outputs/null_match_expectation.json'))
    f['null3_lower'] = _nm['sets']['incumbent3']['p_at_most_observed']
    f['null8_nested_mean'] = _nm['nested_null_8']['null_mean']
    f['null8_nested_p'] = _nm['nested_null_8']['p_at_least_observed']
    _ch = rb['deltaE'].get('channels', {})
    if 'd2.uv' in _ch:
        f['n_eff_d2uv'] = _ch['d2.uv']['n_effective']
        f['n_declared_dyes'] = _ch['d2.uv']['n_declared']
    if 'd3.uv' in _ch:
        f['n_eff_d3uv'] = _ch['d3.uv']['n_effective']
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
    f['mantel_cells'] = rb['deltaE']['mantel_bh']['n_tests']
    f['mantel_nominal'] = sum(1 for _m in rb['deltaE']['mantel']
                              if _m['p_mantel'] < 0.05)

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
    f['dmmp_gd_35'] = _inc['soman (GD)']['d_norm']
    f['dmmp_gf_35'] = _inc['cyclosarin (GF)']['d_norm']
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

    # 공간별 주성분 분산. 본문 4.1절 수치가 파이프라인 밖 폐기본에서 오던 것을
    # 2026-10-24 에 끊었다. 검증기가 이 값을 들고 있지 않아 못 잡았다.
    _pf = json.load(open('outputs/pca_facts.json'))
    f['pc1_struct28'] = round(100 * _pf['panel28']['A_rdkit']['pc1'], 2)
    f['pc12_elec28'] = round(100 * _pf['panel28']['B_electronic']['pc1_pc2'], 2)
    f['pc1_struct35'] = round(100 * _pf['panel35']['A_rdkit']['pc1'], 2)
    f['pc12_elec35'] = round(100 * _pf['panel35']['B_electronic']['pc1_pc2'], 2)
    f['pc12_merged28'] = round(100 * _pf['panel28']['C_merged']['pc1_pc2'], 2)
    f['pc12_merged35'] = round(100 * _pf['panel35']['C_merged']['pc1_pc2'], 2)

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
    # 2026-10-25. 종전에는 13 이 리터럴이었다. 재분석으로 28종이 36좌표가 된
    # 뒤에도 그 리터럴이 남아 본문의 '두 패널에서 같다'는 서술을 받쳐 주었다.
    # 두 패널을 따로 산출물에서 읽는다.
    _c28 = json.load(open('outputs/experiment_summary.json'))['conditions'][0]
    _c35 = json.load(open('outputs/experiment_summary_p35all.json'))['conditions'][0]
    f['struct_candidates'] = int(_c35['n_in'])
    f['struct_kept_35'] = int(_c35['n_kept'])
    f['struct_kept_28'] = int(_c28['n_kept'])
    f['struct_removed'] = int(_c35['n_in'] - _c35['n_kept'])
    f['struct_removed_28'] = int(_c28['n_in'] - _c28['n_kept'])
    if f['struct_removed'] == f['struct_removed_28']:
        raise SystemExit('[verify] 중단: 두 패널의 서술자 감소가 같아졌다. '
                         '본문 서술과 검사 설계를 함께 손봐야 한다.')

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
# 2026-10-24 신설. 표 S1 의 칸 검사는 지금까지 열 간격을 글자 그대로 박아
# 두었다 ('& 38 &', 'structural  &'). 표를 산출물에서 다시 생성하자 간격이
# 달라져 12건이 한꺼번에 '문맥 미발견'으로 꺼졌다. 꺼진 검사는 통과처럼
# 보이므로 가장 위험하다. 행 이름·공간·열 번호로 패턴을 만든다.
def row_pat(label, space, col):
    """표 S1 의 (행, 공간, col 번째 통계) 칸을 잡는 정규식."""
    skip = r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*'
    cap = r'(?:\\textbf\{)?([\d.]+)\}?'
    return (re.escape(label) + r'\s*&\s*' + space
            + r'\s*&\s*\d+\s*&\s*' + skip * col + cap)


# 본문에 인용하지 않기로 한 사실. 키마다 이유를 적는다. 비워 두는 것이 정상이고,
# 여기 들어오는 항목은 '검사를 끈다'는 선언이므로 심사 대상이다.
UNQUOTED = {}

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
    ('pc1_struct28', [r'size axis carrying ([\d.]+)\\% of the variance'], 0.06),
    ('pc12_elec28', [r'the first two\s*\n?carrying ([\d.]+)\\%'], 0.06),
    ('pc1_struct35', [r'The primary panel behaves the same way, at ([\d.]+)\\%'], 0.06),
    ('pc12_elec35', [r'the same way, at [\d.]+\\% and ([\d.]+)\\%'], 0.06),
    ('pc12_merged28', [r'([\d.]+)\\% on the smaller panel'], 0.06),
    ('pc12_merged35', [r'and ([\d.]+)\\% on the primary one'], 0.06),
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
    ('xreader_pairs', [r'agree on all (\d[\d,]*) entries for which both'], 0),
    ('trainset_checked', [r'on all (\d[\d,]*) of\s*\n?those that resolved to a structure'], 0),
    ('comma_rows', [r'(\w+[-\w]*) values reach 1000 or more'], 0),
    ('panel_from_footnote',
     [r'(Five|Nine|Four) of the nineteen values come from the',
      r'(Five|Nine|Four) of the nineteen panel values'], 0),
    ('nullexp_k3',
     [r'designating three\s*\n?non-agents at random matches (\d+\.\d+) agents',
      r'random matches 2\.68 agents on average against (\d+\.\d+) for three'], 0.05),
    ('nullexp_k8', [r'designating eight non-agents at\s*\n?random matches (\d+\.\d+) agents'], 0.05),
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
    ('dcm_panel_rank', [r'no designating\s*\n?document, ranks (\d+) of \d+'], 0),
    ('panel_with_aegl3', [r'no designating\s*\n?document, ranks \d+ of (\d+)'], 0),
    ('heldout_in_domain', [r'(\d+) of the 20 agents pass the domain check'], 0),
    ('heldout_n', [r'and all (\d+) are wrong by'], 0),
    ('heldout_worst_in_domain',
     [r'worst inside the domain is wrong by ([\d.]+) log'], 0.005),
    ('spearman_r', [r'Spearman \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?',
                    r'Spearman \$r=([-+][\d.]+)\$'], 0.0006),
    ('spearman_lo', [r'Spearman.{0,60}?95\\% CI \$?([-+][\d.]+)\$? to'], 0.001),
    ('spearman_hi', [r'Spearman.{0,60}?95\\% CI .{0,20}?to\s*\n?\$?([-+][\d.]+)\$?'], 0.001),
    ('pearson_r', [r'Pearson \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?'], 0.0006),
    # 2026-10-13: 본문이 재집필되어 종전 정규식의 문맥이 사라졌다. kappa 민감도
    # 문장과 Table 1 의 '28종 패널' 블록을 동시에 대조한다. 표에는 패널 3구성이
    # 들어 있으므로 '28 cpd.' 행에 앵커를 걸어야 한다. 앵커 없이 쓰면 35종 행의
    # 값까지 물어와 거짓 불일치가 난다.
    ('p_struct', [r'condition \(structural \$p=([\d.]+)\$',
                  row_pat('28 cpd., incumbent 3', 'structural', 0)], 0.0005),
    ('p_elec', [r'electronic \$p=([\d.]+)\$',
                row_pat('(omits the five)', 'electronic', 0)], 0.0005),
    # 각 블록의 세 번째 행은 행 이름 칸이 비어 있다. 앞 행을 닻으로 쓴다.
    ('p_merged', [r'merged \$p=([\d.]+)\$',
                  r'\(omits the five\)[^\n]*\n\s*&\s*merged\s*&\s*\d+\s*&\s*'
                  r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0005),
    ('pc3_p28', [r'from \$p=[\d.]+\$ to\s*\n?\$p=([\d.]+)\$'], 0.0005),
    ('rand3_frac28', [r'and ([\d.]+)\\% of the subspaces reach'], 0.05),
    ('rand3_med28', [r'panel the median \$p\$ is ([\d.]+)'], 0.0005),
    ('r_lumo_omega', [r'r\(E_\{\\mathrm\{LUMO\}\}, \\omega\) =\s*\n?([-+][\d.]+)'],
     0.0005),
    ('omega_distrank', [r'0\.934 on the primary panel, ([\d.]+) on the smaller one'], 0.0005),
    ('omega_pc3_var', [r'account for ([\d.]+)\\% of its variance'], 0.05),
    ('bh_q28', [r'smallest adjusted value of ([\d.]+)',
                r'Benjamini-Hochberg \$q =\s*\n?([\d.]+)\$'], 0.0005),
    ('bh_q35', [r'Benjamini-Hochberg over 18\s+& ([\d.]+) &'], 0.0005),
    ('bh_q35lit', [r'Benjamini-Hochberg over 18\s+& [\d.]+ & ([\d.]+) &'], 0.0005),
    ('mantel_cells', [r'Of (twelve|eight|ten) panel-channel-space cells',
                      r'(twelve|eight|ten) cells are one\s*\n?family'], 0),
    ('mantel_nominal', [r'(three|two|four|five) of the twelve panel-channel-space cells reach'], 0),
    # rho 를 패턴에 박아 두면 rho 가 바뀔 때마다 검사가 조용히 빠진다
    # (2026-10-24 에 두 번 겪었다). 셀을 가리키는 문구로 건다.
    ('mantel_d2uv_struct',
     [r'that same set and channel \(\$\\rho = [-+][\d.]+\$, Mantel\s*\n?'
      r'\$p = ([\d.]+)\$'], 0.0005),
    ('dft_es_gfn2_p', [r'one-sided \$p\$ of (\d+\.\d+)',
                       r'moves it from (\d+\.\d+) to \d+\.\d+'], 0.006),
    ('dft_es_gfn2_match', [r'(\w+) valid match of eleven'], 0),
    ('dft_es_b3lyp_p', [r'returns four valid matches\s*\n?and (\d+\.\d+)',
                        r'moves it from \d+\.\d+ to (\d+\.\d+)'], 0.006),
    ('dft_es_b3lyp_match', [r'returns (\w+) valid matches',
                            r'from one valid match to (\w+)'], 0),
    ('nopatho_struct_p', [r'to three decimal places \((\d+\.\d+),'], 0.006),
    ('nopatho_elec_p', [r'three decimal places \(\d+\.\d+, (\d+\.\d+) and'], 0.006),
    ('nopatho_merged_p', [r'three decimal places \(\d+\.\d+, \d+\.\d+ and (\d+\.\d+)\)'], 0.006),
    ('gfn2_omega_min', [r'they run from ([\d.]+) to\s*\n?[\d.]+~eV across this panel'], 0.006),
    ('gfn2_omega_max', [r'they run from [\d.]+ to\s*\n?([\d.]+)~eV across this panel'], 0.006),
    ('dft_omega_min', [r'range becomes ([\d.]+) to [\d.]+~eV'], 0.006),
    ('dft_omega_max', [r'range becomes [\d.]+ to ([\d.]+)~eV'], 0.006),
    ('dft_ps', [r'chloropicrin falling from [\d.]+ to ([\d.]+)'], 0.006),
    ('dft_cg', [r'phosgene from [\d.]+ to ([\d.]+)'], 0.006),
    ('dft_dmmp', [r'DMMP at ([\d.]+) against sarin'], 0.006),
    ('dft_gb', [r'DMMP at [\d.]+ against sarin at ([\d.]+)'], 0.006),
    ('dft_itf_above_gb', [r'(\w+) of the fourteen interferents still score above'], 0),
    ('dft_lo1', [r'three lowest agents become VX\s*\n?\(([\d.]+)\)'], 0.006),
    ('dft_lo2', [r'sulfur mustard \(([\d.]+)\) and cyclosarin'], 0.006),
    ('dft_lo3', [r'and cyclosarin\s*\n?\(([\d.]+)\) rather than'], 0.006),
    ('dft_rho_omega', [r'two sets of \$\\omega\$ values is ([\d.]+), against'], 0.006),
    ('dft_rho_lumo', [r'is [\d.]+, against\s*\n?([\d.]+) for the underlying LUMO'], 0.006),
    ('dft_n_dft', [r'recomputed the whole\s*\n?(\w+[-\w]*)-compound panel at'], 0),
    ('ascl3_a28_match', [r'Twenty-eight, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a28_dist', [r'Twenty-eight, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a35_match', [r'Thirty-five, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a35_dist', [r'Thirty-five, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a35l_dist', [r'Thirty-five, plus literature eight & \d+/11 & CEES & ([\d.]+)'], 0.006),
    ('ascl3_a36_match', [r'Thirty-six, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a36_dist', [r'Thirty-six, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a36l_match', [r'Thirty-six, plus literature & (\d+)/11 & AsCl'], 0),
    ('null3_lower', [r'the three at \$p=([\d.]+)\$ over all',
                     r'zero sits at \$p=([\d.]+)\$ against an expected'], 0.006),
    ('null8_nested_mean', [r'gives a mean of ([\d.]+) rather',
                           r'against\s*\n?an expected ([\d.]+) once the null'], 0.006),
    ('null8_nested_p', [r'observed five sits at \$p=([\d.]+)\$ against it',
                        r'the five at \$p=([\d.]+)\$ against'], 0.006),
    ('n_eff_d2uv', [r'leaves (\d+) effective channels',
                    r'of which (\d+) are independent'], 0),
    ('n_eff_d3uv', [r'effective channels in the second\s*\n?set'
                    r"'s ultraviolet data and (\d+) in the third"], 0),
    ('mde_rho_compound', [r'as the unit, which is the unit the Mantel test uses,\s*\n?'
                          r'returns \$([\d.]+)\$'], 0.0015),
    ('mde_rho_mantel', [r'80\\% power is reached at \$\|\\rho\| = ([\d.]+)\$',
                        r'design detects \$\|\\rho\| \\ge ([\d.]+)\$ at 80\\% power'], 0.006),
    ('mantel_power_perm', [r'all \$6! = (\d+)\$ relabellings',
                           r'exact Mantel over all (\d+) relabellings'], 0),
    ('mantel_d2nor_elec', [r'\$\\rho = \+0\.571\$,\s*\n?Mantel \$p = ([\d.]+)\$'], 0.0005),
    ('mantel_d2uv_merged_rho',
     [r"merged space in the\s*\n?second compound set's ultraviolet channel "
      r"\(\$\\rho = ([-+][\d.]+)\$"], 0.0015),
    ('mantel_d2uv_merged', [r'\\rho = \+0\.518\$, Mantel \$p = ([\d.]+)\$'], 0.0005),
    ('mantel_d2uv_struct_rho', [r'that same set and channel \(\$\\rho = ([-+][\d.]+)\$'], 0.0015),
    ('floor_d1uv_obs', [r'The observed means clear that floor by 0.013 to\s*\n?0.066: \$([-+][\d.]+)\$'], 0.0015),
    ('dmmp_ga_uv', [r'tabun, at \$\+?([\d.]+)\$ against sarin'], 0.0015),
    ('dmmp_gb_uv', [r'against sarin\'s \$\+?([\d.]+)\$'], 0.0015),
    # 도판 (d) 에 박힌 p 도 같은 사실표와 대조한다.
    ('p_struct_fig', [r'A_rdkit[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_elec_fig', [r'B_electronic[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_merged_fig', [r'C_merged[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('pct_struct', [r'-dimensional value of [\d.]+ sits at the ([\d.]+)(?:nd|th|st|rd)\s*\n?percentile'], 0.06),
    ('pct_elec', [r'electronic value of [\d.]+ (?:sits )?at the ([\d.]+)(?:nd|th|st|rd)'], 0.06),
    ('pct_pc3', [r'three-principal-component value of [\d.]+ sits at the ([\d.]+)(?:st|nd|rd|th)\s*\n?percentile'], 0.06),
    ('mantel_minq', [r'smallest adjusted value is\s*\n?([\d.]+) and no cell',
                     r'smallest \$q\$ = ([\d.]+)\)'], 0.0005),
    ('white_frac28', [r'the rate falls from [\d.]+\\% to ([\d.]+)\\%'], 0.05),
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
    ('dmmp_gb_35', [r'on that panel sarin stands\s*\n?([\d.]+) from DMMP'], 0.005),
    ('dmmp_gd_35', [r'from DMMP, soman ([\d.]+) and'], 0.005),
    ('dmmp_gf_35', [r'soman [\d.]+ and cyclosarin ([\d.]+)\.'], 0.005),
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
               r'passes its value of (\d+\.\d+)'], 0.0005),
    ('t3_L', [r'L & DCM & ([\d.]+) &'], 0.0005),
    ('t3_PS', [r'PS & DMMP & ([\d.]+) &'], 0.0005),
    ('t3_VX', [r'VX & DPGME & ([\d.]+) &', r'VX at\s*\n?([\d.]+) are not'], 0.0005),
    # res:match 산문의 공간별 수치. 종전 산문은 세 공간의 원거리와 정규화
    # 거리를 섞어 적었고 어느 값도 산출물과 대조되지 않았다.
    ('struct_GD', [r'soman at ([\d.]+)\s*\n?and cyclosarin at'], 0.0005),
    ('struct_GF', [r'cyclosarin at ([\d.]+) (?:falling )?outside'], 0.0005),
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
    ('merged_CG_rival', [r'and in the merged space \(([\d.]+) against'], 0.005),
    ('elec_CG_rival', [r'chloroform in the electronic space \(([\d.]+) against'], 0.005),
    # 블록 가중 민감도 (2026-10-17). SI 표 + SI 산문 + 본문 한 문장 3경로.
    ('bb_nstruct', [r'Weighting the (\d+) structural', r'\$n_\{s\}=(\d+)\$',
                    r'puts (\d+) structural coordinates'], 0),
    ('bb_nelec', [r'and (\d+) electronic descriptors as two equal blocks',
                  r'\$n_\{e\}=(\d+)\$'], 0),
    ('bb_minp35', [r'lowers that value to ([\d.]+)',
                   r'nominal value from 0\.151 to ([\d.]+)',
                   r'Thirty-five, incumbent three & \d+/11 & \d+/11 & [\d.]+ & ([\d.]+) &',
                   r'0\.151 with equal weight per descriptor and ([\d.]+) with'],
     0.0005),
    ('bb_minp35lit',
     [r'Thirty-five, plus literature 8 & \d+/11 & \d+/11 & [\d.]+ & ([\d.]+) &'], 0.0005),
    ('bb_minp28',
     [r'Twenty-eight, incumbent three & \d+/11 & \d+/11 & [\d.]+ & ([\d.]+) &'], 0.0005),
    ('bbq_35', [r'Thirty-five, incumbent three & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_35lit',
     [r'Thirty-five, plus literature 8 & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_28',
     [r'Twenty-eight, incumbent three & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    # 2026-10-24. 유효매칭 열을 per descriptor / per block 두 칸으로 쪼갰다.
    # 둘 다 걸지 않으면 한쪽을 망가뜨리는 주입이 통과한다 (실제로 통과했다).
    ('bb_cur_match35', [r'Thirty-five, incumbent three & (\d+)/11'], 0),
    ('bb_cur_match35lit', [r'Thirty-five, plus literature 8 & (\d+)/11'], 0),
    ('bb_cur_match28', [r'Twenty-eight, incumbent three & (\d+)/11'], 0),
    ('bb_match35', [r'Thirty-five, incumbent three & \d+/11 & (\d+)/11'], 0),
    ('bb_match35lit', [r'Thirty-five, plus literature 8 & \d+/11 & (\d+)/11'], 0),
    ('bb_match28', [r'Twenty-eight, incumbent three & \d+/11 & (\d+)/11'], 0),
    ('bb_meandiff28', [r'there falls from 0\.025 to ([\d.]+)'], 0.0005),
    ('bbpre_28', [r'which stays at ([\d.]+) on the twenty-eight-compound panel'],
     0.0005),
    # 네 번째 구성 (2026-10-20). 표 1 네 번째 블록 + §3.1 산문 + Limitations.
    # 표 3 사건의 교훈대로 표를 먼저 등록한다.
    ('pp_st_md', [row_pat('35 cpd., incumbent 3', 'structural', 0)], 0.0005),
    ('pp_st_lr', [row_pat('35 cpd., incumbent 3', 'structural', 2)], 0.0005),
    ('pp_st_ds', [row_pat('35 cpd., incumbent 3', 'structural', 3)], 0.0005),
    ('pp_el_md', [row_pat('\\quad + literature 5', 'electronic', 0)], 0.0005),
    ('pp_el_mdn', [row_pat('\\quad + literature 5', 'electronic', 1)], 0.0005),
    ('pp_el_nv', [row_pat('\\quad + literature 5', 'electronic', 5)], 0.0005),
    ('pp_mg_md', [row_pat('\\quad + pesticides 2', 'merged', 0)], 0.0005),
    ('pp_mg_lr', [row_pat('\\quad + pesticides 2', 'merged', 2)], 0.0005),
    ('pp_mg_ds', [row_pat('\\quad + pesticides 2', 'merged', 3)], 0.0005),
    ('pp_n_perm', [r'\\binom\{24\}\{10\}=([\d\\,]+)'], 0),
    ('pp_min_p', [r'smallest nominal value of ([\d.]+) and a'], 0.00001),
    ('pf_meff', [r'\$m_\{\\mathrm\{eff\}\}\$ of 18 \(Li-Ji\)     & [\d.]+  & [\d.]+  & [\d.]+  & ([\d.]+)'], 0.005),
    ('pf_prespec', [r'Statistic specified in advance \(3 tests\) & [\d.]+ & [\d.]+ & '
                    r'\\textbf\{[\d.]+\} & \\textbf\{([\d.]+)\}'], 0.0005),
    ('struct_removed', [r'remove (\d+) of the 51 on the primary panel'], 0),
    ('struct_removed_28', [r'of the 51 on the primary panel and\s*(\d+) on the'], 0),
    ('struct_kept_35', [r'leaving (\d+) and \d+ structural coordinates'], 0),
    ('struct_kept_28', [r'leaving \d+ and (\d+) structural coordinates'], 0),
    ('struct_candidates', [r'removes \d+ of the (\d+)',
                           r'descriptors \((\d+) after excluding run metadata\)'], 0),
    ('pp_bh_q', [r'Benjamini-Hochberg value of ([\d.]+) over the same',
                 r'produces an adjusted\s*\n?(?:value of )?(\d+\.\d+)',
                 r'into one with an adjusted value of\s*\n?(\d+\.\d+)'], 0.00005),
    # 인용 문헌의 수치 (2026-10-20). outputs/cite_facts.json 이 단일 출처다.
    ('j18_identified', [r'identified\s*\n?(\d+) simulants'], 0),
    ('j18_unsuitable', [r'judged (\d+) unsuitable'], 0),
    ('j21_identified', [r'a later one screened (\d+) chemicals'], 0),
    ('j21_discounted', [r'discounted (\d+) on the risk'], 0),
    ('j21_suitable', [r'repository of (\d+) suitable'], 0),
    ('j21_possible', [r'suitable and (\d+) provisional'], 0),
]

WORDS = {'zero': 0, 'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6,
         'seven': 7, 'eight': 8, 'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12,
         'thirteen': 13, 'fourteen': 14, 'fifteen': 15, 'sixteen': 16,
         'seventeen': 17, 'eighteen': 18, 'nineteen': 19, 'twenty': 20,
         'twenty-one': 21, 'twenty-two': 22,
         'twenty-six': 26, 'twenty-seven': 27, 'twenty-eight': 28,
         'twenty-nine': 29, 'thirty': 30,
         'twenty-three': 23, 'twenty-four': 24, 'twenty-five': 25,
         'thirty-three': 33, 'thirty-four': 34}


# 폐기값 블랙리스트. 재분석으로 값이 바뀔 때마다 구판을 여기에 넣는다.
# 본문이 구판임을 밝히고 인용하는 문장은 면제 패턴으로 통과시킨다.
SUPERSEDED = [
    ('28종 구조공간 명목 p (정본 0.031)', [r'0\.028\b'], [r'earlier version', r'superseded']),
    ('부분공간 백분위 (정본 25.4)',
     [r'25\.8th percentile', r'\b26th percentile', r'\b28th percentile',
      r'22\.7th percentile', r'\b23rd percentile'], []),
    ('도판2 캡션의 통계량 서술 (실제로는 평균차 하나)',
     [r'Smallest \$p\$ over\s+the six statistics for 2000'], []),
    ('부분공간 0.05 미만 비율 (정본 32.6)',
     [r'36\.3\\?%', r'36\\% of them fall below', r'17\\% after whitening'], []),
    ('사전지정 BH q (정본 0.046)', [r'0\.042\b'], []),
    ('28종 BH over 18 (정본 0.185)',
     [r'adjusted value is 0\.159', r'eighteen tests is 0\.159',
      r'over 18\s+& [\d.]+ & [\d.]+ & 0\.159'], []),
    ('28종 Westfall-Young (정본 0.079)', [r'approximating it, gives 0\.075'], []),
    ('28종 매치 수 (정본 4)',
     [r'restores three matches', r'counted three matches'],
     [r'earlier version', r'superseded', r'descriptor set']),
    ('28종 구조공간 유의 통계량 수 (정본 1)', [r'two of six statistics in the structural'], []),
    ('28종 구조좌표 수 (정본 36)',
     [r'38-dimensional value', r'thirty-eight structural coordinates are not'], []),
    ('tabun 정규화거리 (정본 1.989)', [r'2\.003'],
     [r'earlier', r'superseded', r'put tabun at']),
    ('두 패널 감소가 같다는 서술 (정본: 13 대 15)',
     [r'attrition is the same on\s+both panels'], []),
]


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
    # 공백 정규화본. 문서와 도판 주석 레코드 사이에는 \x00 을 넣는다. 패턴의
    # \s 는 이 문자를 먹지 않으므로 경계를 넘는 오탐이 생기지 않는다. 넣지
    # 않았더니 주석 두 레코드에 걸친 문자열이 한 칸으로 읽혀 세 건이 잘못
    # 걸렸다 (2026-10-24).
    _segs = [re.sub(r'\s+', ' ', seg) for seg in
             ([tex] + [strip_comments(open(x).read()) for x in
                       ('paper/cover_letter.tex', 'paper/highlights.tex',
                        'paper/supplementary.tex') if os.path.exists(x)]
              )]
    # 도판 주석은 레코드가 이미 한 줄이라 줄바꿈 문제가 없다. 정규화본에
    # 넣으면 레코드끼리 붙어 오탐만 생기므로 원문 쪽에만 남긴다.
    body_flat = '\x00'.join(_segs)
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
        # 2026-10-24. 패턴이 리터럴 공백 하나를 요구하는데 그 자리에서 TeX 줄이
        # 바뀌면 문맥을 놓친다. 처음에는 '못 찾았을 때만' 공백을 눌러 다시
        # 보게 했는데, 그러면 같은 수치가 두 군데 적혀 있고 한쪽만 줄바꿈에
        # 걸린 경우 그 한쪽이 검사되지 않는다. 실제로 본문의 백분위 두 개가
        # 그 경로로 구판인 채 통과했다. 항상 둘 다 훑고 중복만 걷어낸다.
        for pat in pats:
            for m in re.finditer(pat, body_flat):
                found.append((tonum(m.group(1)), m.group(0)[:62]))
        seen, uniq = set(), []
        for val, ctx in found:
            k = (round(val, 12), ' '.join(ctx.split()))
            if k in seen:
                continue
            seen.add(k)
            uniq.append((val, ctx))
        found = uniq
        if not found:
            unseen.append(f'{key}: 본문에서 문맥을 찾지 못함 (값 {want:g})')
            continue
        for got, ctx in found:
            checked += 1
            if abs(got - want) > tol:
                fails.append(f'  불일치 {key}: 본문 {got:g} vs 산출물 {want:g}'
                             f'   [{ctx}]')

    # (B2) 폐기값 블랙리스트 (2026-10-25 신설).
    #
    # 왜 필요한가. (B) 의 검사는 어법 단위다. 'from $p=0.031$ to' 라는 패턴
    # 하나가 어딘가에서 맞으면 통과하므로, 같은 양을 다른 문장으로 적은
    # 자리에 구판 값이 남아 있어도 보이지 않는다. 2026-10-25 에 그 경로로
    # 본문 수치 16건이 PASS 를 받은 채 남아 있었다. 여기서는 방향을 뒤집어,
    # 폐기된 값이 원고 어디에든 나타나면 실패로 본다.
    #
    # 각 항목은 (사유, 금지 패턴들, 면제 패턴들)이다. 면제는 구판임을
    # 명시하고 인용하는 문장을 위한 것이다. 금지 패턴이 맞은 자리 앞뒤
    # 240자 안에서 면제 패턴이 맞으면 의도적 인용으로 본다.
    superseded = []
    for why, pats, exempt in SUPERSEDED:
        hits = []
        for src in (body, body_flat):
            for pat in pats:
                for m in re.finditer(pat, src):
                    a, b = max(0, m.start() - 240), m.end() + 240
                    win = src[a:b]
                    if any(re.search(e, win, re.I) for e in exempt):
                        continue
                    hits.append(' '.join(m.group(0).split())[:64])
        if hits:
            uniq = sorted(set(hits))
            superseded.append(f'  폐기값 잔존 {why}: ' + ' / '.join(uniq[:3]))
    if superseded:
        fails.append('폐기값이 원고에 남아 있다 (검사 B2):')
        fails.extend(superseded)

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
    # 2026-10-24. 종전에는 '문맥 미발견' 이 경고였고 판정은 PASS 였다. 그래서
    # 본문을 다시 쓰거나 표를 다시 생성해 정규식이 안 맞게 되면 검사가 조용히
    # 꺼지고 화면에는 통과가 찍혔다. 이 세션에서만 그 경로로 세 번 빠뜨렸다.
    # 꺼진 검사는 없는 검사보다 나쁘다. 실패로 올린다. 본문에 싣지 않기로
    # 한 사실은 UNQUOTED 에 이유와 함께 적어 명시적으로 면제한다.
    unexpected = [u for u in unseen
                  if u.split(':')[0] not in UNQUOTED and '도판' not in u]
    if unseen:
        print(f'[verify] 문맥 미발견 {len(unseen)}건 (검사 누락 — 정규식 보강 필요)')
        for u in unseen:
            print('  ', u)
    if unexpected:
        fails.append(f'  검사 누락 {len(unexpected)}건: '
                     + ', '.join(u.split(':')[0] for u in unexpected)
                     + ' — 정규식을 고치거나 UNQUOTED 에 면제 사유를 적을 것')
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
