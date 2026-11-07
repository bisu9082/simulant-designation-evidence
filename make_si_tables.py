#!/usr/bin/env python3
r"""SI 의 손작성 표 두 개를 산출물에서 생성한다 (2026-11-03 신설).

왜 필요한가. GATE 8 8차의 리뷰어 1 이 SI 의 tab:ablation 차원 열이 40/42 로
구판이라고 지적했고, 사실이었다. 그 표는 10열 24수치를 손으로 유지하고
있었고 어떤 스크립트도 만들지 않았다. 더 나쁜 것은 같은 라벨
tab:ablation 이 make_step4_outputs.py 가 만드는 6열 표에도 붙어 있어
이름이 충돌한다는 점이다.

검증기의 G-D 가드는 '원고가 인용하는 예치 파일에 생성기가 있는가' 를 보지만
표는 보지 않는다. 그래서 이 표는 T2-1(상관필터 교체)과 T2-8(형태이성질체
재생성)을 전부 통과해 살아남았다. 손작성 표가 남아 있으면 같은 일이 또 난다.

여기서 만드는 것:
  outputs/tab_ablation.tex      패널 4구성 x 공간 3종 x 통계량 6종 + 최소 q
  outputs/tab_blockbalance.tex  블록 균등 가중 대조

본문·SI 는 \input 으로 읽는다. 수치를 손으로 적는 자리를 없애는 것이 목적이다.

사용: python3 make_si_tables.py
"""
import json
import os

STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']
SPACE_LABEL = {'structural': 'structural', 'electronic3': 'electronic',
               'merged': 'merged'}

# 표의 행 순서. (SI 가 쓰는 표시 이름, rollback_c 키 또는 'pesticide')
# 2026-11-03 (2차). 종전에는 블록마다 이름 두 개만 두고 셋째 행(병합공간)의
# 이름 칸을 비웠다. 농약 구성은 이름이 세 줄이어서 둘째 칸에 \\ 로 두 줄을
# 넣었고, 그러면 '\quad + literature 5\\' 가 독립된 한 행이 되어 이름이
# 공간과 한 칸씩 어긋나 찍혔다(전자공간 옆에 '+ pesticides 2'). 손표는 세
# 이름을 세 공간에 하나씩 붙였다. 칸을 셋으로 늘려 그 배치를 되살린다.
ROWS = [
    ('\\textbf{35 cpd., incumbent 3}', '\\textbf{(primary)}', '',
     'panel35_designated3'),
    ('35 cpd., eight designated', '\\quad (3 incumbent + 5 literature)', '',
     'panel35_designated_plus_literature8'),
    ('28 cpd., incumbent 3', '\\quad (omits the seven)', '',
     'panel28_designated3'),
    ('35 cpd., ten designated', '\\quad (3 incumbent + 5 literature',
     '\\quad + 2 pesticides)', 'pesticide'),
]


# 2026-11-04 (작업 D). 패널 명명 규칙. 종전에는 같은 구성을 'literature 8',
# 'literature eight', 'the eight literature simulants' 로 적었다. 8 은 문헌이
# 제안한 화합물 수가 아니라 **지정 후보의 총수**다: 기존 3종(DMMP, DPGME,
# 디클로로메탄) + 문헌 5종(CEES, DCNP, DCP, DFP, 살리실산메틸). 아래 라벨은
# 문헌 쪽 수(5)만 적고, 총수 8 은 본문이 문장으로 말한다.
#   28종 구성 = 35종에서 7종을 뺀 것 (문헌 5 + 농약 2)

def fmt(p, dp=4):
    """p 를 SI 관례로 쓴다. 0.001 미만은 <0.001, 0.05 미만은 굵게."""
    if p is None:
        return '--'
    if p < 0.001:
        return '\\textbf{$<$0.001}'
    s = f'{p:.{dp}g}'
    return f'\\textbf{{{s}}}' if p < 0.05 else s


def fmt3(p):
    """소수 3자리 고정. tab:family 가 쓰는 서식."""
    if p is None:
        return '--'
    if p < 0.001:
        return '\\textbf{$<$0.001}'
    s = f'{p:.3f}'
    return f'\\textbf{{{s}}}' if p < 0.05 else s


def load():
    rc = json.load(open('outputs/rollback_c.json'))['panels']
    mq = json.load(open('outputs/multiplicity_across_panels.json'))
    pf = json.load(open('outputs/pesticide_family.json'))
    pp = json.load(open('outputs/pesticide_panel.json'))['spaces']
    # 2026-11-03 (2차). tab:ablation 의 마지막 열은 캡션이 '해당 패널
    # 구성의 18검정에 대한 BH 최소 q' 라고 적는다. 1차에서 내가
    # multiplicity_across_panels(72검정)를 붙여 캡션과 어긋났고, 본문
    # 571/809 행이 0.6299(18검정)를 인용하는 동안 표는 0.3539(72검정)를
    # 찍고 있었다. 캡션이 말하는 양을 쓴다.
    mc = {r[0]: r[4] for r in json.load(open('outputs/multiplicity_c.json'))}
    mc['pesticide10'] = pf['bh_all_min_q']
    return rc, mq, pf, pp, mc


def ablation(rc, mq, pf, pp, mc):
    MQ = mc
    lines = []
    for i, (a, b, c, key) in enumerate(ROWS):
        if key == 'pesticide':
            dims = {s: pp[s]['n_dims'] for s in pp}
            pvals = {}
            for s in SPACE_LABEL:
                src = 'electronic' if s == 'electronic3' else s
                pvals[s] = {st: pf['per_test_p'].get(f'{src}:{st}')
                            for st in STATS}
            q = MQ.get('pesticide10')
        else:
            dims = {s: rc[key][s]['n_descriptors'] for s in rc[key]}
            pvals = {s: {st: rc[key][s]['statistics'][st]['p'] for st in STATS}
                     for s in rc[key]}
            q = MQ.get(key)
        qmin = q['q'] if isinstance(q, dict) else q
        for j, s in enumerate(['structural', 'electronic3', 'merged']):
            left = (a, b, c)[j]
            cells = ' & '.join(fmt(pvals[s][st]) for st in STATS)
            qcell = fmt(qmin) if s == 'merged' and qmin is not None else ''
            lines.append(f"{left:<22} & {SPACE_LABEL[s]:<11} & "
                         f"{dims[s]:>2} & "
                         f"{cells} & {qcell} \\\\")
        lines.append('\\hline')
    return '\n'.join(lines)


# ---------------------------------------------------------------- tab:family
FAM_COLS = [('panel35_incumbent3', '35 cpd.,', 'incumbent 3', '(primary)'),
            ('panel35_incumbent_plus_lit8', '35 cpd.,',
             'eight designated', ''),
            ('panel28_incumbent3', '28 cpd.,', 'incumbent 3',
             '(omits the seven)'),
            ('pesticide', '35 cpd.,', 'ten designated', '(pesticides in)')]
FAM_ROWS = [('Smallest nominal $p$ (18 tests)', 'min_p'),
            ('Benjamini-Hochberg over 18', 'bh_all_min_q'),
            ('Benjamini-Hochberg over 6, per space', 'bh_by_space_min_q'),
            ('Westfall-Young min-$P$ over 18', '_wy'),
            ('$m_{\\mathrm{eff}}$ of 18 (Li-Ji)', 'meff'),
            ('Bonferroni at $m_{\\mathrm{eff}}$', 'meff_adjusted_min_p'),
            ('Statistic specified in advance (3 tests)',
             'prespecified_bh_min_q'),
            ('Benjamini-Hochberg over all 72', '_ap')]


def family(mq):
    ss = json.load(open('outputs/stat_supplement.json'))['family']
    pf = json.load(open('outputs/pesticide_family.json'))
    WY = json.load(open('outputs/stat_supplement2.json'))['westfall_young']
    AP = mq['min_q_by_config']
    APKEY = {'panel35_incumbent3': 'panel35_designated3',
             'panel35_incumbent_plus_lit8':
                 'panel35_designated_plus_literature8',
             'panel28_incumbent3': 'panel28_designated3',
             'pesticide': 'pesticide10'}
    lines = []
    for label, key in FAM_ROWS:
        cells = []
        for cfg, _a, _b, _c in FAM_COLS:
            if key == '_ap':
                v = AP[APKEY[cfg]]['q']
            elif key == '_wy':
                # Westfall-Young 은 stat_supplement2 (패널 3구성) 와
                # pesticide_family (농약 구성) 에 나뉘어 있다.
                v = (pf['westfall_young_adjusted_p'] if cfg == 'pesticide'
                     else WY.get(cfg, {}).get('wy_adjusted_p'))
            else:
                src = pf if cfg == 'pesticide' else ss[cfg]
                v = src.get(key)
            if key == 'meff':
                cells.append(f'{v:.2f} ' if v is not None else '--')
            else:
                cells.append(fmt3(v))
        lines.append(f'{label:<40} & ' + ' & '.join(cells) + ' \\\\')
    head = (r"""\begin{tabular}{lcccc}
\hline
Accounting scheme & """ + ' & '.join(c[1] for c in FAM_COLS) + r""" \\
                  & """ + ' & '.join(c[2] for c in FAM_COLS) + r""" \\
                  & """ + ' & '.join(c[3] for c in FAM_COLS) + r""" \\
\hline
""")
    return head + '\n'.join(lines) + '\n\\hline\n\\end{tabular}'


# ------------------------------------------------------ tab:blockbalance
BB_ROWS = [('Thirty-five, incumbent three', 'panel35_incumbent3'),
           ('Thirty-five, incumbent 3 + literature 5',
            'panel35_incumbent_plus_lit8'),
           ('Twenty-eight, incumbent three', 'panel28_incumbent3')]


S6 = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
      'n_closer', 'n_valid']


def _bh_min(pvals):
    """같은 6검정 가족에 대한 BH 최소 q. 두 가중 방식을 비교 가능하게 한다.

    2026-11-03. 구판 SI 는 '한 쪽은 18검정 BH, 다른 쪽은 미상' 으로 두
    열을 나란히 놓고 있었고 그 q 는 어느 산출물에도 없었다. 양쪽을 같은
    6개 병합공간 통계량 가족으로 계산하고 캡션에 그렇게 적는다.
    """
    import numpy as _np
    o = _np.sort(_np.array(sorted(pvals)))
    m = len(o)
    bh = _np.minimum.accumulate((o * m / _np.arange(1, m + 1))[::-1])[::-1]
    return float(bh[0])


def blockbalance():
    """두 가중 방식 비교. 마지막 두 열은 18검정 가족의 BH 최소 q 다.

    2026-11-03 (2차). 1차에서는 이 두 열을 '병합공간 6검정만의 BH' 로
    계산했는데, verify_consistency 의 정본 사실 bbq_* 는 구조 6 + 전자 6 +
    병합 6 = 18검정 가족의 BH 다. 같은 이름의 서로 다른 양을 만들어 낸
    것이므로 정본 쪽에 맞춘다. 두 열의 차이는 병합공간 6검정을 어느
    가중으로 계산했는가 하나뿐이다.
    """
    bb = json.load(open('outputs/blockbalance.json'))
    rc = json.load(open('outputs/rollback_c.json'))
    RC = {'panel35_incumbent3': 'panel35_designated3',
          'panel35_incumbent_plus_lit8':
              'panel35_designated_plus_literature8',
          'panel28_incumbent3': 'panel28_designated3'}
    lines = []
    for label, key in BB_ROWS:
        d = bb[key]
        e, b = d['equal_per_descriptor'], d['block_balanced']
        sp = rc['panels'][RC[key]]
        fixed = [sp[sx]['statistics'][k]['p']
                 for sx in ('structural', 'electronic3') for k in S6]
        qe = _bh_min(fixed + [e['p'][k] for k in S6])
        qb = _bh_min(fixed + [b['p'][k] for k in S6])
        lines.append(f"{label} & {e['n_valid']}/11 & {b['n_valid']}/11 & "
                     f"{min(e['p'][k] for k in S6):.3f} & "
                     f"{min(b['p'][k] for k in S6):.3f}"
                     f" & {qe:.3f} & {qb:.3f} \\\\")
    head = (r'\begin{tabular}{lcccccc}' '\n' r'\hline' '\n'
            r' & \multicolumn{2}{c}{Valid matches}'
            r' & \multicolumn{2}{c}{Smallest nominal $p$}' '\n'
            r' & \multicolumn{2}{c}{Smallest $q$ over 18} \\' '\n'
            r'Configuration & per desc. & per block & per desc. & per block'
            '\n' r' & per desc. & per block \\' '\n' r'\hline' '\n')
    return head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + r'\end{tabular}'


# ----------------------------------------------------------- tab:sparse
SP_LABEL = {'baseline': 'As reported',
            'drop_binary_nonzero_le_1': 'Drop single-carrier binaries',
            'drop_binary_nonzero_le_2': 'Drop binaries with $\\le 2$',
            'drop_binary_nonzero_le_3': 'Drop binaries with $\\le 3$',
            'drop_binary_nonzero_le_4': 'Drop binaries with $\\le 4$',
            'drop_binary_nonzero_le_5': 'Drop binaries with $\\le 5$',
            'drop_cardinality_le_5': 'Drop cardinality $\\le 5$'}


def sparse():
    sp = json.load(open('outputs/sparse_indicators.json'))
    lines = []
    for v in sp['variants']:
        lab = SP_LABEL.get(v['variant'], v['variant'])
        lines.append(f"    {lab:<28} & {v['n_dropped']:>2} & {v['n_dims']} & "
                     f"{v['n_match']} & {v['tabun_d_norm']:.4f} & "
                     f"{v['p_one_sided']:.4f} \\\\")
    head = (r'  \begin{tabular}{lrrrrr}' '\n' r'    \hline' '\n'
            r'    Encoding & Removed & Dim. & Matches'
            r' & Tabun $d_{\mathrm{norm}}$ & $p$ \\' '\n'
            r'    \hline' '\n')
    return head + '\n'.join(lines) + '\n' + r'    \hline' + '\n' + r'  \end{tabular}' 


# ----------------------------------------------------------- tab:plocal
# 2026-11-03 신설. 이 표의 44개 소수는 손으로 유지되고 있었고 어떤
# 스크립트도 만들지 않았다. 여섯 화합물 x 네 열은 모두
# outputs/_electronic_with_gap.csv 에 있고, 각주의 유도 통계도 거기서
# 계산된다. omega_k = omega * f^+(P) 는 Domingo 2016 의 정의다.
PL_ROWS = [('GB', 'sarin (GB)'), ('GD', 'soman (GD)'),
           ('GF', 'cyclosarin (GF)'), (r'\textbf{DMMP}', 'DMMP'),
           ('GA', 'tabun (GA)'), ('VX', 'VX')]


def plocal():
    import csv as _csv
    with open('outputs/_electronic_with_gap.csv') as fh:
        d = {r['name']: r for r in _csv.DictReader(fh)}
    g = lambda n, c: float(d[n][c])
    wk = lambda n: g(n, 'electrophilicity_index_eV') * g(n, 'fukui_f_plus_P')
    lines = []
    for lab, name in PL_ROWS:
        lines.append(f"{lab} & {g(name,'q_P'):.5f} & "
                     f"{g(name,'q_phosphoryl_O'):.5f} & "
                     f"{g(name,'wiberg_PO'):.4f} & "
                     f"{g(name,'fukui_f_plus_P'):.5f} & "
                     f"{g(name,'electrophilicity_index_eV'):.2f} & "
                     f"{wk(name):.2f} \\\\")
    G = ['sarin (GB)', 'soman (GD)', 'cyclosarin (GF)']
    import statistics as _st
    fm, fs = _st.mean(g(n, 'fukui_f_plus_P') for n in G), \
        _st.stdev([g(n, 'fukui_f_plus_P') for n in G])
    km, ks = _st.mean(wk(n) for n in G), _st.stdev([wk(n) for n in G])
    zf = (g('DMMP', 'fukui_f_plus_P') - fm) / fs
    zk = (wk('DMMP') - km) / ks
    rf = g('DMMP', 'fukui_f_plus_P') / g('VX', 'fukui_f_plus_P')
    rk = wk('DMMP') / wk('VX')
    foot = (
        r'\multicolumn{7}{l}{\footnotesize GB/GD/GF mean '
        r'$f^{+}(\mathrm{P})$ = ' + f'{fm:.4f}' + '\n'
        r'$\pm$ ' + f'{fs:.4f}' + ', DMMP ' + f'{abs(zf):.2f}'
        + r' SD below; mean $\omega_k$ = ' + f'{km:.2f}'
        + r' $\pm$ ' + f'{ks:.2f}' + r',}\\' + '\n'
        r'\multicolumn{7}{l}{\footnotesize DMMP ' + f'{abs(zk):.2f}'
        + r' SD below. $\omega$ and $\omega_k$ are' + '\n'
        r'reported here and used in no distance '
        r'(Section~\ref{M-meth:spaces}).}\\' + '\n'
        r'\multicolumn{7}{l}{\footnotesize DMMP/VX = ' + f'{rf:.2f}'
        + r' by $f^{+}(\mathrm{P})$ but ' + f'{rk:.2f}' + r' by' + '\n'
        r'$\omega_k$. Values are GFN2-derived and relative within one '
        r'protocol.}\\')
    head = (r'\begin{tabular}{lcccccc}' '\n' r'\hline' '\n'
            r'Compound & $q_{\mathrm{P}}$ (e) & $q_{\mathrm{P=O}}$ (e) '
            r'& $W_{\mathrm{P=O}}$ &' '\n'
            r'$f^{+}(\mathrm{P})$ & $\omega$ (eV) & $\omega_k$ (eV) \\'
            '\n' r'\hline' '\n')
    return (head + '\n'.join(lines) + '\n' + r'\hline' + '\n'
            + foot + '\n' + r'\end{tabular}')


# ------------------------------------------------ tab:volatility, tab:vpairing
# 2026-11-03 신설. 휘발성 보조 분석의 두 표. 손으로 유지하면 다른 네 표와
# 같은 전파 실패를 되풀이하므로 처음부터 생성으로 만든다. 출처 문자열은
# volatility_space.VP 에 들어 있는 조회 원문이고 여기서 고치지 않는다.
VOL_ORDER_DESC = True       # 증기압 내림차순


def _vol():
    return json.load(open('outputs/volatility_space.json'))


def volatility():
    """채택 증기압 표. 값·온도·tier·출처."""
    d = _vol()
    vp = d['vapour_pressure_mmHg']
    import csv as _csv
    with open('descriptor_cache_xtb.csv') as fh:
        names = [r['compound_name'] for r in _csv.DictReader(fh)]
    rows = sorted(((vp[n]['value'], n) for n in names), reverse=VOL_ORDER_DESC)
    lines = []
    for val, n in rows:
        e = vp[n]
        t = '' if e['temperature_C'] == 25 else f"{e['temperature_C']}"
        src = e['source'].split('|')[0].strip().rstrip(',')
        if len(src) > 46:
            src = src[:44].rstrip() + '.'
        lines.append(f"{n} & {val:.4g} & {t} & {e['tier']} & {src} \\\\")
    head = (r'\begin{tabular}{lrccl}' '\n' r'\hline' '\n'
            r'Compound & $P$ (mmHg) & $T$ ($^{\circ}$C) & Tier'
            r' & Source \\' '\n' r'\hline' '\n')
    foot = (r'\multicolumn{5}{l}{\footnotesize Blank temperature means '
            r'25~$^{\circ}$C. Tier 1 is a vapour-pressure measurement or}\\'
            '\n'
            r'\multicolumn{5}{l}{\footnotesize compilation at '
            r'25~$^{\circ}$C, tier 2 a government compendium at '
            r'25~$^{\circ}$C, tier 3 the}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize nearest available '
            r'temperature. Source strings are as returned by PubChem; '
            r'the}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize full strings and every '
            r'candidate value are in the deposited JSON. Tier X is a '
            r'compilation}\\' '\n'
            # 2026-11-06 (9차 P-10). TX 가 각주에 정의돼 있지 않았다. 그리고
            # 축 자체에 불확도가 없어서 0.473 대 0.216 의 0.257 여유를
            # 온도차가 만들 수 있는지 독자가 판단할 수 없었다.
            r'\multicolumn{5}{l}{\footnotesize value without a stated '
            r'measurement temperature, used only for solvents. Over this '
            r'compound}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize class the '
            r'Clausius--Clapeyron slope is of order 0.03 to 0.05 '
            r'decades per $^{\circ}$C, so a}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize 5~$^{\circ}$C mismatch '
            r'moves a value by roughly 0.15 to 0.25 decades. No value was '
            r'converted}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize between temperatures.}\\')
    return (head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + foot
            + '\n' + r'\end{tabular}')


def vpairing():
    """작용제별 최근접 모사체 / 최근접 풀 화합물. 기술표."""
    d = _vol()['pairing_table']
    lines = []
    for r in d['rows']:
        ser = r['series'].replace('vomit/irritant', 'irritant')
        lines.append(f"{r['agent']} & {ser} & {r['log10_vp']:.3f} & "
                     f"{r['nearest_simulant']} & "
                     f"{r['gap_to_nearest_simulant']:.3f} & "
                     f"{r['nearest_pool_compound']} & "
                     f"{r['gap_to_nearest_pool']:.3f} \\\\")
    head = (r'\begin{tabular}{llrlrlr}' '\n' r'\hline' '\n'
            r'Agent & Series & $\log_{10}P$ & Nearest simulant & $\Delta$'
            r' & Nearest of pool & $\Delta$ \\' '\n' r'\hline' '\n')
    foot = (r'\multicolumn{7}{l}{\footnotesize $\Delta$ is the gap in '
            r'$\log_{10}(P/\mathrm{mmHg})$. The pool is the seventeen '
            r'non-agents, so}\\' '\n'
            r'\multicolumn{7}{l}{\footnotesize the last column is a '
            r'simulant where one is nearest and an interferent '
            r'otherwise.}\\')
    return (head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + foot
            + '\n' + r'\end{tabular}')


def premisenull():
    """전제검정 라벨순열 귀무. outputs/premise_null.json 에서 생성."""
    d = json.load(open('outputs/premise_null.json'))
    SET = {'d1': 'First', 'd2': 'Second', 'd3': 'Third'}
    CH = {'uv': 'ultraviolet', 'normal': 'reflectance'}
    lines = []
    for c in d['cells']:
        src = (f"{c['n_null']:,} enumerated" if c['enumerated']
               else f"{c['n_null']:,} sampled")
        lines.append(f"{SET[c['dataset']]} & {CH[c['mode']]} & "
                     f"${c['observed']:+.4f}$ & ${c['null_mean']:+.4f}$ & "
                     f"{c['null_sd']:.4f} & {c['p_one_sided']:.4f} & "
                     f"{src} \\\\")
    head = (r'\begin{tabular}{llrrrrl}' '\n' r'\hline' '\n'
            r'Set & Channel & Observed & Null mean & Null s.d. & $p$'
            r' & Null built from \\' '\n' r'\hline' '\n')
    foot = (r'\multicolumn{7}{l}{\footnotesize Adjusted values are '
            r'Benjamini--Hochberg over the six cells; the smallest is }\\'
            '\n'
            r'\multicolumn{7}{l}{\footnotesize '
            + f"{d['bh']['min_q']:.4f}" + r' and all six clear $q=0.05$.}\\')
    return (head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + foot
            + '\n' + r'\end{tabular}')


def registers():
    """세 인컴번트의 규제 레지스터 보유 현황 (2026-11-06 신설).

    GATE 8 9차에서 본문이 'DMMP 는 어떤 종류의 규제 흡입값도 갖지 않는다'
    고 적은 것이 DOE PAC 과 NIOSH IDLH 로 반증됐다. 공백은 규제값의 공백이
    아니라 EPA AEGL 한 레지스터의 공백이다. 수집분은 sources 에 출처 URL
    과 함께 예치했고 이 표는 그 파일에서 생성한다.
    """
    import csv
    SRC = 'sources/registers_20261106.csv'
    KEEP = ('DMMP', 'DPGME', 'dichloromethane')
    rows = []
    with open(SRC, newline='') as fh:
        for r in csv.DictReader(fh):
            if r['compound_name'] not in KEEP:
                continue
            rows.append(r)
    order = {n: i for i, n in enumerate(KEEP)}
    rows.sort(key=lambda r: (order[r['compound_name']], r['register'], r['tier']))
    lines = []
    for r in rows:
        reg = r['register'].replace('_', ' ')
        lines.append(f"{r['compound_name']} & {reg} & {r['tier']} & "
                     f"{r['value']} & {r['unit'].replace('mg/m3', 'mg m$^{-3}$')}"
                     f" & {r['basis'] or '--'} \\\\")
    head = (r'\begin{tabular}{lllrll}' '\n' r'\hline' '\n'
            r'Compound & Register & Tier & Value & Unit & Basis \\' '\n'
            r'\hline' '\n')
    return head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + r'\end{tabular}'


def designations():
    """지정 기록 6열 표 (2026-11-06 신설).

    9차에서 지적받은 자리다. 지정 귀속이 산문 두 단락에만 있었고, 그래서
    'NRC 가 VX 에 지정' 을 'the nerve agent class' 로 요약하는 오기와,
    'designating documents assign these compounds no agent class' 라는
    거짓 서술이 동시에 살아남았다. 문서/작용제/측정/선정근거/출처를
    한 줄에 놓으면 그런 요약이 불가능해진다.
    """
    import csv
    rows = []
    with open('sources/designations_20261106.csv', newline='') as fh:
        for r in csv.DictReader(fh):
            rows.append(r)

    def esc(t):
        return (t.replace('&', r'\&').replace('%', r'\%')
                 .replace('[수집 불가', '[not obtained')
                 .replace('특허 본문 미열람 — Google Patents 및 Justia 접근 거부',
                          'patent text not read; Google Patents and Justia refused')
                 .replace('공개 기록에서 지정 문서를 찾지 못함',
                          'no designating document found in the public record'))

    lines = []
    for r in rows:
        lines.append(' & '.join(esc(r[k]) for k in
                                ('compound', 'document', 'agent',
                                 'measurement', 'stated_basis')) + r' \\')
    head = (r'\begin{tabular}{p{1.9cm}p{3.4cm}p{1.5cm}p{2.3cm}p{4.2cm}}' '\n'
            r'\hline' '\n'
            r'Compound & Document & Agent & Measurement & Basis the document '
            r'states \\' '\n' r'\hline' '\n')
    return head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + r'\end{tabular}'


def dyepanel():
    """29종 어레이 구성 (2026-11-06 신설).

    9차 지적. 본문이 '같은 29종을 같은 순서로' 라고만 적고 어느 화합물인지
    어디에도 없었다. 선행 세 편의 Table 1 을 직독해 목록을 예치하고,
    (a) 염료 15 의 용매가 세 번째 런에서 다르게 적혀 있다는 것과
    (b) 여섯 항목이 발색단이 아니라 형광체 또는 유도체화 시약이라는 것을
    함께 공시한다. 종전에는 둘 다 공시되지 않았다.
    """
    import csv
    rows = list(csv.DictReader(open('sources/dye_panel_20261106.csv')))
    half = (len(rows) + 1) // 2
    left, right = rows[:half], rows[half:]
    lines = []
    for i in range(half):
        a = left[i]
        b = right[i] if i < len(right) else {'index': '', 'name': ''}
        # 2026-11-07 (10차). 종전에는 note 가 있는 행 전부에 $^{a}$ 를
        # 달았다. 염료 15 의 note 는 용매 불일치인데 캡션이 그 마커를
        # '비염료' 로 설명하고 있어 표가 캡션을 거짓으로 만들었다.
        # note_kind 로 비염료만 a, 그 밖은 b 를 단다.
        _m = {'nondye': '$^{a}$', 'other': '$^{b}$'}
        mark = _m.get(a.get('note_kind', ''), '')
        markb = _m.get(b.get('note_kind', ''), '')
        lines.append(f"{a['index']} & {a['name']}{mark} & "
                     f"{b['index']} & {b['name']}{markb} \\\\")
    head = (r'\begin{tabular}{rlrl}' '\n' r'\hline' '\n'
            r'\# & Compound & \# & Compound \\' '\n' r'\hline' '\n')
    notes = [(('$^{a}$ ' if r['note_kind'] == 'nondye' else '$^{b}$ ')
              + f"{r['index']} {r['name']}: {r['note']}")
             for r in rows if r['note']]
    foot = '\n'.join(
        r'\multicolumn{4}{p{13cm}}{\footnotesize ' + n + r'}\\'
        if True else
        r'\multicolumn{4}{p{13cm}}{\footnotesize ' + n + r'}\\'
        for j, n in enumerate(notes))
    return (head + '\n'.join(lines) + '\n' + r'\hline' + '\n' + foot
            + '\n' + r'\end{tabular}')


def main():
    rc, mq, pf, pp, mc = load()
    body = ablation(rc, mq, pf, pp, mc)
    tex = r"""% 자동 생성 — make_si_tables.py. 손으로 고치지 말 것.
% 수치 출처: outputs/rollback_c.json (패널 3구성),
%            outputs/pesticide_family.json + pesticide_panel.json (농약 구성),
%            outputs/multiplicity_across_panels.json (구성별 최소 q).
\begin{table}[htb]
\centering
\caption{Every exact permutation test across panel constructions, descriptor
spaces and statistics. Values are one-sided $p$. Bold marks $p<0.05$ before
correction. The last column is the smallest Benjamini--Hochberg value over the
eighteen tests of that panel construction. Dimensions are the retained
coordinate counts after the correlation filter of
Section~\ref{meth:orderinvariance}.}
\label{tab:ablation}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llcccccccc}
\hline
Panel / simulant set & Space & Dim. & $\Delta\bar{d}$ & $\Delta\tilde{d}$ &
$\overline{\ln r}$ & $\bar{d}_\mathrm{sim}$ & $n_\mathrm{closer}$ &
$n_\mathrm{valid}$ & min $q$ \\
\hline
""" + body + r"""
\hline
\multicolumn{10}{l}{\footnotesize Statistics, in order: mean and median of
$d_\mathrm{sim}-d_\mathrm{intf}$; mean $\ln(d_\mathrm{sim}/d_\mathrm{intf})$;
mean $d_\mathrm{sim}$;}\\
\multicolumn{10}{l}{\footnotesize count of agents with a simulant closer than
every interferent; count of valid matches.}\\
\end{tabular}
\end{table}
"""
    os.makedirs('outputs', exist_ok=True)
    open('outputs/tab_ablation.tex', 'w').write(tex)
    print('-> outputs/tab_ablation.tex')

    hdr = ('% 자동 생성 — make_si_tables.py. 손으로 고치지 말 것.\n')
    open('outputs/tab_family.tex', 'w').write(hdr + family(mq) + '\n')
    open('outputs/tab_registers.tex', 'w').write(
        hdr + '% 출처: sources/registers_20261106.csv (URL 포함)\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{Regulatory acute '
        'inhalation values held for the three incumbent simulants, by '
        'register. Sources and collection dates are in '
        '\\texttt{sources/registers\\_20261106.csv} of the deposit. EPA '
        'AEGL is the axis the agents are placed on in the main text; the '
        'other registers use different averaging times and severity '
        'definitions and are not interchangeable with it.}\n'
        '\\label{tab:registers}\n\\small\n'
        + registers() + '\n\\end{table}\n')
    print('-> outputs/tab_registers.tex')
    open('outputs/tab_designations.tex', 'w').write(
        hdr + '% 출처: sources/designations_20261106.csv (URL 포함)\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{Every simulant '
        'designation we located, with the agent it names, the measurement it '
        'was made against, and the basis the designating document itself '
        'states. Source URLs and the search date are in '
        '\\texttt{sources/designations\\_20261106.csv} of the deposit.}\n'
        '\\label{tab:designations}\n\\footnotesize\n'
        + designations() + '\n\\end{table}\n')
    print('-> outputs/tab_designations.tex')
    open('outputs/tab_dyepanel.tex', 'w').write(
        hdr + '% 출처: sources/dye_panel_20261106.csv\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{The 29-component '
        'array, in the order the three source studies number it. Five entries '
        'marked $^{a}$ are not chromogenic dyes; $^{b}$ marks a component '
        'whose stock solvent differs between the source studies.}\n'
        '\\label{tab:dyepanel}\n\\small\n'
        + dyepanel() + '\n\\end{table}\n')
    print('-> outputs/tab_dyepanel.tex')
    print('-> outputs/tab_family.tex')
    open('outputs/tab_blockbalance.tex', 'w').write(hdr + blockbalance() + '\n')
    print('-> outputs/tab_blockbalance.tex')
    open('outputs/tab_sparse.tex', 'w').write(hdr + sparse() + '\n')
    print('-> outputs/tab_sparse.tex')
    open('outputs/tab_plocal.tex', 'w').write(hdr + plocal() + '\n')
    print('-> outputs/tab_plocal.tex')
    open('outputs/tab_volatility.tex', 'w').write(hdr + volatility() + '\n')
    print('-> outputs/tab_volatility.tex')
    open('outputs/tab_vpairing.tex', 'w').write(hdr + vpairing() + '\n')
    print('-> outputs/tab_vpairing.tex')
    open('outputs/tab_premisenull.tex', 'w').write(hdr + premisenull() + '\n')
    print('-> outputs/tab_premisenull.tex')
    for i, (a, b, c, key) in enumerate(ROWS):
        if key == 'pesticide':
            dims = {s: pp[s]['n_dims'] for s in pp}
        else:
            dims = {s: rc[key][s]['n_descriptors'] for s in rc[key]}
        print(f'   {key:<40} dims {dims}')


if __name__ == '__main__':
    main()
