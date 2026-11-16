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


_NUMWORD = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five', 6: 'Six',
            7: 'Seven', 8: 'Eight', 9: 'Nine', 10: 'Ten'}


def _nondye_count_word():
    """비염료 표식 수를 dye_panel_audit.json 에서 읽는다.

    2026-11-13: 캡션에 'Five' 가 박혀 있었는데 감사 산출물은 6 이었다.
    손으로 센 숫자를 캡션에 쓰지 않는다.
    """
    with open('outputs/dye_panel_audit.json', encoding='utf-8') as fh:
        n = json.load(fh)['n_marked_nondye']
    return _NUMWORD[n]

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
    # 2026-11-12. 휘발성 축을 표에 넣는다. 캡션이 'Every exact permutation
    # test' 라고 적는데 이 축이 빠져 있었다 — 사전지정 전수 검정이고 통계량
    # 여섯 개가 모두 계산돼 있다(outputs/volatility_space.json). 1차원이므로
    # 공간 세 행이 아니라 한 행이다.
    # 2026-11-13 (GATE 8 11차). min q 칸을 비워 두었더니, 독자가 훑는 유일한
    # 표에서 굵은 값 넷이 든 행만 보정이 정지한 모습이 됐다. 이 축 자체의
    # 6통계량 BH 는 이미 family_corrections 에 예치돼 있으므로 기입한다.
    vd = [t for t in _vol()['designation']
          if t.get('label') == 'panel28_incumbent3' and not t.get('skipped')]
    if len(vd) != 1:
        raise SystemExit('[tab_ablation] 휘발성 designation 항목이 1개가 아니다')
    vst = vd[0]['statistics']
    vcells = ' & '.join(fmt(vst[st]['p']) for st in STATS)
    vq = _vol()['family_corrections']['axis_six_statistics_bh']
    lines.append(f"{'28 cpd., incumbent 3':<22} & {'volatility':<11} & "
                 f"{vd[0].get('n_dims', 1):>2} & {vcells} & "
                 f"{vq:.3f}$^{{\\dagger}}$ \\\\")
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


# ------------------------------------------------------- tab:failaudit
# 2026-11-08 신설 (에디터 게이트 1회차, (c)항). designation_failure_audit.json
# 은 판독 3종 x 필터 2종 x 패널 3종 + 희소 인코딩 7종, 모두 25셀에서
# 매칭집합을 집계하는데 원고가 어디에서도 인용하지 않았다. 게이트는 이것을
# "논지를 가장 강하게 지지하는 자료인데 빠져 있다" 로 지적했다. 손표로
# 두지 않고 생성한다.
FA_CODE = {'sarin (GB)': 'GB', 'soman (GD)': 'GD', 'tabun (GA)': 'GA',
           'cyclosarin (GF)': 'GF', 'sulfur mustard (HD)': 'HD',
           'lewisite 1 (L)': 'L', 'hydrogen cyanide (AC)': 'AC',
           'cyanogen chloride (CK)': 'CK', 'phosgene (CG)': 'CG',
           'chloropicrin (PS)': 'PS', 'VX': 'VX'}


def failaudit():
    fa = json.load(open('outputs/designation_failure_audit.json'))
    n = fa['n_cells']
    rows = sorted(fa['agents'].items(), key=lambda kv: -kv[1]['n_match'])
    lines = []
    for name, b in rows:
        code = FA_CODE.get(name, name)
        nm = b['n_match']
        note = 'never' if nm == 0 else ('always' if nm == n else '')
        lines.append(f'    {code:<4} & {name:<24} & {nm:>2} / {n} & '
                     f'{note} \\\\')
    head = (r'  \begin{tabular}{llrl}' '\n' r'    \hline' '\n'
            r'    Agent & & Cells matched & \\' '\n'
            r'    \hline' '\n')
    return (head + '\n'.join(lines) + '\n' + r'    \hline' + '\n'
            + r'  \end{tabular}')


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
        r'(Section~\ref{meth:spaces}).}\\' + '\n'
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
            r'\multicolumn{5}{l}{\footnotesize value the source itself marks as extrapolated rather than}\\' '\n'
            r'\multicolumn{5}{l}{\footnotesize measured, used only for solvents. Over this compound}\\' '\n'
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



CAP_DESIG = (
    '\\caption{The simulant designations we located in the instruments of '
    'Section~\\ref{res:searchscope}, with the agent each names, the '
    'measurement it was made against, and the basis the designating document '
    'itself states. Compounds named in a document without being characterised '
    'there are marked as listed rather than measured. This is what the search '
    'returned, not a census of designations. Source URLs and the date each '
    'document was read are in '
    '\\texttt{sources/designations\\_20261106.csv} of the deposit.}')


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
                 # 2026-11-16. ASTM F2588 행의 사유 세 개. 일반 접두 치환이
                 # 먼저 걸리면 한글이 남아 조판이 깨진다.
                 .replace('[수집 불가: 표준 본문 미열람]',
                          '[not obtained: standard text not read]')
                 .replace('[수집 불가: 공개 URL 없음]',
                          '[not obtained: no public URL]')
                 .replace('[미열람]', '[not read]')
                 .replace('[수집 불가', '[not obtained')
                 .replace('특허 본문 미열람 — Google Patents, Justia, USPTO '
                          '전문 서비스 모두 접근 거부. 서지사항만 확인',
                          'patent text not retrieved; Google Patents, Justia '
                          'and the USPTO full-text service all refused; '
                          'bibliographic record only')
                 .replace('특허 본문 미열람 — Google Patents 및 Justia 접근 거부',
                          'patent text not read; Google Patents and Justia refused')
                 .replace('특허 본문 미열람', 'patent text not retrieved')
                 .replace('공개 기록에서 지정 문서를 찾지 못함',
                          'none in the public record'))

    lines = []
    for r in rows:
        lines.append(' & '.join(esc(r[k]) for k in
                                ('compound', 'document', 'agent',
                                 'measurement', 'stated_basis')) + r' \\')
    # 2026-11-16. 3행(malathion, methyl salicylate 두 건)을 더하자 SI 한
    # 면보다 137pt 길어졌다. longtable 은 float 가 아니므로 면을 넘겨
    # 조판된다. 캡션과 라벨을 longtable 안으로 옮긴다.
    _h = ('Compound & Document & Agent & Measurement & '
          'Basis the document states')
    head = (r'\begin{longtable}{>{\raggedright\arraybackslash}p{2.0cm}>{\raggedright\arraybackslash}p{3.6cm}>{\raggedright\arraybackslash}p{1.6cm}>{\raggedright\arraybackslash}p{2.4cm}>{\raggedright\arraybackslash}p{4.6cm}}'
            '\n' + CAP_DESIG + '\n'
            r'\label{tab:designations}\\' '\n'
            r'\hline' '\n' + _h + r' \\' '\n'
            r'\hline' '\n' r'\endfirsthead' '\n'
            r'\hline' '\n' + _h + r' \\' '\n'
            r'\hline' '\n' r'\endhead' '\n'
            r'\hline' '\n' r'\endfoot' '\n')
    return head + '\n'.join(lines) + '\n' + r'\end{longtable}'


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




def _niosh_nn():
    """표 캡션의 정규화자 두 값을 예치물에서 읽는다 (2026-11-14).

    종전에는 캡션 문자열에 3.399 와 3.042 를 하드코딩했다. 비교군을 네 종에서
    열 종으로 고치자 뒤 값이 3.197 로 바뀌었는데 캡션은 그대로였고, 검증기가
    잡았다. 하드코딩을 제거한다.
    """
    d = json.load(open('outputs/niosh_panel.json'))
    p = d['panels']['panel28_incumbent3']
    return (f"{p['baseline']['nn_median']:.3f} before the addition and "
            f"{p['plus_all_recomputed']['nn_median']:.3f} after it")



def nioshlist():
    """NIOSH 후보목록 11종과 패널 소속 (2026-11-14)."""
    d = json.load(open('outputs/niosh_panel.json'))
    L = ['\\begin{tabular}{llc}', '\\hline',
         'Compound & Entered as & On a panel \\\\', '\\hline']
    for n in d['added_to_comparison_pool']:
        iup = d['iupac_like'][n]
        sm = d['smiles'][n].replace('(', '\\allowbreak (')
        L.append(f"{iup} & \\texttt{{{sm}}} & no \\\\")
    L += ['dimethyl methanephosphonate & \\texttt{COP(C)\\allowbreak (=O)OC} '
          '& yes \\\\',
          '\\hline', '\\end{tabular}']
    return '\n'.join(L)



def volniosh():
    """휘발성 축 비교풀 복원. 세 조건 x 여섯 통계량 (2026-11-16)."""
    d = json.load(open('outputs/volatility_niosh.json'))
    LAB = [('baseline', 'as reported'),
           ('plus3_dated', 'with the three dated'),
           ('plus5_with_undated', 'with all five')]
    STATS = [('mean_diff', 'mean difference'),
             ('median_diff', 'median difference'),
             ('mean_logratio', 'mean log ratio'),
             ('mean_dsim', 'mean distance to a simulant'),
             ('n_closer', 'number closer'),
             ('n_valid', 'valid matches')]
    L = ['\\begin{tabular}{lrrrrrr}', '\\hline',
         ' & \\multicolumn{2}{c}{as reported}'
         ' & \\multicolumn{2}{c}{with the three dated}'
         ' & \\multicolumn{2}{c}{with all five} \\\\',
         '\\cline{2-3}\\cline{4-5}\\cline{6-7}',
         'Statistic & observed & $p$ & observed & $p$ & observed & $p$ \\\\',
         '\\hline']
    for key, name in STATS:
        row = [name]
        for lab, _ in LAB:
            v = d['runs'][lab]['statistics'][key]
            row.append(f"{v['observed']:.5f}")
            row.append(f"{v['p']:.5f}")
        L.append(' & '.join(row) + ' \\\\')
    L.append('\\hline')
    r = [d['runs'][lab] for lab, _ in LAB]
    L.append('compounds in pool & \\multicolumn{2}{c}{'
             + str(r[0]['n_nonagent']) + '} & \\multicolumn{2}{c}{'
             + str(r[1]['n_nonagent']) + '} & \\multicolumn{2}{c}{'
             + str(r[2]['n_nonagent']) + '} \\\\')
    L.append('designations enumerated & \\multicolumn{2}{c}{'
             + str(r[0]['n_permutations']) + '} & \\multicolumn{2}{c}{'
             + str(r[1]['n_permutations']) + '} & \\multicolumn{2}{c}{'
             + str(r[2]['n_permutations']) + '} \\\\')
    L += ['\\hline', '\\end{tabular}']
    return '\n'.join(L)

def niosh():
    """NIOSH 후보 4종을 비교군에 넣을 때 작용제별 최근접 거리 (2026-11-14)."""
    d = json.load(open('outputs/niosh_panel.json'))
    p = d['panels']['panel28_incumbent3']
    b, e = p['baseline'], p['plus_all_recomputed']
    nb, ne = b['nn_median'], e['nn_median']
    L = ['\\begin{tabular}{lrrlrrl}', '\\hline',
         ' & \\multicolumn{3}{c}{twenty-eight-compound pool}'
         ' & \\multicolumn{3}{c}{with the ten added} \\\\',
         '\\cline{2-4}\\cline{5-7}',
         'Agent & $d_{\\mathrm{sim}}$ & $d_{\\mathrm{near}}$ & nearest '
         'non-designated & $d_{\\mathrm{sim}}$ & $d_{\\mathrm{near}}$ & '
         'nearest non-designated \\\\', '\\hline']
    for x, y in zip(b['agents'], e['agents']):
        # 생거리로 적는다. 정규화해 적으면 본문 표 2 의 병합공간 d_norm 과
        # 같은 단위로 보이는데 공간이 다르다 — 이 논문이 경고하는 혼동이다.
        L.append(f"{x['agent']} & {x['d_sim']:.3f} & "
                 f"{x['d_blocker']:.3f} & {x['blocker']} & "
                 f"{y['d_sim']:.3f} & {y['d_blocker']:.3f} & "
                 f"{y['blocker']} \\\\")
    L += ['\\hline',
          f"matched & \\multicolumn{{3}}{{c}}{{{b['n_valid']} of 11}}"
          f" & \\multicolumn{{3}}{{c}}{{{e['n_valid']} of 11}} \\\\",
          'pre-specified $p$ & '
          f"\\multicolumn{{3}}{{c}}{{{b['p']['mean_diff']:.3f}}} & "
          f"\\multicolumn{{3}}{{c}}{{{e['p']['mean_diff']:.3f}}} \\\\",
          '\\hline', '\\end{tabular}']
    return '\n'.join(L)



def voldecomp():
    """휘발성 축 지정 규모별 순위와 피복 사다리 (2026-11-14)."""
    d = json.load(open('outputs/volatility_decomposition.json'))
    sd, sp = d['sidedness'], d['span']
    L = ['\\begin{tabular}{llrrr}', '\\hline',
         'Size & Designated set & $\\overline{\\Delta d}$ & Rank & '
         'One-sided $p$ \\\\', '\\hline']
    for k, n in (('k1', 'single'), ('k2', 'pair')):
        rows = d['by_designation_size'][k]['ranking']
        tot = d['by_designation_size'][k]['n_combinations']
        inc = {'k1': ['DMMP', 'dichloromethane', 'DPGME'],
               'k2': [' + '.join(sorted(x)) for x in
                      (('DMMP', 'dichloromethane'),
                       ('DMMP', 'DPGME'),
                       ('dichloromethane', 'DPGME'))]}[k]
        pick = [rows[0]] + [r for r in rows
                            if ' + '.join(r['compounds']) in inc
                            and r is not rows[0]]
        seen = set()
        for r in pick:
            nm = ' + '.join(r['compounds'])
            if nm in seen:
                continue
            seen.add(nm)
            tag = n if nm == ' + '.join(rows[0]['compounds']) else ''
            L.append(f"{tag} & {nm} & ${r['mean_diff']:+.4f}$ & "
                     f"{r['rank']} of {tot} & {r['p']:.4f} \\\\")
        L.append('\\hline')
    inc3 = d['by_designation_size']['k3']['incumbent']
    L.append(f"triple & {' + '.join(inc3['compounds'])} & "
             f"${inc3['mean_diff']:+.4f}$ & {inc3['rank']} of "
             f"{d['by_designation_size']['k3']['n_combinations']} & "
             f"{inc3['p']:.4f} \\\\")
    best3 = d['by_designation_size']['k3']['ranking'][0]
    L.append(f" & {' + '.join(best3['compounds'])} & "
             f"${best3['mean_diff']:+.4f}$ & {best3['rank']} of "
             f"{d['by_designation_size']['k3']['n_combinations']} & "
             f"{best3['p']:.4f} \\\\")
    L += ['\\hline',
          '\\multicolumn{5}{@{}p{0.92\\textwidth}@{}}{\\footnotesize Two-sided readings of the '
          'incumbent triple: '
          f"{sd['two_sided_2min']:.4f} (twice the smaller tail), "
          f"{sd['two_sided_abs']:.4f} (deviation from the null mean).}}\\\\",
          '\\multicolumn{5}{@{}p{0.92\\textwidth}@{}}{\\footnotesize Spearman correlation between '
          'axis coverage and the statistic over all '
          f"{sp['n']} designations: ${sp['spearman_rho']:.3f}$. The incumbent "
          f"span is {sp['incumbent_span']:.3f}, the largest in the pool.}}"
          '\\\\']
    for r in sp['ladder']:
        L.append('\\multicolumn{5}{@{}p{0.92\\textwidth}@{}}{\\footnotesize Conditioned on span '
                 f"$\\geq$ {r['threshold']:.3f} (top "
                 f"{100 * (1 - r['quantile']):.0f} per cent, $n = {r['n']}$): "
                 f"$p = {r['conditional_p']:.4f}$, rank {r['rank']} of "
                 f"{r['n']}.}}\\\\")
    L.append('\\end{tabular}')
    return '\n'.join(L)


def main():
    rc, mq, pf, pp, mc = load()
    body = ablation(rc, mq, pf, pp, mc)
    tex = r"""% 자동 생성 — make_si_tables.py. 손으로 고치지 말 것.
% 수치 출처: outputs/rollback_c.json (패널 3구성),
%            outputs/pesticide_family.json + pesticide_panel.json (농약 구성),
%            outputs/multiplicity_across_panels.json (구성별 최소 q).
\begin{table}[htb]
\centering
\caption{The exact permutation tests over the four panel and designation
constructions, by descriptor space and statistic. Values are one-sided $p$.
Bold marks $p<0.05$ before correction. The last column is the smallest
Benjamini--Hochberg value over the eighteen tests of that panel construction.
$^{\dagger}$The volatility axis was specified separately and is not a member of
any of those families; the value given for it is the Benjamini--Hochberg
minimum over its own six statistics, and
Section~\ref{M-res:vol} gives what it becomes under wider families.
Dimensions are the retained coordinate counts after the correlation filter of
Section~\ref{meth:orderinvariance}; the volatility axis is a single measured
quantity and carries no filter. Four further enumerations are reported in the
text rather than here: the two-pesticide and acetonitrile relabellings of
Section~\ref{M-res:panels}
(\texttt{outputs/twolabel\_panel.json}, \texttt{outputs/acetonitrile\_panel.json});
the enlargement of the comparison pool by the ten further compounds the
NIOSH candidate list names beside DMMP (Section~\ref{meth:niosh},
\texttt{outputs/niosh\_panel.json}); and the single- and pair-designation
enumerations on the volatility axis (Section~\ref{meth:vol},
\texttt{outputs/volatility\_decomposition.json}).}
\label{tab:ablation}
% 2026-11-13. small/4pt 로는 29.7pt 가 여백을 넘었다 (GATE 8 11차).
\footnotesize
\setlength{\tabcolsep}{3pt}
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
        '\\label{tab:registers}\n\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n'
        + registers() + '\n\\end{table}\n')
    print('-> outputs/tab_registers.tex')
    open('outputs/tab_designations.tex', 'w').write(
        hdr + '% 출처: sources/designations_20261106.csv (URL 포함)\n'
        + '\\begingroup\n\\scriptsize\n'
        # 2026-11-16. 한 셀의 괄호 사유 문구가 2.55pt 넘친다. raggedright
        # 로는 안 잡히므로 단어 간격을 느슨하게 둔다.
        '\\sloppy\n'
        '\\setlength{\\tabcolsep}{3pt}\n'
        '\\renewcommand{\\arraystretch}{1.0}\n'
        + designations() + '\n\\endgroup\n')
    print('-> outputs/tab_designations.tex')
    open('outputs/tab_dyepanel.tex', 'w').write(
        hdr + '% 출처: sources/dye_panel_20261106.csv\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{The 29-component '
        'array, in the order the three source studies number it. '
        + _nondye_count_word() + ' entries '
        'marked $^{a}$ are not chromogenic dyes; $^{b}$ marks a component '
        'whose stock solvent differs between the source studies.}\n'
        '\\label{tab:dyepanel}\n\\small\n'
        + dyepanel() + '\n\\end{table}\n')
    print('-> outputs/tab_dyepanel.tex')
    print('-> outputs/tab_family.tex')
    open('outputs/tab_blockbalance.tex', 'w').write(hdr + blockbalance() + '\n')
    print('-> outputs/tab_blockbalance.tex')
    open('outputs/tab_sparse.tex', 'w').write(hdr + sparse() + '\n')
    open('outputs/tab_failaudit.tex', 'w').write(hdr + failaudit() + '\n')
    print('-> outputs/tab_failaudit.tex')
    print('-> outputs/tab_sparse.tex')
    open('outputs/tab_plocal.tex', 'w').write(hdr + plocal() + '\n')
    print('-> outputs/tab_plocal.tex')
    open('outputs/tab_volatility.tex', 'w').write(hdr + volatility() + '\n')
    print('-> outputs/tab_volatility.tex')
    open('outputs/tab_vpairing.tex', 'w').write(hdr + vpairing() + '\n')
    print('-> outputs/tab_vpairing.tex')
    open('outputs/tab_premisenull.tex', 'w').write(hdr + premisenull() + '\n')
    print('-> outputs/tab_premisenull.tex')
    open('outputs/tab_volniosh.tex', 'w').write(
        hdr + '% 출처: outputs/volatility_niosh.json\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{The designation '
        'test on the volatility axis before and after the compounds the NIOSH '
        'candidate list names beside DMMP enter the comparison pool. The '
        'pre-specified statistic is the mean difference. Three of the ten '
        'further compounds carry a tier-one 25~$^{\\circ}$C vapour pressure; '
        'two more are carried by a source that states no temperature and enter '
        'the third column only; five return no entry and are not estimated. '
        'Section~\\ref{meth:volniosh} gives each adopted value and its '
        'source.}\n'
        '\\label{tab:volniosh}\n\\footnotesize\n'
        '\\setlength{\\tabcolsep}{3pt}\n'
        + volniosh() + '\n\\end{table}\n')
    print('-> outputs/tab_volniosh.tex')
    open('outputs/tab_niosh.tex', 'w').write(
        hdr + '% 출처: outputs/niosh_panel.json\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{Distance from '
        'each agent to its nearest designated simulant and to its nearest '
        'non-designated compound, in the structural space, before and after '
        'the further compounds the NIOSH candidate list names beside DMMP enter '
        'the comparison pool. Distances are unnormalised Euclidean distances '
        'in the standardised structural space, not the merged-space normalised '
        'distances of the main text, and the two are not comparable; the '
        'median nearest-neighbour separation that sets the match threshold is '
        + _niosh_nn() +
        '. The designated set is '
        'the incumbent three throughout. Values are in '
        '\\texttt{outputs/niosh\\_panel.json}.}\n'
        '\\label{tab:niosh}\n\\footnotesize\n'
        '\\setlength{\\tabcolsep}{3pt}\n'
        + niosh() + '\n\\end{table}\n')
    print('-> outputs/tab_niosh.tex')
    open('outputs/tab_nioshlist.tex', 'w').write(
        hdr + '% 출처: outputs/niosh_panel.json\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{The eleven '
        'compounds the NIOSH candidate list for barrier-material permeation '
        'testing names, the structures we entered for the ten that are on '
        'neither panel, and which one is. Names follow the source; its '
        'spelling of triethyl phosphate is corrected here. Values are in '
        '\\texttt{outputs/niosh\\_panel.json}.}\n'
        '\\label{tab:nioshlist}\n\\footnotesize\n'
        '\\setlength{\\tabcolsep}{4pt}\n'
        + nioshlist() + '\n\\end{table}\n')
    print('-> outputs/tab_nioshlist.tex')
    open('outputs/tab_voldecomp.tex', 'w').write(
        hdr + '% 출처: outputs/volatility_decomposition.json\n'
        + '\\begin{table}[htb]\n\\centering\n\\caption{The volatility '
        'axis taken apart. For each designation size, the best set in the '
        'pool and the incumbent sets of that size, ranked on the '
        'pre-specified mean difference (negative is better). Two-sided '
        'readings and the coverage conditioning follow. Values are in '
        '\\texttt{outputs/volatility\\_decomposition.json}.}\n'
        '\\label{tab:voldecomp}\n\\footnotesize\n'
        '\\setlength{\\tabcolsep}{4pt}\n'
        + voldecomp() + '\n\\end{table}\n')
    print('-> outputs/tab_voldecomp.tex')
    for i, (a, b, c, key) in enumerate(ROWS):
        if key == 'pesticide':
            dims = {s: pp[s]['n_dims'] for s in pp}
        else:
            dims = {s: rc[key][s]['n_descriptors'] for s in rc[key]}
        print(f'   {key:<40} dims {dims}')


if __name__ == '__main__':
    main()
