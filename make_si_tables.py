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
    ('35 cpd., incumbent', '\\quad + literature 8', '',
     'panel35_designated_plus_literature8'),
    ('28 cpd., incumbent 3', '\\quad (omits the seven)', '',
     'panel28_designated3'),
    ('35 cpd., incumbent 3', '\\quad + literature 5',
     '\\quad + pesticides 2', 'pesticide'),
]


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
    return rc, mq, pf, pp


def ablation(rc, mq, pf, pp):
    MQ = mq['min_q_by_config']
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
            ('panel35_incumbent_plus_lit8', '35 cpd.,', '+ literature 8', ''),
            ('panel28_incumbent3', '28 cpd.,', 'incumbent 3',
             '(omits the seven)'),
            ('pesticide', '35 cpd.,', 'relabelled 10', '(pesticides in)')]
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
           ('Thirty-five, plus literature 8', 'panel35_incumbent_plus_lit8'),
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


def main():
    rc, mq, pf, pp = load()
    body = ablation(rc, mq, pf, pp)
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
    print('-> outputs/tab_family.tex')
    open('outputs/tab_blockbalance.tex', 'w').write(hdr + blockbalance() + '\n')
    print('-> outputs/tab_blockbalance.tex')
    open('outputs/tab_sparse.tex', 'w').write(hdr + sparse() + '\n')
    print('-> outputs/tab_sparse.tex')
    open('outputs/tab_plocal.tex', 'w').write(hdr + plocal() + '\n')
    print('-> outputs/tab_plocal.tex')
    for i, (a, b, c, key) in enumerate(ROWS):
        if key == 'pesticide':
            dims = {s: pp[s]['n_dims'] for s in pp}
        else:
            dims = {s: rc[key][s]['n_descriptors'] for s in rc[key]}
        print(f'   {key:<40} dims {dims}')


if __name__ == '__main__':
    main()
