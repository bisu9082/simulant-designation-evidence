#!/usr/bin/env python3
"""순열 표(Table S1)를 산출물에서 직접 생성한다 (2026-10-24 신설).

왜 필요한가. 이 표는 4구성 x 3공간 x 6통계 = 72 칸이고 굵은 글씨가 p<0.05 를
뜻한다. 분석을 다시 돌릴 때마다 손으로 고쳐 왔는데, 이 세션에서 차원 수(38 대
36)와 굵은 표시가 값과 어긋난 채로 남는 일이 생겼다. 값과 표시를 같은 곳에서
만든다.

표 본문만 찍는다. 캡션과 tabular 선언은 원고가 들고 있고, splice_tables.py 가
\\hline 사이를 바꿔 끼운다.
"""
import json

ORDER = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']
SPACES = [('structural', 'structural'), ('electronic3', 'electronic'),
          ('merged', 'merged')]

ROWS = [
    ('panel28_designated3',
     ['28 cpd., incumbent 3', '\\quad (omits the five)', '']),
    ('panel35_designated3',
     ['\\textbf{35 cpd., incumbent 3}', '\\textbf{(primary)}', '']),
    ('panel35_designated_plus_literature8',
     ['35 cpd., incumbent', '\\quad + literature 8', '']),
    ('panel35_relabelled10',
     ['35 cpd., incumbent 3', '\\quad + literature 5', '\\quad + pesticides 2']),
]
# 원고의 표 순서: primary 가 먼저다.
ROWS = [ROWS[1], ROWS[2], ROWS[0], ROWS[3]]


def fmt(p):
    if p < 0.001:
        s = '$<$0.001'
    elif p < 0.01:
        s = f'{p:.3f}'
    elif p < 0.1:
        s = f'{p:.4f}'.rstrip('0')
        if len(s.split('.')[1]) < 3:
            s = f'{p:.3f}'
    else:
        s = f'{p:.3f}'
    return f'\\textbf{{{s}}}' if p < 0.05 else s


def main():
    r = json.load(open('outputs/rollback_c.json'))['panels']
    # 농약 재라벨 구성은 별도 스크립트가 낸다. 칸 구조가 달라 맞춰 넣는다.
    pp = json.load(open('outputs/pesticide_panel.json'))
    pf = json.load(open('outputs/pesticide_family.json'))
    r['panel35_relabelled10'] = {
        k: dict(n_descriptors=v.get('n_dims', v.get('n_kept', 0)),
                statistics={s: {'p': v['p'][s]} for s in ORDER})
        for k, v in pp['spaces'].items()}
    mult = {m[0]: m for m in json.load(open('outputs/multiplicity_c.json'))}
    mult['panel35_relabelled10'] = ['panel35_relabelled10', 18, pf['min_p'],
                                    None, pf['bh_all_min_q'], 0]
    out = []
    for tag, labels in ROWS:
        src = r.get(tag)
        if src is None:
            continue
        minq = mult[tag][4] if tag in mult else None
        for k, (skey, sname) in enumerate(SPACES):
            cell = src.get(skey) if isinstance(src, dict) else None
            if cell is None or 'statistics' not in cell:
                continue
            st = cell['statistics']
            cells = '  & '.join(fmt(st[s]['p']) for s in ORDER)
            q = ''
            if k == 2 and minq is not None:
                q = fmt(minq) if minq < 0.05 else f'{minq:.3f}'
            out.append(f"{labels[k]:22s} & {sname:11s} & "
                       f"{cell['n_descriptors']:2d} & {cells} & {q} \\\\")
        out.append('\\hline')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
