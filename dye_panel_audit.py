#!/usr/bin/env python3
"""29 성분 패널의 채널 비대칭 (2026-11-07 신설, 12차 A-9).

왜 필요한가. SI meth:deltaE 가 "Five of the twenty-nine are not chromogenic
dyes" 라고 적고 그 다섯을 sources/dye_panel_20261106.csv 의 note_kind 로
표시해 둔다. 그런데 그 분류는 손으로 적힌 것이고 어느 측정에도 걸려 있지
않았다. 자외 채널과 반사 채널의 응답 크기 비를 재면 형광체는 바로 드러난다.
형광은 여기광으로 켜지므로 자외 쪽에만 큰 응답을 낸다.

재보면 퀴닌이 29개 중 가장 큰 비를 낸다. 퀴닌은 형광 양자수율 표준물질이고,
이미 형광체로 표시된 안트라센·피렌보다도 비가 크다. 표시에서 빠져 있었다.

여기서 세는 것은 '비색 염료가 아닌 성분'의 개수가 아니라 그 비 자체다.
비는 측정이고 분류는 판단이므로, 비를 산출물로 남기고 판단은 SI 가 적는다.

출력: outputs/dye_panel_audit.json
"""
import json
import os

import numpy as np
import pandas as pd

PANEL = 'sources/dye_panel_20261106.csv'
RUNS = ('d1', 'd2', 'd3')
# 한 셀에서 최대 ΔE00 가 이 값 아래면 사실상 응답이 없는 채널로 센다.
DEAD_CUT = 1.0


def main():
    pan = pd.read_csv(PANEL)
    name = dict(zip(pan['index'], pan['name']))
    kind = dict(zip(pan['index'], pan['note_kind'].fillna('')))
    note = dict(zip(pan['index'], pan['note'].fillna('')))

    per = {}
    for tag in RUNS:
        try:
            u = pd.read_csv(f'deltaE/{tag}.uv_deltaE_pivot.csv')
            n = pd.read_csv(f'deltaE/{tag}.normal_deltaE_pivot.csv')
        except OSError as e:
            print(f'  [수집 불가: {tag} {e}]')
            continue
        cols = [c for c in u.columns
                if c in n.columns and c.lower().startswith('dye')]
        for c in cols:
            i = int(c.lower().replace('dye', ''))
            mu, mn = float(u[c].max()), float(n[c].max())
            r = mu / mn if mn > 0 else float('inf')
            d = per.setdefault(i, dict(index=i, name=name.get(i, '?'),
                                       note_kind=kind.get(i, ''),
                                       uv_max=0.0, normal_max=0.0,
                                       ratio_max=0.0, runs=[]))
            d['uv_max'] = max(d['uv_max'], mu)
            d['normal_max'] = max(d['normal_max'], mn)
            d['ratio_max'] = max(d['ratio_max'], r)
            d['runs'].append(tag)

    # 2026-11-07 (12차 R1-F9). 셀별로 거의 응답하지 않는 채널을 센다.
    # 원고는 29 성분을 한 덩어리로 서술하는데, 어떤 셀에서는 그 가운데
    # 여럿이 측정 잡음 수준의 신호만 낸다.
    dead = {}
    for tag in RUNS:
        for mode in ('uv', 'normal'):
            try:
                d = pd.read_csv(f'deltaE/{tag}.{mode}_deltaE_pivot.csv')
            except OSError:
                continue
            cs = [c for c in d.columns if c.lower().startswith('dye')]
            # 정확히 0 인 열은 SI res:arraydetail 이 이미 보고하고
            # 축약이 제거한다. 보탬이 되는 것은 0 은 아니면서 잡음
            # 수준인 열이다. 그쪽은 축약에 걸리지 않고 코사인에 들어간다.
            zz = [int(c.lower().replace('dye', '')) for c in cs
                  if float(d[c].abs().max()) == 0.0]
            dd = [int(c.lower().replace('dye', '')) for c in cs
                  if 0.0 < float(d[c].max()) < DEAD_CUT]
            dead[f'{tag}.{mode}'] = dict(
                n_channels=len(cs), n_zero=len(zz), n_dead=len(dd),
                zero=[name.get(i, str(i)) for i in sorted(zz)],
                dead=[name.get(i, str(i)) for i in sorted(dd)])

    rows = sorted(per.values(), key=lambda d: -d['ratio_max'])
    marked = [d for d in rows if d['note_kind'] == 'nondye']
    # 형광체 판정은 note 열의 문구로 한다. 종전에는 note_kind 를 보고
    # 있었고 그 열에는 'nondye' 만 들어 있어서, 이름 목록 하드코딩에
    # 기대고 있었다 (2026-11-07, 12차 A-9).
    fluor = [d for d in marked
             if 'fluorophore' in str(note.get(d['index'], '')).lower()]

    print('=' * 74)
    print('[패널 감사] 자외 최대 / 반사 최대 (성분별, 세 런의 최대)')
    print('=' * 74)
    print(f'{"번호":>4} {"성분":<34} {"표시":<8} {"자외":>8} {"반사":>8} {"비":>7}')
    for d in rows:
        print(f'{d["index"]:>4} {d["name"][:33]:<34} {d["note_kind"][:7]:<8} '
              f'{d["uv_max"]:>8.2f} {d["normal_max"]:>8.2f} '
              f'{d["ratio_max"]:>7.1f}')

    top = rows[0]
    marked_fluor_max = max((d['ratio_max'] for d in fluor), default=0.0)
    out = dict(
        generated_by='dye_panel_audit.py',
        n_components=len(rows),
        n_marked_nondye=len(marked),
        marked_nondye=[d['name'] for d in marked],
        highest_ratio=dict(index=top['index'], name=top['name'],
                           uv_max=top['uv_max'],
                           normal_max=top['normal_max'],
                           ratio=top['ratio_max'],
                           marked=bool(top['note_kind'] == 'nondye')),
        marked_fluorophore_max_ratio=marked_fluor_max,
        n_fluorophores=len(fluor),
        fluorophores=[d['name'] for d in fluor],
        dead_cut=DEAD_CUT,
        dead_channels=dead,
        n_dead_max=max((v['n_dead'] for v in dead.values()),
                       default=0),
        n_dead_max_cell=max(dead, key=lambda k: dead[k]['n_dead'])
        if dead else None,
        components=rows)
    os.makedirs('outputs', exist_ok=True)
    json.dump(out, open('outputs/dye_panel_audit.json', 'w'),
              indent=1, ensure_ascii=False)

    print(f'\n  표시된 비염료 {len(marked)}개: {", ".join(out["marked_nondye"])}')
    print(f'  가장 큰 비: {top["name"]} {top["ratio_max"]:.1f} '
          f'(표시 {"있음" if top["note_kind"] == "nondye" else "없음"})')
    print(f'  표시된 형광체의 최대 비: {marked_fluor_max:.1f}')
    if not out['highest_ratio']['marked']:
        print('  [주의] 가장 큰 비를 내는 성분이 비염료로 표시돼 있지 않다.')
    print(f'  셀별 채널: 정확히 0 / 0 은 아니고 최대 < {DEAD_CUT}')
    for k, v in dead.items():
        print(f"    {k:<12} 0인열 {v['n_zero']}  잡음열 {v['n_dead']}  "
              f"{', '.join(v['dead'])}")
    print('\n[done] outputs/dye_panel_audit.json')


if __name__ == '__main__':
    main()
