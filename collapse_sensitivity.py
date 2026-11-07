#!/usr/bin/env python3
"""중복 채널 축약의 민감도 (2026-11-07 신설, 10차 R3 지적).

왜 필요한가. SI 가 "the collapse weakens the two cells it touches, from
rho = +0.568 and +0.539 to +0.518 and +0.496, and moves the smallest adjusted
value over the twelve cells from 0.067 to 0.124" 라고 적는데, 저장소의
어느 실행도 그 여섯 수를 내놓지 않는다. 축약은 deltaE_premise.load 에서
무조건 적용되므로 비축약 쪽을 계산하는 코드가 아예 없었다. 즉 원고가
적는 '접기 전' 수치는 구판의 기억이고 [C-3] 위반이다.

여기서 같은 12셀 Mantel 을 두 번 돈다. 한 번은 현행대로(축약), 한 번은
축약을 끄고. 두 결과와 그 차이를 산출물로 남기고 SI 가 그것을 인용한다.

축약을 끄는 방법. deltaE_premise.collapse_duplicate_channels 를 자리에서
바꿔 끼운다. rollback_c.deltaE_corrected 가 m.load 를 부르고 load 가 그
함수를 부르므로, 모듈 속성을 갈아 끼우면 두 경로가 모두 바뀐다.
"""
import itertools
import json
from itertools import permutations as _permutations

import numpy as np
import pandas as pd
from scipy import stats

import deltaE_premise as m
from spaces import ELECTRONIC_28, META_NUM, PROTECTED, pick

# rollback_c.py 와 같은 정의를 쓴다. omega 는 거리공간에서 제외한다.
OMEGA = 'electrophilicity_index_eV'
ELEC3 = [c for c in ELECTRONIC_28 if c != OMEGA]

SEED = 0xC0DE   # rollback_c.py 와 동일


def _zdist(df, cols):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Xz = (X - X.mean(0)) / X.std(0, ddof=1)
    return np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1)), keep


def _bh_min(pv):
    pv = np.sort(np.asarray(pv, float))
    k = len(pv)
    return float(np.minimum.accumulate(
        (pv * k / np.arange(1, k + 1))[::-1])[::-1][0])


def twelve_cells():
    """현행 상태(축약 여부는 호출자가 정한다)에서 12셀 Mantel 을 돈다."""
    df = pd.read_csv('descriptor_cache_xtb.csv')
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == 28 and c not in META_NUM]
    _struct = [c for c in num if c not in ELECTRONIC_28]
    _elec3 = [c for c in num if c in ELEC3]
    SPACES = {'structural': _struct, 'electronic3': _elec3,
              'merged': _struct + _elec3}
    idx = {nm: i for i, nm in enumerate(df.compound_name)}
    DIST = {sp: _zdist(df, cols)[0] for sp, cols in SPACES.items()}

    cells, nchan = [], {}
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                D, dy = m.load(tag, mode)
            except Exception as e:
                print(f'  [수집 불가: {tag} {mode} {e}]')
                continue
            nchan[f'{tag}.{mode}'] = len(dy)
            P = m.consensus(m.centred(m.raw_patterns(D, dy)))
            ags = [a for a in sorted(P) if a in m.MAP]
            if len(ags) < 4:
                continue
            k = len(ags)
            pr = list(itertools.combinations(range(k), 2))
            dis = np.array([1 - float(np.dot(P[ags[i]], P[ags[j]]))
                            for i, j in pr])
            for sp in SPACES:
                M = np.array([[DIST[sp][idx[m.MAP[a]], idx[m.MAP[b]]]
                               for b in ags] for a in ags])
                d = np.array([M[i, j] for i, j in pr])
                r = float(stats.spearmanr(d, dis).statistic)
                cnt = 0
                perms = list(_permutations(range(k)))
                for pmx in perms:
                    dp = np.array([M[pmx[i], pmx[j]] for i, j in pr])
                    rp = stats.spearmanr(dp, dis).statistic
                    if rp == rp and abs(rp) >= abs(r):
                        cnt += 1
                cells.append(dict(dataset=tag, mode=mode, space=sp,
                                  n_compounds=k, rho=r,
                                  p_mantel=cnt / len(perms)))
    return cells, nchan


def main():
    print('=' * 78)
    print('[축약 민감도] 중복 채널을 접은 경우와 접지 않은 경우')
    print('=' * 78)

    on_cells, on_chan = twelve_cells()

    # 축약을 끈다. 선언한 열 전부를 그대로 돌려준다.
    _orig = m.collapse_duplicate_channels
    m.collapse_duplicate_channels = lambda D, dy, key: list(dy)
    try:
        off_cells, off_chan = twelve_cells()
    finally:
        m.collapse_duplicate_channels = _orig

    key = lambda c: (c['dataset'], c['mode'], c['space'])
    on = {key(c): c for c in on_cells}
    off = {key(c): c for c in off_cells}
    shared = sorted(set(on) & set(off))

    # 축약이 건드리는 셀 = 유효 채널 수가 달라진 plate-channel
    touched = sorted(k for k in on_chan if on_chan[k] != off_chan.get(k))
    print(f'\n  축약이 채널 수를 바꾸는 plate-channel: {len(touched)}개')
    for t in touched:
        print(f'    {t:<12} 선언 {off_chan[t]} -> 유효 {on_chan[t]}')

    print('\n  rho 가 가장 많이 움직인 셀 (|차이| 내림차순):')
    diffs = sorted(shared, key=lambda k: -abs(on[k]['rho'] - off[k]['rho']))
    for k in diffs[:6]:
        print(f'    {k[0]} {k[1]:<7} {k[2]:<12} '
              f'rho {off[k]["rho"]:+.4f} -> {on[k]["rho"]:+.4f}   '
              f'p {off[k]["p_mantel"]:.4f} -> {on[k]["p_mantel"]:.4f}')

    q_on = _bh_min([on[k]['p_mantel'] for k in sorted(on)])
    q_off = _bh_min([off[k]['p_mantel'] for k in sorted(off)])
    print(f'\n  12셀 BH 최소 q: 비축약 {q_off:.4f} -> 축약 {q_on:.4f}')

    # SI 가 인용하는 두 셀: 비축약에서 rho 가 가장 큰 두 셀
    top = sorted(shared, key=lambda k: -off[k]['rho'])[:2]
    out = dict(
        generated_by='collapse_sensitivity.py',
        seed=SEED,
        touched_plate_channels={t: dict(declared=off_chan[t],
                                        effective=on_chan[t])
                                for t in touched},
        bh_min_q_uncollapsed=q_off,
        bh_min_q_collapsed=q_on,
        n_cells=len(shared),
        top_two_uncollapsed=[dict(dataset=k[0], mode=k[1], space=k[2],
                                  rho_uncollapsed=off[k]['rho'],
                                  rho_collapsed=on[k]['rho'],
                                  p_uncollapsed=off[k]['p_mantel'],
                                  p_collapsed=on[k]['p_mantel'])
                             for k in top],
        cells=[dict(dataset=k[0], mode=k[1], space=k[2],
                    rho_uncollapsed=off[k]['rho'], rho_collapsed=on[k]['rho'],
                    p_uncollapsed=off[k]['p_mantel'],
                    p_collapsed=on[k]['p_mantel'])
               for k in sorted(shared)])
    json.dump(out, open('outputs/collapse_sensitivity.json', 'w'),
              indent=1, ensure_ascii=False)
    print('\n[done] outputs/collapse_sensitivity.json')
    print('  비축약 최강 두 셀:')
    for t in out['top_two_uncollapsed']:
        print(f"    {t['dataset']} {t['mode']} {t['space']}  "
              f"rho {t['rho_uncollapsed']:+.4f} -> {t['rho_collapsed']:+.4f}")


if __name__ == '__main__':
    main()
