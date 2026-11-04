#!/usr/bin/env python3
"""네 패널 구성을 하나의 검정 family 로 묶어 BH 보정한다 (2026-11-02 신설).

왜 신설했는가. outputs/multiplicity_across_panels.json 이 저장소에 예치되어
있고 SI 가 그것을 인용하는데, **그 파일을 만드는 코드가 저장소에 없었다.**
GATE 8 7차에서 리뷰어가 적발했다. 더 나쁜 것은 그 파일의 값이 현행 산출물과
어긋난다는 점이었다: 35종 패널 merged median_diff 가 그 파일에서는 0.15612,
outputs/permutation_statistics.json 에서는 0.16947 이다. 0.15612 는 어느 현행
산출물에도 없다. conformer_w_max 열을 거리에서 빼기 전(40/42차원) 공간에서
만들어진 잔재로 보인다.

예치물에 생성 코드가 없으면 그 수치는 출처 불명이다. 이 스크립트가 그
공백을 메운다. 모든 p 는 현행 산출물에서만 읽는다.

구성 넷 × 공간 셋 × 통계량 여섯 = 72 검정.

쓰는 법:  python3 multiplicity_panels.py
"""
import json

import numpy as np

OUT = 'outputs/multiplicity_across_panels.json'

# 구성 이름 -> (출처 파일, 그 파일 안의 경로). SI·본문이 쓰는 이름을 키로 둔다.
CONFIGS = {
    'panel28_designated3': ('outputs/permutation_statistics.json',
                            'panel28_designated'),
    'panel35_designated3': ('outputs/permutation_statistics.json',
                            'panel35_designated'),
    'panel35_designated_plus_literature8': (
        'outputs/permutation_statistics.json', 'panel35_all'),
    'pesticide10': ('outputs/pesticide_panel.json', None),
}

# 공간 이름 정규화. permutation_statistics 는 A_rdkit/B_electronic/C_merged,
# pesticide_panel 은 structural/electronic3/merged 를 쓴다. 같은 것이다.
SPACE = {'A_rdkit': 'structural', 'B_electronic': 'electronic',
         'C_merged': 'merged', 'structural': 'structural',
         'electronic3': 'electronic', 'merged': 'merged'}

STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']


def collect():
    """(config, space, statistic, p) 를 전부 모은다. 빠진 칸은 즉시 중단."""
    rows = []
    for cfg, (path, key) in CONFIGS.items():
        d = json.load(open(path))
        node = d if key is None else d[key]
        spaces = node['spaces'] if 'spaces' in node else node
        for raw_space, sd in spaces.items():
            if not isinstance(sd, dict):
                continue
            # permutation_statistics 는 statistics/<stat>/p_one_sided,
            # pesticide_panel 은 p/<stat> 로 담는다.
            if 'statistics' in sd:
                src = {k: v['p_one_sided'] for k, v in sd['statistics'].items()
                       if isinstance(v, dict) and 'p_one_sided' in v}
            elif 'p' in sd:
                src = sd['p']
            else:
                continue
            if raw_space not in SPACE:
                raise SystemExit(f'[중단] 공간 이름 {raw_space!r} 을 모른다')
            for st in STATS:
                if st not in src:
                    raise SystemExit(
                        f'[중단] {cfg}/{raw_space} 에 통계량 {st} 가 없다. '
                        f'산출물이 바뀌었으면 STATS 를 고칠 것.')
                rows.append((cfg, SPACE[raw_space], st, float(src[st])))
    return rows


def bh(p):
    """Benjamini-Hochberg 보정값. 단조 보정까지 포함한다."""
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def main():
    rows = collect()
    n_cfg = len({r[0] for r in rows})
    if len(rows) != n_cfg * 18:
        raise SystemExit(f'[중단] 검정 수가 {len(rows)} 다. '
                         f'구성 {n_cfg}개 × 18 = {n_cfg * 18} 이어야 한다.')
    q = bh([r[3] for r in rows])

    best = {}
    for (cfg, sp, st, p), qi in zip(rows, q):
        if cfg not in best or qi < best[cfg]['q'] - 1e-15:
            best[cfg] = dict(q=float(qi), space=sp, statistic=st, p=float(p))

    out = {
        'n_tests': len(rows),
        'n_configs': n_cfg,
        'generated_by': 'multiplicity_panels.py',
        'sources': {k: v[0] for k, v in CONFIGS.items()},
        'min_q_by_config': best,
    }
    json.dump(out, open(OUT, 'w'), indent=1, ensure_ascii=False)
    print(f'[multiplicity] {OUT} 생성. 검정 {len(rows)}건, 구성 {n_cfg}개')
    for cfg, b in best.items():
        print(f'  {cfg:<38} q={b["q"]:.4f}  '
              f'({b["space"]}/{b["statistic"]}, p={b["p"]:.5f})')


if __name__ == '__main__':
    main()
