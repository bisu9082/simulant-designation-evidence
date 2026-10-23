#!/usr/bin/env python3
"""모사체에서 본 최근접 작용제의 정체 (2026-10-23 신설).

왜 필요한가. 표 3 은 작용제 -> 최근접 모사체 한 방향만 싣는다. 역방향, 즉
각 모사체가 어느 작용제에 가장 가까운가는 다른 질문이고, DPGME 의 지정
문서가 계열을 명시('a preferred G agent simulant')하므로 그 계열이 공간에서
재현되는지가 검증 대상이 된다.

점추정만으로는 쓸 수 없다. 병합공간에서 DPGME 의 최근접은 HD 1.81,
GF 1.89, GB 1.91 로 간격이 5% 안쪽이라 서술자 선택이 순위를 뒤집는다.
그래서 stat_supplement2.distance_bootstrap 과 같은 방법·시드·복제수로
'최근접 작용제가 누구인가' 자체를 재표집해 빈도로 보고한다.

출력: outputs/nearest_agent_bootstrap.json
"""
import json

import numpy as np
import pandas as pd

from spaces import META_NUM, OMEGA_DERIVED, PROTECTED, pick
from stat_supplement import zdist
from stat_supplement2 import N_BOOT, SEED

PANEL = 'descriptor_cache_xtb.csv'
SIMULANTS = ['DMMP', 'dichloromethane', 'DPGME']
G_SERIES = {'sarin (GB)', 'soman (GD)', 'tabun (GA)', 'cyclosarin (GF)'}


def main():
    df = pd.read_csv(PANEL)
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    keep, _ = pick(df, [c for c in num if c != OMEGA_DERIVED], PROTECTED)
    X = df[keep].to_numpy(float)
    nm = df.compound_name.tolist()
    kl = df.compound_class.tolist()
    ia = [i for i, k in enumerate(kl) if k == 'real_agent']
    si = {s: nm.index(s) for s in SIMULANTS}

    rng = np.random.default_rng(SEED)
    tally = {s: {nm[a]: 0 for a in ia} for s in SIMULANTS}
    gcount = {s: 0 for s in SIMULANTS}
    for b in range(N_BOOT):
        c = rng.integers(0, X.shape[1], X.shape[1])
        Xb = X[:, c]
        sd = Xb.std(0, ddof=1)
        ok = sd > 1e-12
        Zb = (Xb[:, ok] - Xb[:, ok].mean(0)) / sd[ok]
        D = zdist(Zb)
        np.fill_diagonal(D, np.inf)
        for s, i in si.items():
            a = min(ia, key=lambda x: D[i, x])
            tally[s][nm[a]] += 1
            if nm[a] in G_SERIES:
                gcount[s] += 1
        if (b + 1) % 200 == 0:
            print(f'  [nearest-agent] {b + 1}/{N_BOOT}')

    out = {'n_boot': N_BOOT, 'seed': SEED, 'panel': PANEL,
           'g_series': sorted(G_SERIES), 'simulants': {}}
    for s in SIMULANTS:
        t = {k: v for k, v in tally[s].items() if v}
        top = max(t.items(), key=lambda kv: kv[1])
        out['simulants'][s] = {
            'counts': dict(sorted(t.items(), key=lambda kv: -kv[1])),
            'top_agent': top[0],
            'top_pct': round(100 * top[1] / N_BOOT, 1),
            'g_series_pct': round(100 * gcount[s] / N_BOOT, 1)}

    with open('outputs/nearest_agent_bootstrap.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(f'[nearest-agent] outputs/nearest_agent_bootstrap.json '
          f'({N_BOOT}회, 시드 {SEED:#x})')
    for s in SIMULANTS:
        o = out['simulants'][s]
        print(f"  {s:16s} {o['top_agent']} {o['top_pct']}% / "
              f"G계열 {o['g_series_pct']}%")


if __name__ == '__main__':
    main()
