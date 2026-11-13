#!/usr/bin/env python3
"""아세토니트릴을 모사체로 라벨한 구성 (2026-11-13 신설).

왜 필요한가. 이 논문의 기존 3종 중 둘 -- 디클로로메탄과 DPGME -- 은
연구문헌 한 편(lee2011sno2)에서 왔다. 그 논문은 네 화합물을 쓰고 한
문장으로 분류한다:

  "chemical agent simulants of blood agents, nerve agents, blister agents,
   and pulmonary agents, respectively"

순서대로 아세토니트릴=혈액, DMMP=신경, DPGME=수포, 디클로로메탄=질식이다.
우리 패널은 그 넷 중 셋을 모사체로 받고 아세토니트릴만 간섭물로 두었다.
그 선택을 기록한 적이 없다.

그런데 아세토니트릴은 결과를 지배한다. threshold_sweep.json 에서
전자공간의 GB·GD·GF·CK 차단자가 아세토니트릴이고, 본문도 "acetonitrile
sits closer to sarin, soman and cyclosarin than DMMP does" 라고 적는다.
즉 같은 출처가 모사체로 지정한 화합물이, 그 모사체들을 떨어뜨리는
비교화합물로 쓰이고 있다.

이 논문의 명제는 "분석자가 인용을 손에 들고 옮길 수 있는 라벨이 답을
정한다" 이다. 여기 옮길 수 있는 라벨이 하나 있고, 우리가 옮기지 않았다.
그래서 옮겨 보고 값을 적는다. 정의는 pesticide_panel.py 를 그대로 따른다.

[C-3] 분류 근거는 lee2011sno2 원문을 직접 읽어 확인했다. 추정하지 않았다.
"""
import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

from perm_statistics import STATS
from pesticide_panel import spaces_of, dmat, bh_min, CHUNK

PANELS = [('panel35_incumbent3_plus_acetonitrile',
           'descriptor_cache_xtb_35.csv'),
          ('panel28_incumbent3_plus_acetonitrile',
           'descriptor_cache_xtb.csv')]
EXTRA = 'acetonitrile'
SOURCE = ('lee2011sno2, Methods: "chemical agent simulants of blood agents, '
          'nerve agents, blister agents, and pulmonary agents, respectively" '
          '-- acetonitrile, DMMP, DPGME, dichloromethane in that order.')


def exact_perm(D, kl, sim_idx, label):
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in sim_idx))
    k = len(true_combo)
    ALL = list(itertools.combinations(range(len(i_pool)), k))
    C = np.asarray(ALL, dtype=np.int32)
    Ds = D[np.ix_(i_real, i_pool)]
    npool = len(i_pool)
    null = {s: np.empty(len(ALL)) for s in STATS}
    t0 = time.time()
    for a0 in range(0, len(C), CHUNK):
        cb = C[a0:a0 + CHUNK]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
        null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
        null['mean_logratio'][a0:a0 + b] = np.log(ds / di).mean(0)
        null['mean_dsim'][a0:a0 + b] = ds.mean(0)
        null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
        null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    el = time.time() - t0
    print(f'    [{label}] {len(ALL):,} 지정 전수, {el:.1f}s')
    ti = ALL.index(true_combo)
    p = {s: float((null[s] <= null[s][ti]).sum() / len(ALL)) for s in STATS}
    return p, len(ALL), nn_med, int(-null['n_valid'][ti])


def one_panel(label, path):
    df = pd.read_csv(path)
    kl = df.compound_class.tolist()
    names = df.compound_name.tolist()
    if EXTRA not in names:
        sys.exit(f'{path}: {EXTRA} 가 없다')
    sp = spaces_of(df)
    sim_idx = ([i for i, k in enumerate(kl) if k == 'simulant']
               + [names.index(EXTRA)])
    print(f'  {label}: {len(df)}종, 작용제 {kl.count("real_agent")}종, '
          f'지정 {len(sim_idx)}종 '
          f'({", ".join(sorted(names[i] for i in sim_idx))})')
    out = {}
    for lab, cols in sp.items():
        p, nperm, nnm, nvalid = exact_perm(dmat(df, cols), kl, sim_idx, lab)
        out[lab] = dict(p=p, n_perm=nperm, nn_median=nnm, n_valid=nvalid,
                        n_dims=len(cols))
        print(f'      {lab:<12} 유효매칭 {nvalid}/11  '
              f'최소p {min(p.values()):.4f}')
    allp = [out[s]['p'][k] for s in sp for k in STATS]
    # 본문이 "여섯 스킴" 을 말하므로 값싼 둘을 함께 낸다. 저장된 p 만으로
    # 계산되는 것은 공간별 BH 와 사전지정 3검정이다. Westfall-Young 과
    # Li-Ji m_eff 는 공유 귀무분포가 필요하므로 여기서 내지 않는다.
    by_space = min(bh_min([out[s]['p'][k] for k in STATS]) for s in sp)
    prespec = bh_min([out[s]['p']['mean_diff'] for s in sp])
    return dict(panel=label, n_simulants=len(sim_idx),
                simulant_names=sorted(names[i] for i in sim_idx),
                spaces=out, min_p=float(min(allp)),
                bh_min_q=float(bh_min(allp)),
                bh_by_space_min_q=float(by_space),
                prespecified_bh_min_q=float(prespec),
                schemes_not_computed=['westfall_young', 'liji_meff',
                                      'bh_over_72'])


def main():
    res = {'generated_by': 'acetonitrile_panel.py',
           'relabel': EXTRA,
           'authority': SOURCE,
           'note': ('기존 3종 중 둘과 같은 출처가 모사체로 지정한 화합물을 '
                    '우리는 간섭물로 두었다. 옮겨 보고 값을 적는다.'),
           'panels': {}}
    for label, path in PANELS:
        res['panels'][label] = one_panel(label, path)
    with open('outputs/acetonitrile_panel.json', 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print()
    for k, v in res['panels'].items():
        print(f"  {k}: 18검정 최소 p {v['min_p']:.6f}  "
              f"BH 최소 q {v['bh_min_q']:.4f}")
    print('  -> outputs/acetonitrile_panel.json')


if __name__ == '__main__':
    main()
