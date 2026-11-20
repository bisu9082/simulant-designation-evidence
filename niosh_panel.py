#!/usr/bin/env python3
"""NIOSH 후보목록의 나머지 네 화합물을 비교군에 넣었을 때 (2026-11-14 신설).

왜 필요한가. 본 논문이 DMMP 의 지정 근거로 쓰는 NIOSH 호흡보호구 문서는
투과시험 후보를 열한 종 적고, 그 목록에 DMMP 와 함께 DIMP·DEEP·TEP·
di-n-butyl sulfide 가 있다. SI meth:niosh 는 그 사실을 적어 두고, Table S6
은 DIMP 를 NRC MIST 의 GB 지정 모사체로 등재해 두고, 그 넷이 왜 우리
패널에 없는지는 적지 않는다. 즉 기록되지 않은 패널 구성 선택이 하나 남아
있다. 이 논문의 명제가 바로 '기록되지 않은 선택이 답을 정한다' 이므로,
그 선택을 열어 값을 적는다.

구성. 패널 그대로, 비교군에 네 종을 추가한다. 지정은 현행 3종 불변이다.
이 넷은 작용제도 지정 모사체도 아니고, 같은 문서가 같은 목적으로 후보로
적은 비교화합물이다. 그 자리가 비교군이다.

[C-3] SMILES 는 IUPAC 명칭에서 직접 적었고, 서술자는 build_descriptors.py
로 생성했다. 같은 생성기가 기존 28종의 RDKit 수치열 53개를 불일치 0 으로
재현한다 (--verify 통과). 추정값은 쓰지 않았다.

범위 한정. 구조공간만 확장할 수 있다. 전자 좌표 세 개는 xtb 반경험 계산이
필요하고 그 바이너리는 이 저장소에 없으며 기존 캐시를 만든 버전도 기록돼
있지 않다 (본문 Reproducibility 가 공시한 그 한계다). 다른 버전으로 네
종만 계산하면 한 거리공간 안에 두 바이너리의 출력이 섞인다. 그래서 전자·
병합 공간은 확장하지 않고, 그 사실을 보고한다.

채택이 아니라 보고다. 헤드라인 패널은 28종·35종 그대로 두고 이 값을 민감도
분석으로 싣는다. acetonitrile_panel.py 와 pesticide_panel.py 가 같은 형태다.
"""
import itertools
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

from perm_statistics import STATS
from pesticide_panel import bh_min, dmat, spaces_of

# NIOSH 후보목록 11종 중 우리 패널에 없는 열 종. 2026-11-14 전면 교정.
#
# 종전 이 목록은 네 종이었다. 1차 출처(원문 PDF)를 다시 열어 11종을 축자
# 대조한 결과 패널에 있는 것은 DMMP 하나이고 나머지 열이 빠져 있었다.
# 누락 기전은 약어 충돌이다. NIOSH 목록의 'DCP' 는 2,4-dichlorophenol 인데
# 이 논문 35종 패널의 'DCP' 는 diethyl chlorophosphate 다. 약어로 대조하면
# "이미 패널에 있다" 로 처리된다. 나머지 다섯은 황화물·이염화물·에터·DEMP 로
# 명칭이 패널의 어느 것과도 겹치지 않는데 함께 빠졌다.
#
# 원문 표기 순서 그대로 적는다 (원문의 'Triethly' 오타는 바로잡아 적었다).
NIOSH_CANDIDATES = [
    ('2-chloroethyl phenyl sulfide', 'ClCCSc1ccccc1',
     '2-chloroethyl phenyl sulfide'),
    ('2-chloroethyl cyclohexyl sulfide', 'ClCCSC1CCCCC1',
     '2-chloroethyl cyclohexyl sulfide'),
    ('di-n-butyl sulfide', 'CCCCSCCCC',
     'di-n-butyl sulfide'),
    ('1,6-dichlorohexane', 'ClCCCCCCCl',
     '1,6-dichlorohexane'),
    ('bis(4-chlorobutyl) ether', 'ClCCCCOCCCCCl',
     'bis(4-chlorobutyl) ether'),
    ('DEMP', 'CCOP(C)(=O)OCC',
     'diethyl methanephosphonate'),
    ('DEEP', 'CCOP(CC)(=O)OCC',
     'diethyl ethanephosphonate'),
    ('DIMP', 'CC(C)OP(C)(=O)OC(C)C',
     'diisopropyl methylphosphonate'),
    ('TEP', 'CCOP(=O)(OCC)OCC',
     'triethyl phosphate'),
    ('2,4-dichlorophenol', 'Oc1ccc(Cl)cc1Cl',
     '2,4-dichlorophenol'),
]
# 하위호환 별칭. 종전 보고(넷만 넣은 값)를 재현해 보려는 독자를 위해 남긴다.
NIOSH4 = [r for r in NIOSH_CANDIDATES
          if r[0] in ('DIMP', 'DEEP', 'TEP', 'di-n-butyl sulfide')]
PANELS = [('panel28_incumbent3', 'descriptor_cache_xtb.csv'),
          ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv')]
SOURCE = ('NIOSH/NPPTL CWA simulant programme (31 July 2002, draft for '
          'discussion), candidate list for barrier-material permeation '
          'testing: eleven compounds, of which DMMP alone is on either of '
          'our panels. See Supporting Information Section meth:niosh.')
CHUNK = 4000


def build_four():
    """build_descriptors.py 로 네 종의 RDKit 서술자를 만든다."""
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, 'niosh4.csv')
        out = os.path.join(td, 'niosh4_desc.csv')
        pd.DataFrame([{'compound_name': n, 'compound_class': 'interferent',
                       'agent_series': '', 'input_smiles': s}
                      for n, s, _ in NIOSH_CANDIDATES]).to_csv(src, index=False)
        r = subprocess.run([sys.executable, 'build_descriptors.py',
                            '--new', src, '--out', out],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f'[niosh] 서술자 생성 실패\n{r.stderr}')
        return pd.read_csv(out)


def exact_perm(D, kl, sim_idx, names):
    """구조공간 전수 순열. 통계량·귀무·매칭규칙은 본문과 동일하다."""
    D = D.copy().astype(float)
    np.fill_diagonal(D, np.inf)
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
    ti = ALL.index(true_combo)
    p = {s: float((null[s] <= null[s][ti]).sum() / len(ALL)) for s in STATS}
    # 작용제별 최근접 지정모사체 / 최근접 비지정과 그 화합물
    m = np.zeros(npool, bool)
    m[list(true_combo)] = True
    ds = np.where(m, Ds, np.inf).min(1)
    di = np.where(~m, Ds, np.inf).min(1)
    blk = np.where(~m, Ds, np.inf)
    agents = [dict(agent=names[i_real[r]], d_sim=float(ds[r]),
                   d_blocker=float(di[r]),
                   blocker=names[i_pool[int(np.argmin(blk[r]))]],
                   matched=bool(ds[r] / nn_med <= 2.0 and ds[r] <= di[r]))
              for r in range(len(i_real))]
    print(f'    {len(ALL):,} 지정 전수, {time.time() - t0:.1f}s')
    # 2026-11-18. k 가 바뀌면 귀무 자체가 이동한다는 것을 수치로 남긴다.
    # mean_diff 귀무 평균과 전체 분포를 함께 돌려준다 (교차 판독용).
    return dict(p=p, n_perm=len(ALL), nn_median=nn_med,
                null_mean={s: float(null[s].mean()) for s in STATS},
                _null_mean_diff=null['mean_diff'],
                n_valid=int(-null['n_valid'][ti]),
                observed={s: float(null[s][ti]) for s in STATS},
                agents=agents)


def one_panel(label, path, four):
    base = pd.read_csv(path)
    names0 = base.compound_name.tolist()
    dup = [n for n, _, _ in NIOSH_CANDIDATES if n in names0]
    if dup:
        sys.exit(f'{path}: {dup} 가 이미 패널에 있다')
    cols0 = spaces_of(base)['structural']
    ext = pd.concat([base, four], ignore_index=True)
    cols1 = spaces_of(ext)['structural']
    rec = {}
    _nulls = {}
    # 2026-11-18 (GATE 8 후속, 리뷰어 D·B 지적). 네 번째 구성을 추가한다.
    #
    # 왜. 복원한 열 종 가운데 DIMP 는 이 논문 자신의 지정표
    # (sources/designations_20261106.csv) 가 NRC 1997 MIST 3장의 GB 직물투과
    # 지정 모사체로 기록한 화합물이다. DMMP(VX)·methyl salicylate(HD) 와
    # 같은 문서, 같은 측정 범주, 같은 근거 서술이다. 그런데 위 세 구성은
    # DIMP 를 비지정 경쟁자로 세우고, 그 상태에서 DIMP 가 소만·사이클로사린·
    # VX 의 매칭을 깬다. 지정집합은 이 논문이 선언을 요구하는 일곱 자유도
    # 가운데 하나이므로, 선언만 하고 민감도를 돌리지 않을 수 없다.
    #
    # k 가 3 에서 4 로 커지면 귀무도 함께 바뀐다. k 개 추출의 최소값은 k 가
    # 커질수록 확률적으로 작아지므로 k=3 귀무는 k=4 관측의 기준분포가
    # 아니다. 이것이 2.6 절이 전수열거를 도입한 군크기 교란 그 자체다.
    # exact_perm 은 sim_idx 길이로 k 를 정하므로 재열거가 자동으로 된다.
    DIMP_DESIG = 'plus_all_dimp_designated'
    for tag, df, cols in (('baseline', base, cols0),
                          ('plus_all_recomputed', ext, cols1),
                          ('plus_all_columns_fixed', ext, cols0),
                          (DIMP_DESIG, ext, cols1)):
        kl = df.compound_class.tolist()
        nm = df.compound_name.tolist()
        sim_idx = [i for i, k in enumerate(kl) if k == 'simulant']
        if tag == DIMP_DESIG:
            if 'DIMP' not in nm:
                sys.exit('DIMP 가 확장 프레임에 없다')
            sim_idx = sorted(sim_idx + [nm.index('DIMP')])
        print(f'  {label}/{tag}: {len(df)}종, 좌표 {len(cols)}개, '
              f'지정 {len(sim_idx)}종')
        r = exact_perm(dmat(df, cols), kl, sim_idx, nm)
        r['n_compounds'] = len(df)
        r['n_dims'] = len(cols)
        r['designated'] = [nm[i] for i in sim_idx]
        _nulls[tag] = r.pop('_null_mean_diff')
        rec[tag] = r
        print(f'      사전지정 mean_diff p = {r["p"]["mean_diff"]:.5f}  '
              f'유효매칭 {r["n_valid"]}/11  '
              f'6통계량 BH 최소 q {bh_min(list(r["p"].values())):.4f}')
    # k=3 귀무에 k=4 관측을 대고 읽으면 어떤 값이 나오는가. 이 수치는
    # 채택값이 아니라 오독의 크기를 적기 위한 것이다.
    if DIMP_DESIG in rec and 'plus_all_recomputed' in rec:
        _obs4 = rec[DIMP_DESIG]['observed']['mean_diff']
        _n3 = _nulls['plus_all_recomputed']
        rec[DIMP_DESIG]['misread_against_k3_null'] = {
            'note': ('k=3 귀무에 k=4 관측을 대고 읽은 값. 군크기가 다르면 '
                     'k 개 추출의 최소값 분포가 달라지므로 기준분포가 아니다. '
                     '보고값이 아니라 오독의 크기를 적는다.'),
            'observed_k4': float(_obs4),
            'null_mean_k3': float(_n3.mean()),
            'null_mean_k4': float(rec[DIMP_DESIG]['null_mean']['mean_diff']),
            'p_if_read_against_k3': float((_n3 <= _obs4).sum() / len(_n3)),
            'p_correct_k4': float(rec[DIMP_DESIG]['p']['mean_diff']),
        }
    rec['columns_only_in_baseline'] = sorted(set(cols0) - set(cols1))
    rec['columns_only_in_extended'] = sorted(set(cols1) - set(cols0))
    return rec


def main():
    four = build_four()
    print(f'[niosh] 서술자 생성 {len(four)}종 '
          f'({", ".join(four.compound_name)})')
    res = {'generated_by': 'niosh_panel.py',
           'added_to_comparison_pool': [n for n, _, _ in NIOSH_CANDIDATES],
           'smiles': {n: s for n, s, _ in NIOSH_CANDIDATES},
           'iupac_like': {n: i for n, _, i in NIOSH_CANDIDATES},
           'authority': SOURCE,
           'designated_set': ('incumbent three in the baseline and the two '
                              'all-added runs; a fourth run designates DIMP as '
                              'well — see each run\'s own designated field'),
           'spaces_extended': ['structural'],
           'spaces_not_extended': ['electronic3', 'merged'],
           'why_not_extended': ('전자 좌표는 xtb 반경험 계산이 필요하고 그 '
                                '바이너리가 저장소에 없으며 기존 캐시의 '
                                '생성 버전도 기록돼 있지 않다. 다른 버전으로 '
                                '계산하면 한 거리공간에 두 바이너리의 출력이 '
                                '섞인다.'),
           'n_candidates_in_source': 11,
           'n_on_our_panels': 1,
           'n_added': len(NIOSH_CANDIDATES),
           'abbreviation_collision': ("NIOSH 'DCP' = 2,4-dichlorophenol; "
                                      "this paper's panel 'DCP' = diethyl "
                                      "chlorophosphate. The earlier version "
                                      "of this analysis added only four "
                                      "compounds because the match was made "
                                      "on abbreviations."),
           'status': 'reported rather than adopted',
           'panels': {}}
    for label, path in PANELS:
        res['panels'][label] = one_panel(label, path, four)
    with open('outputs/niosh_panel.json', 'w', encoding='utf-8') as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print()
    for k, v in res['panels'].items():
        b = v['baseline']['p']['mean_diff']
        e = v['plus_all_recomputed']['p']['mean_diff']
        print(f"  {k}: 사전지정 p {b:.5f} -> {e:.5f}  "
              f"유효매칭 {v['baseline']['n_valid']} -> "
              f"{v['plus_all_recomputed']['n_valid']}")
    print('  -> outputs/niosh_panel.json')


if __name__ == '__main__':
    main()
