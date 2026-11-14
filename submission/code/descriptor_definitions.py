#!/usr/bin/env python3
"""서술자 정의가 특정 화합물에 주는 값이 화학적으로 틀린 자리 (2026-11-14 신설).

왜 필요한가. 외부 검토에서 두 자리가 지적됐고 둘 다 사실이다. 이 파일은 그
사실과 그것이 결론을 얼마나 움직이는지를 수치로 예치한다.

(1) 물의 수소결합 공여자·수용자 수가 0 이다. RDKit 의 Lipinski 계열 정의가
    `[O;H1]` 만 공여자로 보고, 수용자 쪽도 물을 세지 않는다. SI 가 보호집합을
    "hydrogen-bond acceptor count" 라고 부르므로, 물에 0 을 주는 값은 그
    이름대로 읽으면 화학적으로 거짓이다. 이 값은 외형상의 문제로 끝나지
    않는다. 물이 비공여자로 기록되면 n_OH 열이 HBD 와 28종 패널에서 완전히
    같아지고(둘 다 메탄올·에탄올·DPGME 에서만 1), |r| = 1.00 이므로 상관필터가
    n_OH 를 삭제한다. 물을 고치면 그 우연한 동일성이 깨져 좌표가 하나 늘고
    헤드라인 p 가 움직인다.

(2) 선형 분자의 InertialShapeFactor 가 0 이다. ISF = PMI2/(PMI1*PMI3) 이고
    사이안화수소·염화사이안은 PMI1 이 수치영이므로 이 양은 발산한다. RDKit 은
    영 나눗셈을 0 으로 돌려주는데, 0 은 발산의 정반대 극한이다. 같은 준선형
    간섭물인 아세토니트릴은 0.316 을 받으므로, 준선형 세 종을 한 끝에 모아야
    하는 축이 그 셋을 갈라놓고 부호도 뒤집혀 있다.

무엇을 하고 무엇을 하지 않는가. 캐시를 손으로 고치지 않는다. 캐시는
build_descriptors.py 가 생성하고 관문이 53열 전량을 재현하는 대상이며, 그
재현성이 이 논문의 예치물 중 가장 단단한 부분이다. 손으로 고친 값은 생성기가
재현하지 못한다. 그래서 정의를 바꾼 대안 구성을 여기서 따로 계산해 보고하고,
본문·SI 는 채택이 아니라 선언된 민감도로 싣는다. acetonitrile_panel.py 와
niosh_panel.py 가 같은 형태다.

[C-3] 모든 값은 아래 코드가 실제로 돌려 얻은 것이다.
"""
import itertools
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick

PATH = 'descriptor_cache_xtb.csv'
ISF = 'InertialShapeFactor'


def spaces_of(df):
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    keep, _ = pick(df, [c for c in num if c != OMEGA_DERIVED], PROTECTED)
    e3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]
    return [c for c in keep if c not in e3], keep


def dmat(df, cols):
    X = df[cols].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    return D


def exact(D, kl):
    nn = float(np.median(D.min(1)))
    ir = [i for i, k in enumerate(kl) if k == 'real_agent']
    ip = [i for i, k in enumerate(kl) if k != 'real_agent']
    pos = {g: j for j, g in enumerate(ip)}
    sim = [i for i, k in enumerate(kl) if k == 'simulant']
    true = tuple(sorted(pos[s] for s in sim))
    ALL = np.asarray(list(itertools.combinations(range(len(ip)), 3)), np.int32)
    Ds = D[np.ix_(ir, ip)]
    npl = len(ip)
    md = np.empty(len(ALL))
    nv = np.empty(len(ALL))
    for a0 in range(0, len(ALL), 4000):
        cb = ALL[a0:a0 + 4000]
        b = len(cb)
        m = np.zeros((b, npl), bool)
        m[np.arange(b)[:, None], cb] = True
        A = Ds[:, None, :]
        ds = np.where(m[None], A, np.inf).min(2)
        di = np.where(~m[None], A, np.inf).min(2)
        md[a0:a0 + b] = (ds - di).mean(0)
        nv[a0:a0 + b] = -((ds / nn <= 2.0) & (ds <= di)).sum(0)
    ti = ALL.tolist().index(list(true))
    return float((md <= md[ti]).mean()), int(-nv[ti])


def one(df, drop_isf=False):
    d = df.drop(columns=[ISF]) if drop_isf else df
    st, mg = spaces_of(d)
    ps, nv = exact(dmat(d, st), d.compound_class.tolist())
    pm, _ = exact(dmat(d, mg), d.compound_class.tolist())
    return dict(n_structural=len(st), n_merged=len(mg),
                p_structural=ps, p_merged=pm, n_valid=nv,
                n_OH_retained=('n_OH' in st))


def main():
    base = pd.read_csv(PATH)
    w = base.index[base.compound_name == 'water'][0]
    out = {'generated_by': 'descriptor_definitions.py',
           'status': 'reported rather than adopted',
           'cache_unmodified': True,
           'derived_columns_kept_consistent': ['lser_A_proxy = HBD',
                                               'lser_B_proxy = HBA'],
           'observed': {
               'water_HBA': int(base.at[w, 'HBA']),
               'water_HBD': int(base.at[w, 'HBD']),
               'water_TPSA': float(base.at[w, 'TPSA']),
               'n_OH_equals_HBD_on_panel28': bool((base.n_OH == base.HBD).all()),
               'isf_zero_compounds': base.loc[base[ISF] == 0,
                                              'compound_name'].tolist(),
               'isf_pmi1_of_those': base.loc[base[ISF] == 0,
                                             'PMI1'].tolist(),
               'isf_acetonitrile': float(
                   base.loc[base.compound_name == 'acetonitrile', ISF].iloc[0]),
           },
           'configurations': {}}

    out['configurations']['as_published'] = one(base)

    # 파생열을 함께 고친다. build_descriptors.py 가 lser_A_proxy = HBD,
    # lser_B_proxy = HBA 로 정의하므로, HBA/HBD 만 고치면 그 둘이 정의와
    # 어긋나 상관·종속 필터의 판단이 달라진다. 실제로 그렇게 했을 때
    # lser_A_proxy 가 인공적으로 살아남아 좌표가 하나 늘었다. 정의를 고치는
    # 민감도는 정의에 묶인 열을 전부 고쳐야 뜻이 있다.
    for hba, hbd, tag in ((1, 1, 'water_HBA1_HBD1'), (1, 2, 'water_HBA1_HBD2')):
        d = base.copy()
        d.at[w, 'HBA'] = hba
        d.at[w, 'HBD'] = hbd
        d.at[w, 'lser_B_proxy'] = hba          # 정의: lser_B_proxy = HBA
        d.at[w, 'lser_A_proxy'] = hbd          # 정의: lser_A_proxy = HBD
        out['configurations'][tag] = one(d)

    out['configurations']['isf_dropped'] = one(base, drop_isf=True)

    d = base.copy()
    d.at[w, 'HBA'] = 1
    d.at[w, 'HBD'] = 1
    d.at[w, 'lser_B_proxy'] = 1
    d.at[w, 'lser_A_proxy'] = 1
    out['configurations']['both_corrections'] = one(d, drop_isf=True)

    with open('outputs/descriptor_definitions.json', 'w',
              encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    print(f"[정의] 물 HBA={out['observed']['water_HBA']} "
          f"HBD={out['observed']['water_HBD']} (TPSA "
          f"{out['observed']['water_TPSA']})")
    print(f"[정의] ISF=0: {out['observed']['isf_zero_compounds']} "
          f"(아세토니트릴 {out['observed']['isf_acetonitrile']:.3f})")
    print(f"{'구성':<22}{'구조좌표':>8}{'구조 p':>10}{'병합 p':>10}"
          f"{'매칭':>6}{'n_OH':>6}")
    for k, v in out['configurations'].items():
        print(f"{k:<22}{v['n_structural']:>8}{v['p_structural']:>10.4f}"
              f"{v['p_merged']:>10.4f}{v['n_valid']:>6}"
              f"{'Y' if v['n_OH_retained'] else 'n':>6}")
    print('  -> outputs/descriptor_definitions.json')


if __name__ == '__main__':
    main()
