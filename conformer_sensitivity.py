#!/usr/bin/env python3
"""3D 서술자의 형태이성질체 의존성 대조 (2026-10-25 신설).

왜. build_descriptors.py --verify 가 52개 수치열 중 11개를 재현 불가로
판정한다. 그중 9개(Asphericity, Eccentricity, InertialShapeFactor, NPR1,
NPR2, PMI1, PMI2, RadiusOfGyration, SpherocityIndex)가 두 패널의 구조공간에
그대로 들어간다. 원인은 best_conf() 가 ETKDG 앙상블에서 MMFF 최저에너지
형태이성질체 하나를 골라 거기서 Descriptors3D 를 읽기 때문이다. 앙상블
평균이 아니므로 임베딩 시드가 바뀌면 값이 바뀐다.

여기서 재는 것. 시드만 바꿔 9개 3D 열을 다시 계산하고 나머지 열은 배포
캐시를 그대로 둔 뒤, 본문과 같은 pick/표준화/매칭 규칙으로 28종 패널의
유효매칭과 35종 패널의 유효매칭, 그리고 구조·병합공간 전수 순열 p 를 낸다.
본문이 '4건'이라고 적은 28종 매칭이 시드에 얼마나 민감한지 보는 것이 목적이다.

출력: outputs/conformer_sensitivity.json
"""
import json
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D

from spaces import pick, PROTECTED
from match_by_panel import merged_space

RDLogger.DisableLog('rdApp.*')

THREE_D = ['Asphericity', 'Eccentricity', 'InertialShapeFactor', 'NPR1',
           'NPR2', 'PMI1', 'PMI2', 'PMI3', 'RadiusOfGyration',
           'SpherocityIndex']
NMAX = 20
SEEDS = [0xf00d, 1, 2, 7, 1234, 99991]


RT_KCAL = 0.5922


def conf_descriptors(smiles, seed, method='ensemble'):
    """build_descriptors.py 의 임베딩 절차를 시드만 바꿔 되풀이한다.

    method='single' 은 구판(MMFF 최저에너지 형태이성질체 하나),
    method='ensemble' 은 현행(Boltzmann 가중 평균)이다.
    """
    m = Chem.MolFromSmiles(smiles)
    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.pruneRmsThresh = 0.5
    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=NMAX, params=ps))
    if not cids:
        ps.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=NMAX, params=ps))
    if AllChem.MMFFHasAllMoleculeParams(mh):
        AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
        props = AllChem.MMFFGetMoleculeProperties(mh)
        en = np.array([AllChem.MMFFGetMoleculeForceField(
            mh, props, confId=c).CalcEnergy() for c in cids], float)
    else:
        AllChem.UFFOptimizeMoleculeConfs(mh, maxIters=2000)
        en = np.array([AllChem.UFFGetMoleculeForceField(
            mh, confId=c).CalcEnergy() for c in cids], float)
    if method == 'single':
        best = cids[int(np.argmin(en))]
        return {n: float(getattr(Descriptors3D, n)(mh, confId=best))
                for n in THREE_D}
    w = np.exp(-(en - en.min()) / RT_KCAL)
    w = w / w.sum()
    out = {}
    for n in THREE_D:
        v = np.array([getattr(Descriptors3D, n)(mh, confId=c)
                      for c in cids], float)
        out[n] = float((w * v).sum())
    return out


def matches(df, cols, sim_cls=('simulant',)):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nn = float(np.median(D.min(1)))
    kl = df.compound_class.tolist()
    nm = df.compound_name.tolist()
    ia = [i for i, k in enumerate(kl) if k == 'real_agent']
    ism = [i for i, k in enumerate(kl) if k in sim_cls]
    iot = [i for i, k in enumerate(kl)
           if k != 'real_agent' and k not in sim_cls]
    rows = {}
    n = 0
    for i in ia:
        js = min(ism, key=lambda j: D[i, j])
        dn = D[i, js] / nn
        ok = bool(dn <= 2.0 and D[i, js] < min(D[i, j] for j in iot))
        rows[nm[i]] = round(float(dn), 4)
        n += ok
    return n, rows, len(keep), D


def perm_p(df, D):
    kl = df.compound_class.to_numpy()
    ir = np.flatnonzero(kl == 'real_agent')
    ip = np.flatnonzero(kl != 'real_agent')
    isim = set(np.flatnonzero(kl == 'simulant').tolist())
    true = tuple(sorted(list(ip).index(j) for j in isim))
    Ds = D[np.ix_(ir, ip)]
    allc = list(combinations(range(len(ip)), 3))
    diff = np.empty(len(allc))
    ti = None
    for t, c in enumerate(allc):
        m = np.zeros(len(ip), bool)
        m[list(c)] = True
        diff[t] = (Ds[:, m].min(1) - Ds[:, ~m].min(1)).mean()
        if c == true:
            ti = t
    return float((diff <= diff[ti]).mean())


def main():
    global METHOD
    METHOD = 'single' if '--single' in sys.argv else 'ensemble'
    print(f'[민감도] 방식 {METHOD}, 시드 {len(SEEDS)}개')
    out = {'method': METHOD, 'seeds': SEEDS, 'n_confs': NMAX,
           'cols_3d': THREE_D, 'runs': []}
    base = {}
    for tag, path in (('panel28', 'descriptor_cache_xtb.csv'),
                      ('panel35', 'descriptor_cache_xtb_35.csv')):
        base[tag] = pd.read_csv(path)
    for seed in SEEDS:
        rec = {'seed': seed}
        for tag, df0 in base.items():
            df = df0.copy()
            new = [conf_descriptors(s, seed, METHOD)
                   for s in df.input_smiles]
            for c in THREE_D:
                if c in df.columns:
                    df[c] = [r[c] for r in new]
            cls = ('simulant',)
            nm_, rows, nk, D = matches(df, merged_space(df), cls)
            rec[tag] = {'n_match': nm_, 'n_dims': nk, 'd_norm': rows}
            if tag == 'panel28':
                rec[tag]['perm_p_merged'] = perm_p(df, D)
        out['runs'].append(rec)
        t28 = rec['panel28']
        print(f"  시드 {seed:>6}  28종 매칭 {t28['n_match']}/11  "
              f"tabun {t28['d_norm'].get('tabun (GA)')}  "
              f"p {t28['perm_p_merged']:.4f}   "
              f"35종 매칭 {rec['panel35']['n_match']}/11", flush=True)
    fn = f'outputs/conformer_sensitivity_{METHOD}.json'
    json.dump(out, open(fn, 'w'), indent=1, ensure_ascii=False)
    print(f'\n[done] {fn}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
