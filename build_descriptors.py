#!/usr/bin/env python3
"""RDKit 서술자 생성 (2026-10-12 재작성).

원래 이 표(descriptor_cache.csv)를 만든 스크립트가 세션 사이에 유실됐다.
패널을 확장하려면 신규 화합물에 **같은 정의**로 서술자를 붙여야 하므로
재작성한다. 재작성본을 믿기 전에 기존 28종을 재현하는지 먼저 검증한다
(--verify). 재현되지 않는 열은 그대로 보고하고, 재현되는 열만 신규
화합물에 쓴다.

컬럼 정의는 기존 표에서 역산해 확인했다.
  lser_A_proxy = HBD
  lser_B_proxy = HBA
  lser_L_proxy = MR_crippen / 10
  lser_S_proxy = (TPSA + MR_crippen) / 100
이들은 Abraham 용매화 파라미터가 아니라 TPSA/MR/HBA/HBD 의 선형 재척도다.
명칭이 오해를 부르므로 본문에서도 그렇게 서술해야 한다 (리뷰어 C Minor 2).
"""
import argparse
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Descriptors3D, rdMolDescriptors

RDLogger.DisableLog('rdApp.*')

SMA = {
    'n_P': '[P]',
    'n_PdO': '[PX4]=[OX1]',
    'n_P_F': '[PX4][F]',
    'n_P_S': '[PX4][SX2]',
    'n_P_CN': '[PX4]C#N',
    'n_P_alkoxy': '[PX4][OX2][CX4]',
    'n_P_N': '[PX4][NX3]',
    'n_P_C': '[PX4][CX4]',
    'n_As': '[As]',
    'n_S_divalent': '[SX2]',
    'n_halogen': '[F,Cl,Br,I]',
    'n_nitrile': '[CX2]#[NX1]',
    'n_nitro': '[N+](=O)[O-]',
    'n_carbonyl': '[CX3]=[OX1]',
    'n_OH': '[OX2H]',
}
PAT = {k: Chem.MolFromSmarts(v) for k, v in SMA.items()}


def leaving_group(m):
    # 'Cl' 은 2026-10-12 추가. 기존 28종에는 P-Cl 화합물이 없어 범주가 없었고,
    # 확장 패널의 DCP(diethyl chlorophosphate)가 처음으로 해당한다. 이 열은
    # 문자열 메타데이터이며 수치 서술자 집합에는 들어가지 않는다.
    for key, sma in [('F', '[PX4][F]'), ('Cl', '[PX4][Cl]'), ('CN', '[PX4]C#N'),
                     ('S-alkyl', '[PX4][SX2][#6]'), ('alkoxy', '[PX4][OX2][#6]')]:
        p = Chem.MolFromSmarts(sma)
        if p is not None and m.HasSubstructMatch(p):
            return key
    return 'none'


def embed(m, seed=0xf00d, nmax=20):
    """수소 추가 -> ETKDGv3 다중 임베딩 -> MMFF(불가 시 UFF) 최적화 -> RMSD 가지치기."""
    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.pruneRmsThresh = 0.5
    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=nmax, params=ps))
    if not cids:
        ps.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=nmax, params=ps))
    ff = 'MMFF94'
    if AllChem.MMFFHasAllMoleculeParams(mh):
        AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
    else:
        ff = 'UFF'
        AllChem.UFFOptimizeMoleculeConfs(mh, maxIters=2000)
    return mh, cids, ff


def best_conf(mh, ff):
    if ff == 'MMFF94' and AllChem.MMFFHasAllMoleculeParams(mh):
        props = AllChem.MMFFGetMoleculeProperties(mh)
        en = [AllChem.MMFFGetMoleculeForceField(mh, props, confId=c).CalcEnergy()
              for c in range(mh.GetNumConformers())]
    else:
        en = [AllChem.UFFGetMoleculeForceField(mh, confId=c).CalcEnergy()
              for c in range(mh.GetNumConformers())]
    return int(np.argmin(en))


def descriptors(name, smiles, klass='', series=''):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    d = dict(compound_name=name, compound_class=klass, agent_series=series,
             input_smiles=smiles, canonical_smiles=Chem.MolToSmiles(m),
             formula=rdMolDescriptors.CalcMolFormula(m))
    d['MW'] = Descriptors.MolWt(m)
    d['heavy_atoms'] = m.GetNumHeavyAtoms()
    d['logP_crippen'] = Crippen.MolLogP(m)
    d['MR_crippen'] = Crippen.MolMR(m)
    d['TPSA'] = rdMolDescriptors.CalcTPSA(m)
    d['HBA'] = rdMolDescriptors.CalcNumHBA(m)
    d['HBD'] = rdMolDescriptors.CalcNumHBD(m)
    d['rot_bonds'] = rdMolDescriptors.CalcNumRotatableBonds(m)
    d['n_arom_rings'] = rdMolDescriptors.CalcNumAromaticRings(m)
    d['fraction_csp3'] = rdMolDescriptors.CalcFractionCSP3(m)
    d['labute_asa'] = rdMolDescriptors.CalcLabuteASA(m)
    d['n_heteroatoms'] = rdMolDescriptors.CalcNumHeteroatoms(m)
    d['balaban_j'] = Descriptors.BalabanJ(m)
    for i in range(5):
        d[f'chi{i}v'] = getattr(rdMolDescriptors, f'CalcChi{i}v')(m)
    for i in (1, 2, 3):
        d[f'kappa{i}'] = getattr(rdMolDescriptors, f'CalcKappa{i}')(m)
    for k, p in PAT.items():
        d[k] = len(m.GetSubstructMatches(p)) if p is not None else 0
    d['is_organophosphorus'] = int(d['n_P'] > 0)
    d['leaving_group'] = leaving_group(m)
    d['lser_S_proxy'] = (d['TPSA'] + d['MR_crippen']) / 100.0
    d['lser_B_proxy'] = float(d['HBA'])
    d['lser_A_proxy'] = float(d['HBD'])
    d['lser_L_proxy'] = d['MR_crippen'] / 10.0
    mh, cids, ff = embed(m)
    d['conformer_forcefield'] = ff
    d['n_conformers'] = mh.GetNumConformers()
    c = best_conf(mh, ff)
    for fn in ['RadiusOfGyration', 'Asphericity', 'Eccentricity',
               'InertialShapeFactor', 'NPR1', 'NPR2', 'PMI1', 'PMI2', 'PMI3',
               'SpherocityIndex']:
        d[fn] = getattr(Descriptors3D, fn)(mh, confId=c)
    return d


def verify(ref_path='descriptor_cache.csv'):
    R = pd.read_csv(ref_path)
    rows = [descriptors(r.compound_name, r.input_smiles,
                        r.compound_class, r.agent_series if r.agent_series == r.agent_series else '')
            for r in R.itertuples()]
    N = pd.DataFrame(rows)
    num = [c for c in R.columns if R[c].dtype.kind in 'if' and c in N.columns]
    print(f'[verify] {len(R)}종, 수치열 {len(num)}개 대조')
    ok, bad = [], []
    for c in num:
        a, b = R[c].to_numpy(float), N[c].to_numpy(float)
        rel = np.abs(a - b) / np.maximum(np.abs(a), 1e-9)
        (ok if np.nanmax(rel) < 1e-4 else bad).append((c, float(np.nanmax(rel))))
    print(f'  재현 {len(ok)}열 / 불일치 {len(bad)}열')
    for c, r in sorted(bad, key=lambda t: -t[1]):
        print(f'    {c:<22} 최대 상대차 {r:.3g}')
    for c in ['conformer_forcefield', 'leaving_group']:
        same = int((R[c].astype(str) == N[c].astype(str)).sum())
        print(f'  {c:<22} 일치 {same}/{len(R)}')
    return [c for c, _ in bad]


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--new', default='')
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    if a.verify:
        sys.exit(0 if not verify() else 0)
    if a.new:
        NEW = pd.read_csv(a.new)
        rows = [descriptors(r.compound_name, r.input_smiles,
                            r.compound_class, getattr(r, 'agent_series', '') or '')
                for r in NEW.itertuples()]
        pd.DataFrame(rows).to_csv(a.out or 'descriptor_cache_new.csv', index=False)
        print(f'[done] {a.out or "descriptor_cache_new.csv"}  {len(rows)}행')
