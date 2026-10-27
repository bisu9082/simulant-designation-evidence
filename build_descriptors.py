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
import os as _os
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


D3_NAMES = ['RadiusOfGyration', 'Asphericity', 'Eccentricity',
            'InertialShapeFactor', 'NPR1', 'NPR2', 'PMI1', 'PMI2', 'PMI3',
            'SpherocityIndex']

# 298.15 K 에서의 RT. MMFF/UFF 에너지가 kcal/mol 이므로 같은 단위를 쓴다.
RT_KCAL = 0.5922


def conf_energies(mh, ff):
    if ff == 'MMFF94' and AllChem.MMFFHasAllMoleculeParams(mh):
        props = AllChem.MMFFGetMoleculeProperties(mh)
        return np.array([
            AllChem.MMFFGetMoleculeForceField(mh, props, confId=c).CalcEnergy()
            for c in range(mh.GetNumConformers())], float)
    return np.array([
        AllChem.UFFGetMoleculeForceField(mh, confId=c).CalcEnergy()
        for c in range(mh.GetNumConformers())], float)


def ensemble_3d(mh, ff):
    """3D 서술자를 형태이성질체 앙상블의 Boltzmann 가중 평균으로 낸다.

    2026-10-25 변경. 종전에는 best_conf() 로 MMFF 최저에너지 형태이성질체
    하나를 골라 거기서 읽었다. 그 값은 ETKDG 시드에 의존하며, --verify 가
    52개 수치열 중 11개를 재현 불가로 판정하는 원인이었다. 그중 9개가 본문
    구조공간에 들어가고, 28종 패널의 네 번째 유효매칭(tabun)이 그 자유도
    위에서 역치를 넘나들었다 (시드 6개 중 1개에서만 성립).

    앙상블 평균은 단일 선택을 없애므로 시드 의존을 크게 줄인다. 완전히
    없애지는 않는다. 앙상블 자체가 시드의 함수이기 때문이다. 남은 폭은
    conformer_sensitivity.py 가 잰다.
    """
    e = conf_energies(mh, ff)
    w = np.exp(-(e - e.min()) / RT_KCAL)
    w = w / w.sum()
    out = {}
    for fn in D3_NAMES:
        v = np.array([getattr(Descriptors3D, fn)(mh, confId=c)
                      for c in range(mh.GetNumConformers())], float)
        out[fn] = float((w * v).sum())
    # 앙상블이 한 형태이성질체에 몰리면 평균이 단일 선택과 같아진다.
    # 보고용으로 최대 가중치를 남긴다.
    out['_w_max'] = float(w.max())
    return out


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
    ens = ensemble_3d(mh, ff)
    d['conformer_w_max'] = ens.pop('_w_max')
    d.update(ens)
    return d


def gate_verify():
    """재현 검사를 실제 관문으로 쓴다 (2026-10-25).

    종전에는 sys.exit(0 if not verify() else 0) 이라 어느 쪽이든 0 을
    반환했다. 그래서 52열 중 11열이 재현되지 않는 상태로 파이프라인이
    계속 돌았고, 그 열 9개가 본문의 모든 거리에 들어가 있었다.
    """
    import sys as _s
    refs = ['descriptor_cache.csv', 'descriptor_cache_35.csv',
            'descriptor_cache_ascl3.csv']
    for _a in _s.argv[1:]:
        if _a.startswith('--ref='):
            refs = [_a.split('=', 1)[1]]
    bad = []
    for _r in refs:
        if not _os.path.exists(_r):
            continue
        print(f'\n[gate] {_r}')
        bad += [f'{_r}:{c}' for c in verify(_r)]
    if bad:
        print(f'\n[중단] 재현되지 않는 열 {len(bad)}개. 거리 계산에 들어가는 '
              '열이 코드에서 재생성되지 않으면 분석을 신뢰할 수 없다.')
        return 1
    print('\n[gate] 전 수치열 재현. 통과.')
    return 0


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
        sys.exit(gate_verify())
    if a.new:
        NEW = pd.read_csv(a.new)
        rows = [descriptors(r.compound_name, r.input_smiles,
                            r.compound_class, getattr(r, 'agent_series', '') or '')
                for r in NEW.itertuples()]
        pd.DataFrame(rows).to_csv(a.out or 'descriptor_cache_new.csv', index=False)
        print(f'[done] {a.out or "descriptor_cache_new.csv"}  {len(rows)}행')
