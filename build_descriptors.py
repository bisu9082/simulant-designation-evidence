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
    # 2026-11-07 (10차 R5-B5). 종전에는 없는 캐시를 continue 로 넘겼다.
    # 캐시가 없는 패키지에서 관문이 조용히 건너뛰고도 '전 수치열 재현' 을
    # 찍었고, 풀어 돌린 사람에게는 통과로 보였다. 없는 것은 중단 사유이고,
    # 재계산(28종 x 수십 초)보다 먼저 확인한다. 반경험 캐시는 이 관문이
    # 재계산하지 못하지만, 없으면 거리 계산이 전자공간을 조용히 빼고
    # 돌아가므로 존재만이라도 확인한다.
    missing = [_r for _r in refs if not _os.path.exists(_r)]
    missing += [_x for _x in ('descriptor_cache_xtb.csv',
                              'descriptor_cache_xtb_35.csv',
                              'descriptor_cache_xtb_36.csv',
                              'descriptor_cache_xtb_ascl3.csv')
                if not _os.path.exists(_x)]
    if missing:
        print(f'[중단] 서술자 캐시 {len(missing)}개가 없다: '
              + ', '.join(missing)
              + '\n  관문이 건너뛰면 통과로 보이지만 아무것도 대조하지 않은 '
                '것이다. 예치 패키지에서 누락된 파일을 받을 것.')
        return 1
    bad = []
    for _r in refs:
        print(f'\n[gate] {_r}')
        bad += [f'{_r}:{c}' for c in verify(_r)]
    if bad:
        print(f'\n[중단] 재현되지 않는 열 {len(bad)}개. 거리 계산에 들어가는 '
              '열이 코드에서 재생성되지 않으면 분석을 신뢰할 수 없다.')
        return 1
    # 2026-11-07. 관문의 범위를 산출물로 남긴다. SI 가 '53 of the 61' 이라고
    # 적으므로 그 두 수가 손으로 적힌 값이 아니라 이 실행의 결과여야 한다.
    # verify_consistency.py 가 이 파일을 읽어 본문과 대조한다.
    import json as _j
    _R = pd.read_csv(refs[0])
    _META = ('compound_name', 'compound_class', 'agent_series', 'input_smiles',
             'canonical_smiles', 'formula', 'leaving_group',
             'conformer_forcefield')
    _all = [c for c in _R.columns
            if _R[c].dtype.kind in 'if' and c not in _META]
    _os.makedirs('outputs', exist_ok=True)
    _j.dump({'ref': refs[0], 'n_numeric': len(_all),
             'n_checked': len([c for c in _all if c not in XTB_EXEMPT]),
             'exempt': list(XTB_EXEMPT)},
            open('outputs/descriptor_gate_scope.json', 'w'), indent=1)
    print('\n[gate] RDKit 서술자 전 수치열 재현. 통과.')
    print('[gate] 범위: 반경험 전자 좌표 8열은 이 관문이 재계산하지 않는다 '
          '(xtb 바이너리 필요, 버전 미기록). SI meth:dft 의 B3LYP 대조가 '
          '그 축을 따로 검사한다.')
    return 0


# 2026-11-07 (10차 R5-B5). 서술자 관문이 재계산하지 못하는 열. 전부
# GFN2-xTB 또는 그 파생이고, descriptors() 는 RDKit 경로만 돈다. 여기에
# 열을 더하려면 왜 재계산할 수 없는지를 함께 적을 것.
XTB_EXEMPT = ('E_HOMO_eV', 'E_LUMO_eV', 'dipole_moment_D',
              'electrophilicity_index_eV', 'fukui_f_plus_P',
              'q_P', 'q_phosphoryl_O', 'wiberg_PO')


def verify(ref_path='descriptor_cache.csv'):
    R = pd.read_csv(ref_path)
    rows = [descriptors(r.compound_name, r.input_smiles,
                        r.compound_class, r.agent_series if r.agent_series == r.agent_series else '')
            for r in R.itertuples()]
    N = pd.DataFrame(rows)
    # 2026-11-02. 종전에는 R 쪽 dtype 만 보았다. descriptor_cache_ascl3.csv 는
    # 1행이고 agent_series 가 비어 있어 pandas 가 그 문자열 열을 float NaN 으로
    # 읽는다. 그러면 num 에 들어가고, 새로 계산한 값은 '' 이므로 to_numpy(float)
    # 가 ValueError 로 죽는다. 빈 outputs/ 로 클론해 돌리면 관문이 그 자리에서
    # 멈춘다. 양쪽 모두에서 수치인 열만 대조한다.
    num = [c for c in R.columns
           if R[c].dtype.kind in 'if' and c in N.columns
           and pd.to_numeric(N[c], errors='coerce').notna().any()]
    # 2026-11-07 (10차 R5-B5). 종전에는 대조한 열 수만 찍고 '전 수치열 재현'
    # 이라고 통과를 선언했다. 그런데 캐시의 수치열은 61개이고 descriptors()
    # 가 내놓는 것은 53개다. 빠진 8개는 전부 반경험·전자 좌표이고, 그 셋은
    # 본문 electronic space 의 전부다. 즉 관문은 '전 수치열' 을 본 적이
    # 없는데 그렇게 적고 있었다. 빠진 열을 이름으로 찍고, 선언한 면제
    # 목록에 없는 열이 빠지면 중단한다.
    # 식별·메타 열은 서술자가 아니다. descriptor_cache_ascl3.csv 는 1행이고
    # agent_series 가 비어 있어 pandas 가 float NaN 열로 읽는다(위 주석과 같은
    # 자리). 그것까지 '미대조 서술자' 로 세면 거짓 중단이 난다.
    META = ('compound_name', 'compound_class', 'agent_series', 'input_smiles',
            'canonical_smiles', 'formula', 'leaving_group',
            'conformer_forcefield')
    allnum = [c for c in R.columns
              if R[c].dtype.kind in 'if' and c not in META]
    uncov = sorted(set(allnum) - set(num))
    print(f'[verify] {len(R)}종, 수치열 {len(allnum)}개 중 {len(num)}개 대조')
    if uncov:
        print(f'  미대조 {len(uncov)}열 (재계산에 xtb 바이너리가 필요): '
              + ', '.join(uncov))
        extra = [c for c in uncov if c not in XTB_EXEMPT]
        if extra:
            print(f'  [중단] 면제 목록에 없는 미대조 열: {", ".join(extra)}')
            return extra
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
