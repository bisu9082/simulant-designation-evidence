#!/usr/bin/env python3
"""ppm -> mg/m3 환산의 단일 기준 (2026-10-11 신설).

배경. 같은 EPA 값(메틸렌 클로라이드 AEGL-3 10분 = 12,000 ppm)을 두 파이프라인이
각자 다른 분자량으로 환산해 본문에 41,683 과 41,685 가 동시에 나타났다.
  - 28종 경로: descriptor_cache_xtb.csv 의 RDKit MolWt  84.933
  - 247종 경로: PubChem property MolecularWeight         84.93
차이는 반올림이며 전 항목 최대 0.03% 이지만, 같은 수치가 본문에서 두 값으로
나타나는 것 자체가 결함이다. 분자량은 원고가 발표하는 SMILES 에서 재계산
가능해야 하므로, 양쪽 모두 RDKit MolWt 를 쓴다. PubChem 값은 교차검증용으로만
보존한다.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog('rdApp.*')

# 25 C, 1 atm 에서의 기체 몰부피 (L/mol)
MOLAR_VOL = 24.45


def mol_weight(smiles):
    """표준 분자량. 발표된 SMILES 에서 재계산 가능한 유일한 값."""
    m = Chem.MolFromSmiles(str(smiles))
    return float(Descriptors.MolWt(m)) if m is not None else None


def ppm_to_mgm3(value_ppm, smiles=None, mw=None):
    """ppm 농도를 mg/m3 로 환산한다. mw 를 직접 주면 그 값을 쓴다."""
    if value_ppm is None or value_ppm != value_ppm:
        return None
    if mw is None:
        mw = mol_weight(smiles)
    if mw is None or mw != mw:
        return None
    return float(value_ppm) * float(mw) / MOLAR_VOL


def convert(value, unit, smiles=None, mw=None):
    """단위를 보고 환산한다. 이미 mg/m3 면 그대로 돌려준다."""
    if value is None or value != value:
        return None
    u = str(unit).strip().lower()
    if u == 'mg/m3':
        return float(value)
    if u == 'ppm':
        return ppm_to_mgm3(value, smiles=smiles, mw=mw)
    return None
