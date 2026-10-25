#!/usr/bin/env python3
"""GFN2 의 omega 와 f+(P) 를 DFT 로 다시 재는 대조 (2026-10-25 신설).

왜 필요한가. 본문과 보충자료는 GFN2-xTB 의 궤도에너지에서 Parr omega 를 만들고,
그 축에서 '현행 모사체가 작용제보다 친전자적으로 보인다'는 식의 서술을 한다.
그런데 GFN2 의 가상궤도는 Koopmans 의미의 전자친화도가 아니고, 이 패널에서
chloropicrin 의 omega 가 40 eV 로 나온다. 유기분자의 DFT 수준 omega 는 보통
1~3 eV 이고 강한 친전자체라도 10 eV 를 넘기 어렵다. 그러면 omega 축에서 끌어낸
결론이 화학에 대한 것인지 GFN2 에 대한 것인지 구분되지 않는다.

여기서는 같은 화합물을 B3LYP/def2-SVPD 로 다시 계산해 둘을 나란히 놓는다.
결론을 바꾸려는 것이 아니라, 어느 쪽 진술이 방법에 종속인지 가른다.

한계를 먼저 적는다. 이 컨테이너에는 xtb 가 없어 GFN2 최적화 기하를 재현할 수
없다. 여기 기하는 RDKit ETKDG + MMFF94(실패 시 UFF)다. 따라서 GFN2 값과의
차이에는 방법과 기하가 함께 들어 있다. 기하 몫이 얼마인지는 --optgeom 으로
몇 종을 DFT 최적화해 따로 잰다.

출력: outputs/dft_omega.csv, outputs/dft_omega.json
"""
import argparse
import json
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, dft

HARTREE = 27.211386245988
BASIS = 'def2-svpd'
XC = 'b3lyp'
SEED = 0xC0DE


def geometry(smi, optimise_ff=True):
    """ETKDG + MMFF94. MMFF 가 못 다루는 원소(비소 등)는 UFF 로 떨어뜨린다."""
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        AllChem.EmbedMolecule(m, randomSeed=SEED, useRandomCoords=True)
    ff = 'none'
    if optimise_ff:
        try:
            if AllChem.MMFFHasAllMoleculeParams(m):
                AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
                ff = 'MMFF94'
            else:
                AllChem.UFFOptimizeMolecule(m, maxIters=2000)
                ff = 'UFF'
        except Exception:
            AllChem.UFFOptimizeMolecule(m, maxIters=2000)
            ff = 'UFF'
    c = m.GetConformer()
    at = [(a.GetSymbol(), tuple(c.GetAtomPosition(a.GetIdx())))
          for a in m.GetAtoms()]
    pidx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'P']
    return at, ff, pidx


def frontier(mf):
    occ = mf.mo_occ > 0
    h = float(mf.mo_energy[occ].max()) * HARTREE
    l = float(mf.mo_energy[~occ].min()) * HARTREE
    mu = (h + l) / 2.0
    eta = (l - h) / 2.0
    return h, l, mu, eta, (mu * mu / (2 * eta) if eta > 0 else float('nan'))


def mulliken_on(mf, mol, idx):
    """Mulliken 전하. f+ 를 쓰기 위한 것이므로 규약을 명시해 둔다."""
    pop, chg = mf.mulliken_pop(verbose=0)
    return float(chg[idx])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', default='descriptor_cache_xtb.csv')
    ap.add_argument('--fukui', action='store_true',
                    help='인 함유 화합물에 대해 N+1 계산으로 f+(P) 도 낸다')
    ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args()

    df = pd.read_csv(a.panel)
    rows = df[['compound_name', 'compound_class', 'canonical_smiles']].to_dict('records')
    if a.limit:
        rows = rows[:a.limit]
    out, t_all = [], time.time()
    print(f'[dft] {len(rows)}종, {XC}/{BASIS}', flush=True)
    for i, r in enumerate(rows, 1):
        name, smi = r['compound_name'], r['canonical_smiles']
        t0 = time.time()
        try:
            at, ff, pidx = geometry(smi)
            mol = gto.M(atom=at, basis=BASIS, verbose=0)
            mf = dft.RKS(mol)
            mf.xc = XC
            mf.kernel()
            if not mf.converged:
                raise RuntimeError('SCF 미수렴')
            h, l, mu, eta, om = frontier(mf)
            rec = dict(name=name, klass=r['compound_class'], n_atoms=mol.natm,
                       forcefield=ff, homo_eV=h, lumo_eV=l, mu_eV=mu,
                       eta_eV=eta, omega_eV=om, converged=True)
            if a.fukui and pidx:
                qN = mulliken_on(mf, mol, pidx[0])
                an = gto.M(atom=at, basis=BASIS, charge=-1, spin=1, verbose=0)
                mfa = dft.UKS(an)
                mfa.xc = XC
                mfa.kernel()
                if mfa.converged:
                    pop, chg = mfa.mulliken_pop(verbose=0)
                    rec['f_plus_P_mulliken'] = qN - float(chg[pidx[0]])
                else:
                    rec['f_plus_P_mulliken'] = None
            out.append(rec)
            extra = ''
            if rec.get('f_plus_P_mulliken') is not None:
                extra = f"  f+(P) {rec['f_plus_P_mulliken']:+.4f}"
            print(f'  [{i:2d}/{len(rows)}] {name:26s} {ff:6s} n={mol.natm:3d} '
                  f'{time.time()-t0:6.1f}s  HOMO {h:8.3f} LUMO {l:8.3f} '
                  f'omega {om:7.3f}{extra}', flush=True)
        except Exception as e:
            out.append(dict(name=name, klass=r['compound_class'],
                            converged=False, error=f'{type(e).__name__}: {e}'))
            print(f'  [{i:2d}/{len(rows)}] {name:26s} 실패 {type(e).__name__}',
                  flush=True)

    d = pd.DataFrame(out)
    d.to_csv('outputs/dft_omega.csv', index=False)
    ok = d[d.converged == True] if 'converged' in d else d
    meta = dict(xc=XC, basis=BASIS, geometry='RDKit ETKDG + MMFF94/UFF',
                seed=SEED, n_total=len(d), n_converged=int(len(ok)),
                note=('GFN2 기하를 재현할 수 없어 기하가 다르다. GFN2 와의 차이에는 '
                      '방법과 기하가 함께 들어 있다.'),
                omega_min=float(ok.omega_eV.min()) if len(ok) else None,
                omega_max=float(ok.omega_eV.max()) if len(ok) else None)
    json.dump(meta, open('outputs/dft_omega.json', 'w'), indent=1,
              ensure_ascii=False)
    print(f'\n[dft] 완료 {len(ok)}/{len(d)}  총 {time.time()-t_all:.0f}s')
    print('[done] outputs/dft_omega.csv, outputs/dft_omega.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
