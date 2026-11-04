#!/usr/bin/env python3
"""
빠른 전자 서술자 (스크리닝용): ETKDG 1-conformer -> MMFF 최적화 -> GFN2 단일점.
HOMO / LUMO / omega / dipole 만 뽑는다. xtb 기하 최적화·Fukui 생략.

Ku 데스크탑의 정밀 모드(--opt tight + 5컨포머 Boltzmann)와의 일치도는
동일 28종을 양쪽으로 돌려 상관으로 검증한다 (compare 모드).
"""
import argparse, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

XTB = os.environ.get('XTB_BIN', 'xtb')
XTB_SHARE = os.environ.get('XTB_SHARE', '')
OPT = os.environ.get('XTB_OPT', '1') == '1'
HARTREE_EV = 27.211386245988


def parse_orbitals(txt):
    h = re.findall(r'(-?\d+\.\d+)\s*\(HOMO\)', txt)
    l = re.findall(r'(-?\d+\.\d+)\s*\(LUMO\)', txt)
    return (float(h[-1]) if h else None), (float(l[-1]) if l else None)


def parse_dipole(txt):
    m = re.search(r'molecular dipole:.*?\n.*?\n.*?full:\s*'
                  r'[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+([\d.]+)', txt, re.S)
    return float(m.group(1)) if m else None


def one(args):
    name, smiles = args
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return dict(name=name, ok=False, err='SMILES 파싱 실패')
        m = Chem.AddHs(m)
        p = AllChem.ETKDGv3(); p.randomSeed = 42
        if AllChem.EmbedMolecule(m, p) != 0:
            p.useRandomCoords = True
            if AllChem.EmbedMolecule(m, p) != 0:
                return dict(name=name, ok=False, err='3D 임베딩 실패')
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        except Exception:
            AllChem.UFFOptimizeMolecule(m, maxIters=2000)
        chg = Chem.GetFormalCharge(m)
        nrad = sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())
        d = tempfile.mkdtemp(prefix='xtbf_')
        try:
            xyz = os.path.join(d, 'in.xyz')
            Chem.MolToXYZFile(m, xyz)
            env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
            if XTB_SHARE:
                env['XTBPATH'] = XTB_SHARE
            base = [XTB, '--gfn', '2', '--chrg', str(chg), '--uhf', str(nrad)]
            if OPT:
                subprocess.run(base[:1] + ['in.xyz'] + base[1:] + ['--opt', 'normal'],
                               cwd=d, env=env, capture_output=True, text=True,
                               timeout=600)
                geom = 'xtbopt.xyz' if os.path.exists(os.path.join(d, 'xtbopt.xyz')) \
                    else 'in.xyz'
            else:
                geom = 'in.xyz'
            r = subprocess.run(base[:1] + [geom] + base[1:], cwd=d, env=env,
                               capture_output=True, text=True, timeout=600)
            txt = r.stdout + r.stderr
            ho, lu = parse_orbitals(txt)
            dip = parse_dipole(txt)
            if ho is None or lu is None:
                return dict(name=name, ok=False, err='궤도 파싱 실패')
            # Parr 정의: eta 는 간극의 절반이다 (2026-10-12 교정).
            # 종전에는 간극 전체를 넣어 omega 가 전부 절반이었다.
            mu, eta = (ho + lu) / 2.0, (lu - ho) / 2.0
            return dict(name=name, ok=True, E_HOMO_eV=ho, E_LUMO_eV=lu,
                        gap_eV=lu - ho, dipole_moment_D=dip,
                        electrophilicity_index_eV=(mu * mu) / (2 * eta) if eta > 0 else None,
                        n_heavy=m.GetNumHeavyAtoms() - sum(a.GetSymbol() == 'H'
                                                           for a in m.GetAtoms()) * 0)
        finally:
            shutil.rmtree(d, ignore_errors=True)
    except subprocess.TimeoutExpired:
        return dict(name=name, ok=False, err='xtb 시간초과')
    except Exception as e:
        return dict(name=name, ok=False, err=f'{type(e).__name__}: {e}')


def run(pairs, nproc):
    out, n = [], len(pairs)
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for i, r in enumerate(ex.map(one, pairs), 1):
            out.append(r)
            if i % 10 == 0 or i == n:
                nf = sum(not x['ok'] for x in out)
                print(f'  ... {i}/{n}  실패 {nf}', flush=True)
    return pd.DataFrame(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['aegl', 'compare'], required=True)
    ap.add_argument('--nproc', type=int, default=8)
    A = ap.parse_args()

    if A.mode == 'compare':
        C = pd.read_csv('descriptor_cache_xtb.csv')
        pairs = list(zip(C.compound_name, C.canonical_smiles))
        print(f'[compare] 28종 빠른 모드 재계산 (정밀 모드와 대조)')
        R = run(pairs, A.nproc)
        M = R.merge(C[['compound_name', 'E_HOMO_eV', 'E_LUMO_eV',
                       'electrophilicity_index_eV', 'dipole_moment_D']]
                    .rename(columns=lambda c: c if c == 'compound_name' else c + '_slow'),
                    left_on='name', right_on='compound_name')
        M.to_csv('outputs/xtb_fast_vs_slow.csv', index=False)
        print(f'\n[compare] 성공 {int(M.ok.sum())}/28')
        for c in ['E_HOMO_eV', 'E_LUMO_eV', 'electrophilicity_index_eV',
                  'dipole_moment_D']:
            s = M[M.ok][[c, c + '_slow']].dropna()
            r = np.corrcoef(s[c], s[c + '_slow'])[0, 1]
            mad = float(np.abs(s[c] - s[c + '_slow']).mean())
            print(f'  {c:<28} Pearson r = {r:.4f}   평균절대차 {mad:.4f}  (n={len(s)})')
    else:
        D = pd.read_csv('outputs/aegl_trainset_raw.csv')
        pairs = list(zip(D.cas, D.smiles))
        print(f'[aegl] {len(pairs)}종 빠른 모드')
        R = run(pairs, A.nproc).rename(columns={'name': 'cas'})
        M = D.merge(R, on='cas', how='left')
        M.to_csv('outputs/aegl_electronic.csv', index=False)
        print(f'\n[aegl] 성공 {int(M.ok.fillna(False).sum())}/{len(M)}')
        bad = M[~M.ok.fillna(False)]
        if len(bad):
            print('  실패:', ', '.join(f'{r.cas}({r.err})' for r in bad.itertuples()))
    print('[done]')
