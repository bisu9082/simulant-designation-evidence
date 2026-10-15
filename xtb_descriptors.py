#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulant Grammar — GFN2-xTB 전자 서술자 산출 (Ubuntu / WSL).

입력 : descriptor_cache.csv   (compound_name, input_smiles 열 사용)
출력 : descriptor_cache_xtb.csv   — 아래 8열을 채운 사본

  dipole_moment_D            쌍극자 모멘트 (Debye)
  E_HOMO_eV / E_LUMO_eV      프런티어 궤도 에너지
  electrophilicity_index_eV  ω = μ² / (2η),  μ=(E_H+E_L)/2,  η=(E_L−E_H)/2
  q_P                        인 원자 Mulliken 전하
  q_phosphoryl_O             P=O 산소 Mulliken 전하
  wiberg_PO                  P=O Wiberg 결합차수
  fukui_f_plus_P             인 원자 친핵성공격 감수성 Fukui 지수 f(+)

형태 처리: ETKDGv3 5 conformer → xtb --opt tight → 전자 물성은
          xtb 전체에너지 기반 298 K Boltzmann 가중 평균.

사용법
  python3 xtb_descriptors.py                      # 전체 28종
  python3 xtb_descriptors.py --limit 3            # 연기 테스트
  python3 xtb_descriptors.py --nproc 8            # 병렬 스레드
"""
import argparse, csv, math, os, re, shutil, subprocess, sys, tempfile, time
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

HARTREE_EV = 27.211386245988
R_HARTREE  = 3.166811563e-6      # Hartree/K
T          = 298.15
N_CONF     = 5
SEED       = 42

PATT_PO = Chem.MolFromSmarts('[P]=[O]')
PATT_P  = Chem.MolFromSmarts('[P]')


# ───────────────────────── xtb 실행 ─────────────────────────
def run(cmd, cwd, timeout=1800):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                          timeout=timeout)


def write_xyz(mol, confId, path):
    conf = mol.GetConformer(confId)
    with open(path, 'w') as f:
        f.write('%d\n\n' % mol.GetNumAtoms())
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write('%-3s %18.10f %18.10f %18.10f\n'
                    % (a.GetSymbol(), p.x, p.y, p.z))


def parse_energy(txt):
    m = re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh', txt)
    return float(m[-1]) if m else None


def parse_homo_lumo(txt):
    """'(HOMO)' / '(LUMO)' 태그 바로 앞 숫자가 eV 값. 점유수 열 유무에 무관."""
    h = re.findall(r'(-?\d+\.\d+)\s*\(HOMO\)', txt)
    l = re.findall(r'(-?\d+\.\d+)\s*\(LUMO\)', txt)
    return (float(h[-1]) if h else None), (float(l[-1]) if l else None)


def parse_dipole(txt):
    m = re.search(r'molecular dipole:.*?full:\s+'
                  r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)',
                  txt, re.S)
    return float(m.group(4)) if m else None


def parse_fukui(txt):
    """xtb --vfukui 표. 6.7.x 는 번호와 원소기호가 '1O' 처럼 붙어서 출력된다."""
    out = {}
    lines = txt.splitlines()
    start = None
    for i, ln in enumerate(lines):
        if 'f(+)' in ln and 'f(-)' in ln:
            start = i + 1
            break
    if start is None:
        return out
    pat = re.compile(r'^\s*(\d+)\s*([A-Za-z]{1,2})\s+'
                     r'(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$')
    for ln in lines[start:]:
        m = pat.match(ln)
        if not m:
            if ln.strip() == '' and out:
                break
            continue
        out[int(m.group(1)) - 1] = float(m.group(3))   # 0-based → f(+)
    return out


def read_charges(d, n):
    p = os.path.join(d, 'charges')
    if not os.path.exists(p):
        return None
    v = [float(x) for x in open(p).read().split()]
    return v if len(v) == n else None


def read_wbo(d):
    p = os.path.join(d, 'wbo')
    if not os.path.exists(p):
        return {}
    w = {}
    for ln in open(p):
        t = ln.split()
        if len(t) >= 3:
            i, j = int(t[0]) - 1, int(t[1]) - 1
            w[(min(i, j), max(i, j))] = float(t[2])
    return w


# ───────────────────────── 분자 1종 ─────────────────────────
def one_molecule(name, smi, nproc, keep=None):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'SMILES 파싱 실패'
    mh = Chem.AddHs(mol)
    p = AllChem.ETKDGv3(); p.randomSeed = SEED; p.pruneRmsThresh = 0.3
    cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=N_CONF, params=p))
    if not cids:
        p.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mh, numConfs=N_CONF, params=p))
    if not cids:
        return None, 'conformer 생성 실패'
    nat = mh.GetNumAtoms()

    # P / phosphoryl O 인덱스 (수소 추가 후에도 heavy atom 순서는 보존됨)
    iP = mh.GetSubstructMatches(PATT_P)
    iP = iP[0][0] if iP else None
    po = mh.GetSubstructMatches(PATT_PO)
    iPO = po[0][1] if po else None

    recs = []
    base = keep or tempfile.mkdtemp(prefix='xtb_')
    for k, cid in enumerate(cids):
        d = os.path.join(base, '%s_c%d' % (re.sub(r'\W', '_', name)[:24], k))
        os.makedirs(d, exist_ok=True)
        xyz = os.path.join(d, 'in.xyz')
        write_xyz(mh, cid, xyz)
        env = dict(os.environ, OMP_NUM_THREADS=str(nproc),
                   MKL_NUM_THREADS=str(nproc), OMP_STACKSIZE='4G')
        r1 = subprocess.run(['xtb', 'in.xyz', '--gfn', '2', '--opt', 'tight',
                             '--chrg', '0', '--uhf', '0'],
                            cwd=d, capture_output=True, text=True, env=env,
                            timeout=3600)
        opt = os.path.join(d, 'xtbopt.xyz')
        geom = 'xtbopt.xyz' if os.path.exists(opt) else 'in.xyz'
        # 중성 단일점 — 궤도·쌍극자·전하·WBO 는 전부 여기서만 읽는다.
        # (--vfukui 런은 N+1 / N-1 전자 SCF 를 함께 돌려 궤도표가 오염된다)
        rsp = subprocess.run(['xtb', geom, '--gfn', '2', '--chrg', '0', '--uhf', '0'],
                             cwd=d, capture_output=True, text=True, env=env,
                             timeout=3600)
        E = parse_energy(rsp.stdout) or parse_energy(r1.stdout)
        if E is None:
            continue
        hE, lE = parse_homo_lumo(rsp.stdout)
        dp = parse_dipole(rsp.stdout)
        qs = read_charges(d, nat)
        wb = read_wbo(d)
        # Fukui 전용 런
        rfk = subprocess.run(['xtb', geom, '--gfn', '2', '--chrg', '0',
                              '--uhf', '0', '--vfukui'],
                             cwd=d, capture_output=True, text=True, env=env,
                             timeout=3600)
        recs.append(dict(E=E, homo=hE, lumo=lE, dip=dp, q=qs, wbo=wb,
                         fk=parse_fukui(rfk.stdout)))
    if not recs:
        return None, 'xtb 전체 실패'

    E = np.array([r['E'] for r in recs]); E -= E.min()
    w = np.exp(-E / (R_HARTREE * T)); w /= w.sum()

    def avg(key):
        v = np.array([r[key] if r[key] is not None else np.nan for r in recs],
                     dtype=float)
        ok = np.isfinite(v)
        return float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else ''

    out = {'dipole_moment_D': avg('dip'),
           'E_HOMO_eV': avg('homo'), 'E_LUMO_eV': avg('lumo')}
    if out['E_HOMO_eV'] != '' and out['E_LUMO_eV'] != '':
        mu  = (out['E_HOMO_eV'] + out['E_LUMO_eV']) / 2.0
        # Parr 정의: eta = 간극의 절반 (2026-10-12 교정)
        eta = (out['E_LUMO_eV'] - out['E_HOMO_eV']) / 2.0
        out['electrophilicity_index_eV'] = mu * mu / (2 * eta) if eta > 0 else ''
    else:
        out['electrophilicity_index_eV'] = ''

    def avg_at(getter):
        v, ww = [], []
        for r, wi in zip(recs, w):
            x = getter(r)
            if x is not None and np.isfinite(x):
                v.append(x); ww.append(wi)
        return float(np.dot(v, ww) / sum(ww)) if ww else ''

    out['q_P'] = avg_at(lambda r: r['q'][iP] if (r['q'] and iP is not None) else None) if iP is not None else ''
    out['q_phosphoryl_O'] = avg_at(lambda r: r['q'][iPO] if (r['q'] and iPO is not None) else None) if iPO is not None else ''
    if iP is not None and iPO is not None:
        key = (min(iP, iPO), max(iP, iPO))
        out['wiberg_PO'] = avg_at(lambda r: r['wbo'].get(key))
    else:
        out['wiberg_PO'] = ''
    out['fukui_f_plus_P'] = avg_at(lambda r: r['fk'].get(iP)) if iP is not None else ''
    out['_n_conf_xtb'] = len(recs)
    if keep is None:
        shutil.rmtree(base, ignore_errors=True)
    return out, None


# ───────────────────────── main ─────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input',  default='descriptor_cache.csv')
    ap.add_argument('--output', default='descriptor_cache_xtb.csv')
    ap.add_argument('--nproc',  type=int, default=4)
    ap.add_argument('--limit',  type=int, default=0)
    ap.add_argument('--keep',   default='')
    a = ap.parse_args()

    if shutil.which('xtb') is None:
        sys.exit('xtb 를 PATH 에서 찾을 수 없습니다.  conda install -c conda-forge xtb')
    v = subprocess.run(['xtb', '--version'], capture_output=True, text=True)
    print('xtb:', (v.stdout + v.stderr).strip().splitlines()[-1] if (v.stdout or v.stderr) else '?')

    rows = list(csv.DictReader(open(a.input, encoding='utf-8-sig')))
    if a.limit: rows = rows[:a.limit]
    COLS = ['dipole_moment_D','E_HOMO_eV','E_LUMO_eV','electrophilicity_index_eV',
            'q_P','q_phosphoryl_O','wiberg_PO','fukui_f_plus_P']
    fields = list(rows[0].keys())
    for c in COLS + ['_n_conf_xtb']:
        if c not in fields: fields.append(c)

    t0 = time.time(); fails = []
    for i, r in enumerate(rows, 1):
        nm, smi = r['compound_name'], r['input_smiles']
        t1 = time.time()
        try:
            res, err = one_molecule(nm, smi, a.nproc, a.keep or None)
        except subprocess.TimeoutExpired:
            res, err = None, '시간 초과'
        except Exception as e:
            res, err = None, repr(e)[:70]
        if res is None:
            fails.append((nm, err))
            print('[%2d/%2d] %-26s  실패: %s' % (i, len(rows), nm[:26], err), flush=True)
            continue
        r.update({k: ('' if res[k] == '' else round(res[k], 5)) for k in COLS})
        r['_n_conf_xtb'] = res['_n_conf_xtb']
        print('[%2d/%2d] %-26s conf=%d  HOMO %7s  LUMO %7s  w %6s  dip %6s  q(P) %7s  f+(P) %7s  %5.1fs'
              % (i, len(rows), nm[:26], res['_n_conf_xtb'],
                 r['E_HOMO_eV'], r['E_LUMO_eV'], r['electrophilicity_index_eV'],
                 r['dipole_moment_D'], r['q_P'], r['fukui_f_plus_P'],
                 time.time() - t1), flush=True)

    with open(a.output, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

    print('\n%s — %d행 / 총 %.1f분' % (a.output, len(rows), (time.time()-t0)/60))
    got = sum(1 for r in rows if str(r.get('E_LUMO_eV','')) != '')
    print('전자 서술자 확보 %d / %d' % (got, len(rows)))
    nP = sum(1 for r in rows if str(r.get('q_P','')) != '')
    print('인 함유(q_P 산출) %d종' % nP)
    if fails:
        print('실패:', fails)

if __name__ == '__main__':
    main()
