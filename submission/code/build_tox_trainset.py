#!/usr/bin/env python3
"""
AEGL-3 (10분) 예측 모델용 학습셋 구축.

1) EPA Compiled AEGL PDF 전 항목에서 CAS + 단위 + AEGL-3 10분 수치 추출
2) PubChem PUG REST 로 CAS -> SMILES, MW 해석 (초당 ~4건)
3) ppm 항목은 MW 로 mg/m3 환산 (25 C, 1 atm)

이중용도 가드: 화학작용제(CWA)는 학습에서 제외하고 held-out 평가점으로만 쓴다.
여기서는 섹션·CAS 만 기록하고, 제외는 모델 적합 단계에서 수행한다.
"""
import json, re, time
import pandas as pd

import aegl_extract
import mw_std
import requests

PDF = 'sources/compiled_aegls_27jul2018.pdf'
MOLAR_VOL = 24.45

recs, _diag = aegl_extract.extract(PDF)
aegl_extract.report(_diag, recs)

rows = []
for r in recs:
    e3 = r['eff'].get('AEGL3')
    e2 = r['eff'].get('AEGL2')
    if e3 is None:
        continue
    rows.append(dict(cas=r['cas'], aegl_name=r['name'], unit=e3[1],
                     section=r['section'],
                     aegl3_10min_raw=e3[0], aegl3_provenance=e3[2],
                     aegl2_10min_raw=(e2[0] if e2 else None)))
D = pd.DataFrame(rows)
_nocas = D[D.cas.astype(str).str.strip() == '']
if len(_nocas):
    print(f'[extract] CAS 가 없는 항목 {len(_nocas)}건 제외 '
          '(원표에 등록번호가 없어 구조를 일의적으로 해석할 수 없다): '
          + ', '.join(_nocas.aegl_name.astype(str).str[:40]))
# 빈 CAS 끼리는 서로 다른 물질이므로 중복제거 대상이 아니다. 먼저 제외한 뒤
# 남은 항목만 CAS 로 중복제거한다.
D = D[D.cas.astype(str).str.strip() != ''].drop_duplicates('cas').reset_index(drop=True)
print(f'[extract] AEGL-3 10분 수치 보유 {len(D)}건')
print('  섹션별:', D.section.value_counts().to_dict())
print('  단위별:', D.unit.value_counts().to_dict())
print('  출처별:', D.aegl3_provenance.value_counts().to_dict())

# ------------------------------------------------------- PubChem CAS -> SMILES
S = requests.Session()
BASE = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}'
        '/property/SMILES,MolecularWeight/JSON')

# 이전 실행 결과를 CAS 캐시로 재사용한다. 파서 교정으로 신규 CAS 가 늘었으므로
# 전건 재조회는 불필요하고, 기존 해석을 바꿀 이유도 없다.
import os
CACHE = 'outputs/pubchem_cas_cache.csv'
cache = {}
if os.path.exists(CACHE):
    c = pd.read_csv(CACHE)
    cache = {str(r.cas).strip(): (r.smiles, r.pubchem_mw) for r in c.itertuples()}
elif os.path.exists('outputs/aegl_trainset_raw.csv'):
    c = pd.read_csv('outputs/aegl_trainset_raw.csv')
    cache = {str(r.cas).strip(): (r.smiles, r.pubchem_mw) for r in c.itertuples()}
print(f'[pubchem] 캐시 {len(cache)}건 적재, 신규 조회 대상 '
      f'{sum(1 for c_ in D.cas if str(c_).strip() not in cache)}건')

smi, mw, fails = [], [], []
for i, cas in enumerate(D.cas):
    key = str(cas).strip()
    if key in cache:
        s_, m_ = cache[key]
        smi.append(s_ if isinstance(s_, str) else None)
        mw.append(m_ if m_ == m_ else None)
        continue
    got = (None, None)
    for attempt in range(3):
        try:
            r = S.get(BASE.format(key), timeout=25)
            if r.status_code == 200:
                pr = r.json()['PropertyTable']['Properties'][0]
                got = (pr.get('SMILES'),
                       float(pr['MolecularWeight']) if pr.get('MolecularWeight') else None)
            break
        except Exception:
            time.sleep(1.5)
    if got[0] is None:
        fails.append(key)
    smi.append(got[0]); mw.append(got[1])
    time.sleep(0.25)
    if (i + 1) % 40 == 0:
        print(f'  ... {i+1}/{len(D)} 해석, 실패 누적 {len(fails)}')
print(f'[pubchem] 미해석 {len(fails)}건: {fails[:12]}')
D['smiles'], D['pubchem_mw'] = smi, mw
print(f'[pubchem] 성공 {D.smiles.notna().sum()}/{len(D)}  실패 {len(fails)}')

# 분자량은 mw_std 가 단일 기준이다 (RDKit MolWt, 발표된 SMILES 에서 재계산 가능).
# PubChem 값은 pubchem_mw 로 보존해 교차검증에만 쓴다.
D['mw_rdkit'] = D.smiles.map(mw_std.mol_weight)
# RDKit 이 구조를 파싱하지 못하는 항목(할로겐간화합물 등)은 PubChem 분자량으로
# 대체해 원표 수집 단계에서는 보존한다. 서술자를 만들 수 없으므로 모델에는
# 들어가지 못하지만, 그 제외는 사유와 함께 qsar_tox.py 에서 일어나야 하고
# 수집 단계에서 조용히 사라져서는 안 된다.
D['mw_source'] = ['rdkit' if m == m and m is not None else 'pubchem_fallback'
                  for m in D.mw_rdkit]
D['mw'] = [r if (r == r and r is not None) else p_
           for r, p_ in zip(D.mw_rdkit, D.pubchem_mw)]
D['aegl3_10min_mgm3'] = [
    mw_std.convert(v, u, mw=m)
    for v, u, m in zip(D.aegl3_10min_raw, D.unit, D.mw)]
_fb = int((D.mw_source == 'pubchem_fallback').sum())
print(f'[mw] RDKit {int((D.mw_source=="rdkit").sum())} / '
      f'PubChem 대체 {_fb}  (대체분은 서술자 생성 불가로 모델에서 제외된다)')
D = D.dropna(subset=['smiles', 'aegl3_10min_mgm3']).reset_index(drop=True)
D.to_csv('outputs/aegl_trainset_raw.csv', index=False)
print(f'[done] 구조+엔드포인트 확보 {len(D)}건 -> outputs/aegl_trainset_raw.csv')
print(D[['cas', 'aegl_name', 'section', 'unit', 'aegl3_10min_mgm3']]
      .head(12).to_string(index=False, float_format=lambda v: f'{v:.4g}'))
json.dump(dict(n_extracted=int(len(recs)), n_with_aegl3=int(len(D)),
               pubchem_failures=fails), open('outputs/trainset_build.json', 'w'),
          indent=2, ensure_ascii=False)
