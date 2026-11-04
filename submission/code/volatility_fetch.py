#!/usr/bin/env python3
"""증기압 후보값 수집 (2026-11-03 신설). 수집 전용. 판정하지 않는다.

왜 분리하는가. 리뷰어 1 이 '서론이 휘발성을 제일 먼저 말하는데 서술자
패널에 휘발성이 없다'고 지적했다. 축을 넣으려면 측정된 증기압이 필요하고,
[C-3] 에 따라 추정값은 쓸 수 없다. 그런데 공개 DB 는 같은 화합물에 온도와
출처가 다른 값을 여러 개 담고 있다 (사린 2.86 mmHg @25C Epstein 1974 /
2.9 mmHg @25C EPA 1998, 루이사이트는 이성질체별로 0.4 와 1.562 로 4배).
따라서 '수집'과 '채택'을 분리한다.

- 이 스크립트: PubChem 의 Vapor Pressure 항목 전문과 출처 이름을 그대로
  긁어 outputs/volatility_candidates.json 에 담는다. 아무 값도 고르지
  않는다.
- volatility_space.py: 사람이 읽고 확정한 채택값을 출처 인용과 함께
  리터럴로 들고 있다. 그 파일이 정본이다.

출처 우선순위는 volatility_space.py 에 적는다.
"""
import json
import os
import time
import urllib.parse
import urllib.request

PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
VIEW = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound'

# SMILES 조회가 실패하는 화합물은 CAS 로 받는다. 루이사이트는 PubChem 의
# SMILES 정규화가 trans/cis 를 받아들이지 않아 400 을 돌려준다.
CAS_OVERRIDE = {'lewisite 1 (L)': '541-25-3'}


def _get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'simulant/1.0'})
    with urllib.request.urlopen(req, timeout=45) as fh:
        return fh.read().decode('utf-8', 'replace')


def cid_for(name, smiles):
    if name in CAS_OVERRIDE:
        q = urllib.parse.quote(CAS_OVERRIDE[name], safe='')
        return _get(f'{PUG}/compound/name/{q}/cids/TXT').strip().split('\n')[0]
    q = urllib.parse.quote(smiles, safe='')
    return _get(f'{PUG}/compound/smiles/{q}/cids/TXT').strip().split('\n')[0]


def vapor_pressure(cid):
    """(출처, 값 문자열) 목록. 항목이 없으면 빈 목록."""
    try:
        d = json.loads(_get(f'{VIEW}/{cid}/JSON?heading=Vapor+Pressure'))
    except Exception:
        return []
    refs = {r.get('ReferenceNumber'): r.get('SourceName', '')
            for r in d.get('Record', {}).get('Reference', [])}
    out = []

    def walk(o):
        if isinstance(o, dict):
            if o.get('TOCHeading') == 'Vapor Pressure':
                for inf in o.get('Information', []):
                    cite = inf.get('Reference')
                    if cite:
                        src = str(cite[0])
                    else:
                        src = refs.get(inf.get('ReferenceNumber'), '?')
                    for x in inf.get('Value', {}).get('StringWithMarkup', []):
                        out.append({'source': src, 'value': x.get('String', '')})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(d)
    return out


def main():
    import csv
    with open('panel_all.csv') as fh:
        rows = list(csv.DictReader(fh))
    res = []
    for r in rows:
        name = r['compound_name']
        try:
            cid = cid_for(name, r['input_smiles'])
        except Exception as e:
            cid = f'ERR {e}'
        time.sleep(0.4)
        vals = [] if cid.startswith('ERR') else vapor_pressure(cid)
        time.sleep(0.4)
        res.append({'name': name, 'klass': r['compound_class'], 'cid': cid,
                    'candidates': vals})
        print(f'{name:<24} cid={cid:<12} 후보 {len(vals)}')
    os.makedirs('outputs', exist_ok=True)
    with open('outputs/volatility_candidates.json', 'w') as fh:
        json.dump({'generated_by': 'volatility_fetch.py',
                   'source': 'PubChem PUG-View, heading=Vapor Pressure',
                   'note': '수집 전용. 채택값은 volatility_space.py 가 정본.',
                   'compounds': res}, fh, indent=1, ensure_ascii=False)
    miss = [r['name'] for r in res if not r['candidates']]
    print(f'\n후보 있음 {len(res) - len(miss)}/{len(res)}')
    print('후보 없음:', miss)


if __name__ == '__main__':
    main()
