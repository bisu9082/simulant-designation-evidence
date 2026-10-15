#!/usr/bin/env python3
"""화합물 동정 검증 — 이름이 함의하는 원소가 저장 구조에 실제로 있는가 (2026-10-11 신설).

배경. AEGL 학습세트의 SMILES 는 PubChem 을 화합물 '이름'으로 조회해 얻었고,
"Pentaborane" 이 EPN(유기인계 살충제)으로 해석되어 AEGL-3 라벨이 5.1배 틀렸다.
ppm 토큰 왕복 검증은 환산 '입력'(동정·분자량)을 보지 않으므로 이를 탐지하지 못한다.

CAS 역색인(xref/RN) 재조회도 단독으로는 신뢰할 수 없음을 확인했다.
247건 중 37건이 왕복 실패했는데, 그중 36건은 저장 구조가 오히려 옳았고
역색인 쪽이 엉뚱한 CID 를 돌려준 경우였다 (styrene -> fructose 등).

따라서 여기서는 외부 조회에 의존하지 않는 결정적 검사를 쓴다.
화합물 이름의 어간이 함의하는 원소가 저장 구조의 분자식에 존재하는지 본다.
이름이 borane 이면 B 가 있어야 하고, silane 이면 Si 가, phosph 면 P 가 있어야 한다.
전혀 다른 화합물로 해석된 경우는 이 검사를 통과하지 못한다.
"""
import json
import re
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog('rdApp.*')

# 이름 어간 -> 반드시 존재해야 하는 원소
IMPLIES = [
    (r'bor(ane|on)|borate', 'B'),
    (r'silane|silic|siloxane', 'Si'),
    (r'phosph(ine|ide|ate|onate|oro)', 'P'),
    (r'chlor', 'Cl'),
    (r'fluor', 'F'),
    (r'brom', 'Br'),
    (r'iod', 'I'),
    (r'sulf|thio|mercaptan|sulph', 'S'),
    (r'nitr(o|ic|ous|ogen)|amine|amide|cyan|azine|hydrazine|isocyanate', 'N'),
    (r'arsen|arsine', 'As'),
    (r'mercur', 'Hg'),
    (r'nickel', 'Ni'),
    (r'zinc', 'Zn'),
    (r'alumin', 'Al'),
    (r'calcium', 'Ca'),
    (r'sodium', 'Na'),
    (r'potassium', 'K'),
    (r'selen', 'Se'),
    (r'tellur', 'Te'),
]
# 어간이 있어도 해당 원소가 없을 수 있는 예외 (관용명 등)
EXEMPT_NAME = re.compile(r'HFC|HCFC|CFC|halon|freon', re.I)


# 배수 접두사 + 할로겐 어간 -> 해당 원소 개수 일치 확인
MULT = {'mono': 1, 'di': 2, 'tri': 3, 'tetra': 4, 'penta': 5, 'hexa': 6}
HAL = {'chloro': 'Cl', 'fluoro': 'F', 'bromo': 'Br', 'iodo': 'I'}


def primary_name(name):
    """괄호 안 동의어를 떼어낸다.

    'Ethylene chlorohydrin (2-Chloroethanol)' 는 같은 화합물을 두 번 부르므로
    할로겐 어간이 두 번 집계되어 Cl 2개를 요구했다. 구조는 Cl 1개로 옳았는데
    검증기가 FAIL 로 잘못 판정했다. 괄호 앞의 주명칭만 센다.
    """
    return re.sub(r'\s*\([^)]*\)', '', str(name)).strip() or str(name)


def count_check(name, formula):
    """이름에 언급된 할로겐 개수를 모두 합산해 구조와 대조한다.

    이름 한 개에 할로겐 어간이 여러 번 나올 수 있다
    ("Chloromethyl trichlorosilane" = 1 + 3 = Cl 4개).
    구조가 이름보다 '적으면' 다른 화합물일 가능성이 높으므로 FAIL,
    '많으면' 치환체 누락 표기일 수 있으므로 ADVISORY 로 분리한다.
    """
    if formula is None:
        return None, None
    want = {}
    for m in re.finditer(
            r'(mono|di|tri|tetra|penta|hexa)?\s*(chloro|fluoro|bromo|iodo|'
            r'chloride|fluoride|bromide|iodide)', primary_name(name), re.I):
        stem = m.group(2).lower()
        el = HAL.get(stem) or HAL.get(stem.replace('ide', 'o'))
        if el is None:
            continue
        want[el] = want.get(el, 0) + MULT.get((m.group(1) or '').lower(), 1)
    low, high = [], []
    for el, w in want.items():
        g = re.search(el + r'(\d*)', formula)
        have = (int(g.group(1)) if g.group(1) else 1) if g else 0
        if have < w:
            low.append(f'{el}:{have}<{w}')
        elif have > w:
            high.append(f'{el}:{have}>{w}')
    return (','.join(low) or None), (','.join(high) or None)


def elements(smiles):
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None, None
    f = rdMolDescriptors.CalcMolFormula(m)
    els = set(re.findall(r'[A-Z][a-z]?', re.sub(r'[+\-\d]', '', f)))
    return f, els


def check(name, smiles):
    f, els = elements(smiles)
    if els is None:
        return dict(status='UNPARSEABLE', formula=None, missing=None)
    if EXEMPT_NAME.search(str(name)):
        return dict(status='EXEMPT', formula=f, missing=None)
    miss = [el for pat, el in IMPLIES
            if re.search(pat, primary_name(name), re.I) and el not in els]
    low, high = count_check(name, f)
    # 다성분 레코드 여부. 염·수화물처럼 이름이 부가 성분을 함의하는 경우가 아니라면,
    # 부가 성분은 ppm->mg/m3 환산에 쓰는 분자량을 그만큼 틀리게 만든다.
    multi = '.' in str(smiles)
    if miss or low:
        why = ','.join([x for x in [','.join(miss) if miss else None, low] if x])
        return dict(status='FAIL', formula=f, missing=why)
    if high and multi:
        # 예: 27137-85-5 Trichlorophenylsilane 의 저장 구조는 PhSiCl3 + 2 HCl 이고
        # 분자량 284.47 로 실제 211.55 보다 34% 크다. 같은 화합물이 98-13-5 에
        # 정상 구조로 따로 등재돼 있다. 단성분 레코드의 명명 축약(Chloroform 등)과
        # 달리 환산 입력이 틀리므로 FAIL 이다.
        return dict(status='FAIL', formula=f,
                    missing=f'{high} (다성분 레코드 — 환산 분자량 오류)')
    if high:
        return dict(status='ADVISORY', formula=f, missing=high)
    return dict(status='PASS', formula=f, missing=None)


if __name__ == '__main__':
    T = pd.read_csv('outputs/aegl_trainset_raw.csv')
    out = []
    for r in T.itertuples():
        c = check(r.aegl_name, r.smiles)
        c.update(cas=r.cas, aegl_name=r.aegl_name, mw=r.pubchem_mw)
        out.append(c)
    R = pd.DataFrame(out)
    R.to_csv('outputs/identity_validation.csv', index=False)
    n = len(R)
    fail = R[R.status == 'FAIL']
    print(f'[identity] {n}건 검사  PASS {int((R.status=="PASS").sum())}  '
          f'EXEMPT {int((R.status=="EXEMPT").sum())}  '
          f'UNPARSEABLE {int((R.status=="UNPARSEABLE").sum())}  '
          f'ADVISORY {int((R.status=="ADVISORY").sum())}  FAIL {len(fail)}')
    for r in R[R.status == 'ADVISORY'].itertuples():
        print(f'   ADVISORY {r.cas:12s} {str(r.aegl_name)[:34]:34s} '
              f'구조 {r.formula:18s} 이름보다 많음 {r.missing}')
    for r in fail.itertuples():
        print(f'   FAIL {r.cas:12s} {str(r.aegl_name)[:34]:34s} '
              f'구조 {r.formula:18s} 누락원소 {r.missing}')
    json.dump(dict(n=n, n_fail=int(len(fail)),
                   failures=[dict(cas=r.cas, name=r.aegl_name, formula=r.formula,
                                  missing=r.missing) for r in fail.itertuples()]),
              open('outputs/identity_validation.json', 'w'),
              ensure_ascii=False, indent=1)
    sys.exit(1 if len(fail) else 0)
