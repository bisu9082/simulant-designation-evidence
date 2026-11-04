#!/usr/bin/env python3
"""
EPA Compiled AEGL (2018-07-27 판) 파싱 → 28종 독성 축.

원칙
- ppm 로 고시된 값만 읽고 mg/m3 로는 내가 환산한다 (MW/24.45, 25 C, 1 atm).
  PDF 의 대괄호 mg/m3 값은 줄바꿈으로 깨져 있어 파싱이 불안정하므로 읽지 않되,
  읽을 수 있는 것은 **환산 검증용**으로만 쓴다.
- mg/m3 로 고시된 항목은 그대로 쓴다.
- NR(Not recommended) 은 결측으로 둔다. 추정 금지.
"""
import json, re
import pandas as pd

import aegl_extract
import mw_std

PDF = 'sources/compiled_aegls_27jul2018.pdf'
MOLAR_VOL = 24.45                      # L/mol, 25 C, 1 atm

recs, _diag = aegl_extract.extract(PDF)
aegl_extract.report(_diag, recs)

# ------------------------------------------------- 28종 CAS 명부 (2026-10-11 신설)
# 종전 매칭은 이름 문자열 기준이었다. 원자료 1132행의 항목명은 'GF Agent' 인데
# 패턴은 'agent gf' 였으므로 매칭이 실패하고 cyclosarin 이 '미등재'로 기록됐다.
# 같은 원고의 QSAR 은 CAS 기준으로 같은 표를 읽어 0.38 mg/m3 을 쓰고 있었으므로
# 내부 모순이 발생했다. CWA 제외 규칙에서 이미 한 번 겪은 실패 모드다.
# 이제 CAS 를 1차 키로 쓰고, 이름 패턴은 교차검증으로만 쓴다.
PANEL_CAS = {
    'DMMP': '756-79-6',                 'acetone': '67-64-1',
    'methanol': '67-56-1',              'acetonitrile': '75-05-8',
    'chloroform': '67-66-3',            'n-hexane': '110-54-3',
    'THF': '109-99-9',                  'ethyl acetate': '141-78-6',
    'benzene': '71-43-2',               'water': '7732-18-5',
    'dichloromethane': '75-09-2',       'DPGME': '34590-94-8',
    'ethanol': '64-17-5',               'toluene': '108-88-3',
    'DMSO': '67-68-5',                  'p-dichlorobenzene': '106-46-7',
    'p-xylene': '106-42-3',             'sarin (GB)': '107-44-8',
    'soman (GD)': '96-64-0',            'tabun (GA)': '77-81-6',
    'sulfur mustard (HD)': '505-60-2',  'lewisite 1 (L)': '541-25-3',
    'hydrogen cyanide (AC)': '74-90-8', 'cyanogen chloride (CK)': '506-77-4',
    'phosgene (CG)': '75-44-5',         'chloropicrin (PS)': '76-06-2',
    'cyclosarin (GF)': '329-99-7',      'VX': '50782-69-9',
}

# ---------------------- 이름 패턴 (교차검증 전용)
MATCH = {
    'DMMP': ['dimethyl methylphosphonate'],
    'acetone': ['acetone'],
    'methanol': ['methanol'],
    'acetonitrile': ['acetonitrile'],
    'chloroform': ['chloroform'],
    'n-hexane': ['hexane'],
    'THF': ['tetrahydrofuran'],
    'ethyl acetate': ['ethyl acetate'],
    'benzene': ['benzene'],
    'water': [],
    'dichloromethane': ['methylene chloride', 'dichloromethane'],
    'DPGME': ['dipropylene glycol monomethyl ether', 'dipropylene glycol methyl ether'],
    'ethanol': ['ethanol', 'ethyl alcohol'],
    'toluene': ['toluene'],
    'DMSO': ['dimethyl sulfoxide'],
    'p-dichlorobenzene': ['1,4-dichlorobenzene', 'p-dichlorobenzene'],
    'p-xylene': ['p-xylene', 'xylenes'],
    'sarin (GB)': ['agent gb'],
    'soman (GD)': ['agent gd'],
    'tabun (GA)': ['agent ga'],
    'sulfur mustard (HD)': ['agent hd', 'sulfur mustard', 'mustard'],
    'lewisite 1 (L)': ['lewisite 1'],
    'hydrogen cyanide (AC)': ['hydrogen cyanide'],
    'cyanogen chloride (CK)': ['cyanogen chloride'],
    'phosgene (CG)': ['phosgene'],
    'chloropicrin (PS)': ['chloropicrin'],
    'cyclosarin (GF)': ['agent gf'],
    'VX': ['agent vx'],
}
low = [(r, r['name'].lower()) for r in recs]
df = pd.read_csv('descriptor_cache_xtb.csv')
MW = dict(zip(df.compound_name, df.MW))

# EPA 'Holding AEGLs' 절 등재 물질 (sources/compiled_aegls_27jul2018.pdf).
# 값표에는 없지만 심의 대상이고, 기관이 '현재 값을 개발하지 않는다'고 적은 것.
# 2026-10-23 신설. 종전에는 값표에 없으면 전부 '미등재'로 적었는데, 보충자료는
# 염화시안을 보류목록이라 쓰고 도판은 '미등재'를 찍어 1차 출처와 어긋났다.
# 본 패널 9종을 CAS 로 직접 대조한 결과 보류목록은 염화시안 하나뿐이다.
AEGL_HOLDING = {'506-77-4'}

out, unmatched, conv_check, xcheck = [], [], [], []
for _, row in df.iterrows():
    nm = row.compound_name
    # 1차: CAS 완전일치
    want_cas = PANEL_CAS[nm]
    by_cas = [r for r in recs if str(r.get('cas', '')).strip() == want_cas]
    hit = by_cas[0] if by_cas else None

    # 2차: 이름 패턴 (교차검증 전용 — CAS 결과와 다르면 기록)
    by_name = None
    for pat in MATCH[nm]:
        exact = [r for r, l in low if l == pat]
        if exact:
            by_name = exact[0]; break
        rx = re.compile(r'(?<![a-z0-9-])' + re.escape(pat) + r'(?![a-z0-9])')
        cands = [r for r, l in low if rx.search(l)]
        if cands:
            by_name = sorted(cands, key=lambda r: len(r['name']))[0]
            break
    if want_cas and (hit is not None) and (by_name is None):
        xcheck.append((nm, want_cas, hit['name'], 'CAS만 적중 — 이름 패턴 실패'))
    elif (hit is not None) and (by_name is not None) and \
            str(by_name.get('cas', '')).strip() != want_cas:
        xcheck.append((nm, want_cas, hit['name'],
                       f"이름 패턴은 다른 행 적중: {by_name['name']}"))
    elif (hit is None) and (by_name is not None):
        xcheck.append((nm, want_cas, by_name['name'], 'CAS 미적중인데 이름만 적중'))
        hit = by_name

    if hit is None:
        unmatched.append(nm)
        note = ('[수집 불가: EPA AEGL 보류목록, 값 미개발]'
                if want_cas in AEGL_HOLDING
                else '[수집 불가: EPA AEGL 미등재]')
        out.append(dict(compound_name=nm, compound_class=row.compound_class,
                        aegl_cas='', aegl_name='', aegl_unit='',
                        AEGL3_10min_mgm3=None, AEGL2_10min_mgm3=None,
                        AEGL3_provenance='', note=note))
        continue
    vals, prov = {}, {}
    for lvl in ('AEGL2', 'AEGL3'):
        e = hit['eff'].get(lvl)
        if e is None:
            vals[lvl] = None
            prov[lvl] = ''
            continue
        v10, unit, how = e
        prov[lvl] = how          # 'table' | 'footnote' — 도판·본문이 세는 근거
        if unit == 'ppm':
            vals[lvl] = mw_std.ppm_to_mgm3(v10, mw=MW[nm])
            br = hit['bracket'].get(lvl)
            if br and how == 'table':
                conv_check.append((nm, lvl, round(vals[lvl], 5), br[0]))
        else:
            vals[lvl] = v10
    out.append(dict(compound_name=nm, compound_class=row.compound_class,
                    aegl_cas=hit['cas'], aegl_name=hit['name'],
                    aegl_unit=hit['unit'],
                    AEGL3_10min_mgm3=vals['AEGL3'],
                    AEGL2_10min_mgm3=vals['AEGL2'],
                    AEGL3_provenance=prov.get('AEGL3', ''),
                    note='' if vals['AEGL3'] is not None
                         else '[수집 불가: AEGL-3 = NR]'))

T = pd.DataFrame(out)
T.to_csv('outputs/toxicity_aegl.csv', index=False)
print(f"\n[match] AEGL 등재 {int(T.aegl_name.astype(bool).sum())}/28, "
      f"AEGL-3 10분 확보 {int(T.AEGL3_10min_mgm3.notna().sum())}/28")
if unmatched:
    print('[match] 미등재:', ', '.join(unmatched))

print('\n[검증] 내 ppm->mg/m3 환산 vs EPA 자체 대괄호값')
for nm, lvl, mine, epa in conv_check:
    d = 100 * abs(mine - epa) / epa
    print(f'  {nm:<24} {lvl}  mine {mine:<10.5g} EPA [{epa}]  차이 {d:5.1f}%')

if xcheck:
    print('\n[교차검증] CAS 기준과 이름 패턴 기준의 불일치')
    for nm, cas, src, why in xcheck:
        print(f'  {nm:<22} CAS {cas:<12} 소스명 {src[:30]:<30} {why}')
else:
    print('\n[교차검증] CAS 기준 = 이름 패턴 기준, 불일치 0건')

print('\n[결과] AEGL-3 10분 (mg/m3) 오름차순 = 독성 강한 순')
S = T[T.AEGL3_10min_mgm3.notna()].sort_values('AEGL3_10min_mgm3')
print(S[['compound_name', 'compound_class', 'aegl_cas',
         'AEGL3_10min_mgm3']].to_string(index=False,
                                        float_format=lambda v: f'{v:.4g}'))
json.dump(dict(source='EPA Compiled AEGL values, 27 July 2018',
               source_url='https://www.epa.gov/sites/default/files/2018-08/'
                          'documents/compiled_aegls_update_27jul2018.pdf',
               metric='AEGL-3, 10-minute, mg/m3 (lower = more acutely toxic)',
               ppm_to_mgm3='value_ppm * MW / 24.45 (25 C, 1 atm)',
               n_listed=int(T.aegl_name.astype(bool).sum()),
               n_with_aegl3=int(T.AEGL3_10min_mgm3.notna().sum()),
               n_from_footnote=int((T.AEGL3_provenance == 'footnote').sum()),
               n_from_table=int((T.AEGL3_provenance == 'table').sum()),
               footnote_compounds=T[T.AEGL3_provenance == 'footnote']
                   .compound_name.tolist(),
               unmatched=unmatched,
               cas_name_crosscheck=[dict(compound=c, cas=k, source_name=s, note=w)
                                    for c, k, s, w in xcheck],
               conversion_check=[dict(compound=c, level=l, computed=m, epa_bracket=e)
                                 for c, l, m, e in conv_check]),
          open('outputs/toxicity_source.json', 'w'), indent=2, ensure_ascii=False)
print('\n[done] outputs/toxicity_aegl.csv, toxicity_source.json')
