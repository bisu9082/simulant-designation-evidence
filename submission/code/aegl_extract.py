#!/usr/bin/env python3
"""EPA Compiled AEGL PDF 의 레코드 추출 — 단일 경로 (2026-10-11 신설).

이 파일이 생긴 이유. 같은 PDF 를 세 곳에서 따로 파싱하고 있었다.
  parse_aegl.py          28종 패널용
  build_tox_trainset.py  247종 학습세트용
  validate_aegl_parse.py 검증용
세 경로가 각자 다른 결함을 가졌고, 서로의 결함을 덮어 주지 못했다. cyclosarin
이름 매칭 실패, 분자량 84.93/84.933 분기, 각주 행 오인이 모두 이 중복에서 나왔다.

이제 추출은 여기 하나뿐이다. parse_aegl.py 와 build_tox_trainset.py 가 이 함수를
쓴다. validate_aegl_parse.py 는 독립 검증이 목적이므로 일부러 공유하지 않고
다른 엔진(pdftotext -layout)과 다른 알고리즘(문자 열 위치)을 쓴다.
"""
import re

import pdfplumber

import aegl_table

PDF = 'sources/compiled_aegls_27jul2018.pdf'

SEC = re.compile(r'^(Final|Interim|Proposed|Holding) AEGLs\s*\((\d+)\)')
HDR = re.compile(r'^(?P<cas>[\d\-]{6,}(?:\s+and\s+[\d\-]{6,})?)?\s*\*?'
                 r'(?P<name>.+?)\s*\((?P<unit>ppm|mg/m3)\)')


def read_lines(pdf_path=PDF, first_page=7):
    with pdfplumber.open(pdf_path) as pdf:
        pages = [p.extract_text() or '' for p in pdf.pages[first_page:]]
    out = []
    for t in pages:
        for ln in t.split('\n'):
            ln = ln.strip()
            if ln and not ln.startswith('Page ') and not ln.startswith('July 27, 2018'):
                out.append(ln)
    return out


def extract(pdf_path=PDF, first_page=7):
    """레코드 리스트와 진단 정보를 돌려준다.

    각 레코드:
      cas, name, unit, section
      levels[AEGLn]  표 행의 5개 칸 (없으면 None)
      bracket[AEGLn] 기관이 병기한 mg/m3 등가값 (환산 검증용)
      foot[AEGLn]    각주에서 읽은 10분 값 (value, unit)
      merged         줄바꿈 병합이 일어난 레벨 목록
      eff[AEGLn]     (value, unit, 'table'|'footnote') 유효 10분 값
    """
    lines = read_lines(pdf_path, first_page)
    recs, cur, prev, section = [], None, None, 'unknown'
    fail_rows, foot_used = [], []

    def flush():
        if cur and (cur['levels'] or cur['foot']):
            recs.append(cur)

    for ln in lines:
        m = SEC.match(ln)
        if m:
            section = m.group(1)
            prev = None
            continue
        f = aegl_table.parse_footnote(ln)
        if f is not None and cur is not None:
            lvl, val, unit = f
            cur['foot'][f'AEGL{lvl}'] = (val, unit)
            prev = None
            continue
        r = aegl_table.parse_row(ln, prev)
        if r is not None and cur is not None:
            lvl, vals, status = r
            key = f'AEGL{lvl}'
            if vals is None:
                fail_rows.append((cur['cas'], cur['name'], key, status, ln[:70]))
                cur['levels'][key] = None
            else:
                cur['levels'][key] = vals
                if status == 'merged':
                    cur['merged'].append(key)
            br = aegl_table_bracket(ln)
            if br:
                cur['bracket'].setdefault(key, br)
            prev = None
            continue
        h = HDR.match(ln)
        if h and 'AEGL' not in ln and len(h.group('name')) < 95:
            flush()
            cur = dict(cas=(h.group('cas') or '').strip(),
                       name=h.group('name').strip(), unit=h.group('unit'),
                       section=section, levels={}, bracket={}, foot={},
                       merged=[], eff={})
            prev = None
            continue
        prev = ln if aegl_table.numeric_only(ln) else None
    flush()

    for r in recs:
        for key in ('AEGL1', 'AEGL2', 'AEGL3'):
            v = r['levels'].get(key)
            v10 = v[0] if (v and v[0] is not None) else None
            if v10 is not None:
                r['eff'][key] = (v10, r['unit'], 'table')
            elif key in r['foot']:
                val, u = r['foot'][key]
                r['eff'][key] = (val, u or r['unit'], 'footnote')
                if key == 'AEGL3':
                    foot_used.append((r['cas'], r['name'], val, u or r['unit']))

    return recs, dict(fail_rows=fail_rows, foot_used=foot_used,
                      merged=[(r['cas'], r['merged']) for r in recs if r['merged']])


def aegl_table_bracket(line):
    from aegl_num import bracket_nums
    return bracket_nums(line)


def report(diag, recs):
    print(f'[extract] 화학물질 항목 {len(recs)}건 '
          f"(ppm {sum(r['unit']=='ppm' for r in recs)} / "
          f"mg/m3 {sum(r['unit']=='mg/m3' for r in recs)})")
    n_tab = sum(1 for r in recs if r['eff'].get('AEGL3', (0, 0, ''))[2] == 'table')
    print(f"[extract] 10분 AEGL-3 유효값 "
          f"{sum(1 for r in recs if 'AEGL3' in r['eff'])}건 "
          f'(표 칸 {n_tab} / 각주 {len(diag["foot_used"])})')
    if diag['foot_used']:
        print('[extract] 10분 AEGL-3 을 각주에서 읽은 물질 '
              '— 표 칸은 see below 등으로 비어 있다')
        for c, n, v, u in diag['foot_used']:
            print(f'   {c:12s} {n[:40]:40s} {v:>10,.0f} {u}')
    if diag['merged']:
        print(f'[extract] 줄바꿈으로 쪼개진 행을 합친 항목 {len(diag["merged"])}건: '
              + ', '.join(f"{c}({','.join(k)})" for c, k in diag['merged']))
    if diag['fail_rows']:
        print(f'[extract] 셀 5개를 만들지 못해 값을 쓰지 않은 행 '
              f'{len(diag["fail_rows"])}건')
        for c, n, key, st, raw in diag['fail_rows']:
            print(f'   {c:12s} {n[:28]:28s} {key} {st:10s} | {raw}')


if __name__ == '__main__':
    recs, diag = extract()
    report(diag, recs)
