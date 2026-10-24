#!/usr/bin/env python3
"""AEGL 10분 값 전표 독립 검증 (2026-10-11 전면 재작성).

종전 검증기는 'AEGL 3' 줄의 숫자를 순서대로 읽어 첫 숫자를 10분 값으로 삼았다.
그래서 10분 칸이 '**see below' 로 채워진 물질에서 30분 값을 10분 값으로 당겨
읽었다 (아세톤 8,600 / 벤젠 5,600 / 메탄올 14,000). 주 파서가 이들을 수집하지
못했으므로 대조 대상에서 빠져 '불일치 0건' 판정이 나왔지만, 검증기 자체가
쉼표 결함과 같은 계열의 칼럼 밀림을 갖고 있었다.

따라서 이번에는 추출 엔진과 알고리즘을 모두 주 파서와 다르게 한다.
  - 엔진: 주 파서는 pdfplumber, 검증기는 pdftotext -layout
  - 알고리즘: 주 파서는 토큰열 문법, 검증기는 문자 열 위치
    '   10 min    30 min   60 min    4 hr    8 hr' 헤더에서 각 지속시간 라벨의
    문자 구간을 구하고, 'AEGL 3' 줄의 같은 구간을 잘라 읽는다. 칼럼이 밀리면
    잘라낸 구간에 다른 값이 들어오므로 반드시 드러난다.
  - 각주: 10분 칸이 비어 있는 물질은 각주를 읽되, 주 파서와 달리 등호 뒤
    전체를 다시 토큰화하지 않고 '10 min' 라벨이 등호 왼쪽에 있는지 먼저 확인한다.
"""
import json
import re
import subprocess
import sys

import pandas as pd

import mw_std

PDF = 'sources/compiled_aegls_27jul2018.pdf'
HDR_DUR = re.compile(r'10\s*min.*30\s*min.*60\s*min.*4\s*hr.*8\s*hr', re.I)
CAS_HDR = re.compile(r'^\s*(?P<cas>\d{2,7}-\d{2}-\d(?:\s+and\s+\d{2,7}-\d{2}-\d)?)\s+'
                     r'(?P<name>.+?)\s*\((?P<unit>ppm|mg/m3)\)')
A3 = re.compile(r'^\s*AEGL\s*3\s+(?![-–])', re.I)
FOOT3 = re.compile(r'^\s*AEGL\s*3\s*[-–]\s*(?P<spec>[^=]*)=\s*(?P<rest>.*)$', re.I)
PROSE3 = re.compile(r'10[-\s]*min\s+AEGL[-\s]*3\s+value\s+of\s+'
                    r'(?P<v>[\d,]+(?:\.\d+)?)\s*(?P<u>ppm|mg/m3)', re.I)
NUMTOK = re.compile(r'^\**\s*([\d,]+(?:\.\d+)?)\s*\**$')


def layout_text():
    r = subprocess.run(['pdftotext', '-layout', PDF, '-'],
                       capture_output=True, text=True, check=True)
    return r.stdout.split('\n')


def ten_min_span(header_line):
    """헤더에서 10분 칸의 끝 문자 위치를 구한다.

    시작은 행마다 'AEGL n' 라벨을 공백으로 지운 뒤 0 으로 둔다. 라벨 들여쓰기가
    물질마다 다르기 때문이다 (' AEGL 3' / '  AEGL 3'). 종전에는 10분 라벨
    위치에서 8칸을 당겨 시작점을 잡았는데, 들여쓰기가 한 칸 더 깊은 물질에서
    구간이 라벨과 겹쳐 'AEGL 3' 의 3 을 값으로 읽었다
    (Vinyl trichlorosilane: 실제 210 ppm 을 3 ppm 으로).
    30분 라벨 시작 직전까지가 10분 칸이다.
    """
    m10 = re.search(r'10\s*min', header_line, re.I)
    m30 = re.search(r'30\s*min', header_line, re.I)
    if not (m10 and m30):
        return None
    return 0, m30.start() - 1


def read_cell(line, span):
    s = line[span[0]:span[1]]
    if re.search(r'see\s+below', s, re.I):
        return None, 'see_below'
    if re.search(r'\b(NR|ND)\b', s, re.I):
        return None, 'nr'
    toks = [t for t in s.split() if t.strip()]
    toks = [t for t in toks if not re.fullmatch(r'\*+', t)]
    if not toks:
        return None, 'empty'
    m = NUMTOK.match(toks[0])
    if not m:
        return None, f'nonnumeric:{toks[0][:14]}'
    return float(m.group(1).replace(',', '')), 'value'


def collect():
    """CAS -> dict(name, unit, value, how)"""
    L = layout_text()
    out, cur, span = {}, None, None
    for i, ln in enumerate(L):
        h = CAS_HDR.match(ln)
        # 목차 줄은 점선 리더로 끝난다. 2026-10-24 까지 이 줄들도 레코드를
        # 열었고, 목차 마지막 항목(Ricin)은 뒤이은 본문 표에 CAS 헤더가 없어
        # 다음 화합물(Magnesium aluminum phosphide)의 AEGL-3 2.4 ppm 을
        # 받아 갔다. 본문의 진짜 Ricin 헤더는 중복 회피 규칙에 걸려 영구히
        # 건너뛰어졌다. 레코드 경계 혼입이며, 공저자 심사에서 적발됐다.
        if h and '....' in ln:
            cur, span = None, None      # 열린 레코드도 닫는다. 닫지 않으면
            continue                    # 목차 항목이 본문 표까지 살아남는다.
        if h and 'AEGL' not in ln:
            cas = h.group('cas').strip()
            was = out.get(cas)
            if was is not None and was['value'] is not None:
                # 문서 말미의 변경이력 부록이 같은 CAS 를 'Replaced ... / Added ...'
                # 설명과 함께 재등장시킨다. 그 줄로 레코드를 새로 만들면 앞서
                # 확보한 값이 지워진다 (Vinyl acetate, Ethylenimine, Methanol).
                cur, span = None, None
                continue
            cur = dict(cas=cas, name=h.group('name').strip(),
                       unit=h.group('unit'), value=None, how=None)
            out[cas] = cur
            span = None
            continue
        if cur is None:
            continue
        if HDR_DUR.search(ln) and span is None:
            span = ten_min_span(ln)
            continue
        if A3.match(ln) and span and cur['how'] is None:
            # 'AEGL 3' 라벨을 같은 길이의 공백으로 치환해 문자 위치를 보존한다.
            lab = re.match(r'^\s*AEGL\s*3', ln, re.I)
            masked = ' ' * lab.end() + ln[lab.end():]
            v, how = read_cell(masked, span)
            cur['value'], cur['how'] = v, how
            continue
        f = FOOT3.match(ln)
        if f and re.search(r'(?<!\d)10\s*min', f.group('spec'), re.I):
            rest = f.group('rest').split(';')[0]
            g = re.search(r'([\d,]+(?:\.\d+)?)', rest)
            if g:
                u = re.search(r'(ppm|mg/m3)', rest, re.I)
                cur['value'] = float(g.group(1).replace(',', ''))
                cur['how'] = 'footnote'
                if u:
                    cur['unit'] = u.group(1).lower()
            continue
        pm = PROSE3.search(ln)
        if pm and cur['value'] is None:
            cur['value'] = float(pm.group('v').replace(',', ''))
            cur['how'] = 'footnote_prose'
            cur['unit'] = pm.group('u').lower()
    return out


if __name__ == '__main__':
    src = collect()
    have = {c: r for c, r in src.items() if r['value'] is not None}
    print(f'[검증] pdftotext -layout 문자열위치 방식으로 읽은 항목 {len(src)}건, '
          f'10분 AEGL-3 보유 {len(have)}건')
    byhow = {}
    for r in have.values():
        byhow[r['how']] = byhow.get(r['how'], 0) + 1
    print(f'[검증] 출처별: {byhow}')

    T = pd.read_csv('outputs/aegl_trainset_raw.csv')
    T['cas'] = T.cas.astype(str).str.strip()
    checked = mism = notfound = 0
    bad = []
    for r in T.itertuples():
        s = have.get(r.cas)
        if s is None:
            notfound += 1
            continue
        checked += 1
        exp = mw_std.convert(s['value'], s['unit'], mw=r.mw)
        if exp is None or abs(exp - r.aegl3_10min_mgm3) > max(
                0.001, 0.002 * abs(r.aegl3_10min_mgm3)):
            mism += 1
            bad.append((r.cas, str(r.aegl_name)[:30], s['value'], s['unit'],
                        s['how'], exp, r.aegl3_10min_mgm3))
    print(f'[검증] 학습세트 {checked}건 대조, 원전 미발견 {notfound}건, '
          f'불일치 {mism}건')
    for c, n, v, u, how, exp, got in bad:
        print(f'   불일치 {c:12s} {n:30s} 원전 {v:>10,.4g} {u} ({how}) '
              f'-> {exp:,.2f}  저장 {got:,.2f}')

    # 28종 패널 대조
    P = pd.read_csv('outputs/toxicity_aegl.csv')
    P = P[P.AEGL3_10min_mgm3.notna()]
    pchecked = pmism = 0
    pbad = []
    MWP = dict(zip(pd.read_csv('descriptor_cache_xtb.csv').compound_name,
                   pd.read_csv('descriptor_cache_xtb.csv').MW))
    print(f'\n[검증] 28종 중 공시값 보유 {len(P)}건 — 원전 10분 칸 직접 대조')
    for r in P.itertuples():
        s = have.get(str(r.aegl_cas).strip())
        if s is None:
            pbad.append((r.compound_name, '原典미발견', None, None))
            continue
        pchecked += 1
        exp = mw_std.convert(s['value'], s['unit'], mw=MWP[r.compound_name])
        ok = exp is not None and abs(exp - r.AEGL3_10min_mgm3) <= max(
            0.001, 0.002 * abs(r.AEGL3_10min_mgm3))
        flag = '' if ok else '  <<< 불일치'
        if not ok:
            pmism += 1
        print(f'   {r.compound_name:<24} 원전 {s["value"]:>10,.4g} {s["unit"]:6s} '
              f'({s["how"]:14s}) -> {exp:>12,.2f}  저장 {r.AEGL3_10min_mgm3:>12,.2f}{flag}')

    # 2026-10-24. 종전에는 how=='value'(표 칸 경로)로 한정했다. 천단위 쉼표가
    # 실제로 산문 안에 박히는 것은 각주 경로이므로, 탐지하려는 실패 모드를
    # 표본에서 제외하고 있었다. 경로와 무관하게 전수로 넓힌다.
    comma = {c for c, r in have.items() if r['value'] >= 1000}
    comma_cell = {c for c in comma if have[c]['how'] == 'value'}
    cov = sum(1 for r in T.itertuples() if r.cas in comma)
    print(f'\n[검증] 천단위 쉼표 구간(>=1000) 원전 {len(comma)}건 '
          f'(표 칸 {len(comma_cell)}, 각주 {len(comma) - len(comma_cell)}) 중 '
          f'학습세트 대조 {cov}건')

    # 2026-10-24 신설. 종전에는 학습세트 257행을 단방향으로만 훑었고,
    # 학습세트에 없는 레코드(구조 미해결 등)는 대조 대상에 아예 없었다.
    # 그 결과 검증 파서 자신의 오류 한 건이 구조적으로 보이지 않았다.
    import io as _io, contextlib as _cl
    import aegl_extract as _E
    with _cl.redirect_stdout(_io.StringIO()):
        _recs, _ = _E.extract()
    prim = {}
    for _r in _recs:
        _c = str(_r.get('cas') or '').strip()
        _L = _r.get('levels', {}).get('AEGL3')
        if _c and _L and _L[0] is not None:
            prim[_c] = float(_L[0])
    pair = set(prim) & set(have)
    xmis = [(c, prim[c], have[c]['value']) for c in pair
            if abs(prim[c] - have[c]['value']) > max(1e-6, 0.002 * abs(prim[c]))]
    print(f'\n[검증] 두 리더 양방향 대조 — 양쪽이 표 칸에서 값을 얻은 {len(pair)}건, '
          f'불일치 {len(xmis)}건 (1차만 {len(set(prim) - set(have))}, '
          f'검증기만 {len(set(have) - set(prim))})')
    for c, x, y in sorted(xmis):
        print(f'   불일치 {c:14s} 1차 {x:>12,.4g}   검증기 {y:>12,.4g}')

    ok = (mism == 0 and pmism == 0 and len(xmis) == 0)
    json.dump(dict(engine='pdftotext -layout', method='character column position',
                   source_items=len(src), source_with_aegl3=len(have),
                   by_provenance=byhow,
                   trainset_checked=checked, trainset_notfound=notfound,
                   trainset_mismatch=mism,
                   panel_checked=pchecked, panel_mismatch=pmism,
                   comma_rows_source=len(comma),
                   comma_rows_cell_path=len(comma_cell),
                   comma_rows_footnote_path=len(comma) - len(comma_cell),
                   comma_rows_checked=cov,
                   cross_reader_pairs=len(pair), cross_reader_mismatch=len(xmis),
                   cross_reader_primary_only=len(set(prim) - set(have)),
                   cross_reader_validator_only=len(set(have) - set(prim))),
              open('outputs/aegl_parse_validation.json', 'w'),
              indent=2, ensure_ascii=False)
    print(f'\n[검증] 판정: {"PASS" if ok else "FAIL"}')
    sys.exit(0 if ok else 1)
