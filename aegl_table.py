#!/usr/bin/env python3
"""EPA Compiled AEGL 표의 셀 단위 파서 (2026-10-11 전면 재작성).

종전 파서는 'AEGL N' 으로 시작하는 줄에서 숫자를 모두 긁어 앞 5개를 취했다.
그 설계 때문에 결함 5종이 생겼다.

  D1 각주 행이 표 행으로 인식됐다. 각주도 'AEGL 3' 으로 시작하기 때문이다.
     'AEGL 3 - 10 min/30 min/60 min/4 hr/8 hr = *** 33,000 ppm' 에서
     숫자 [10, 30, 60, 4, 8, 33000] 의 앞 5개를 취해, 프로판의 10분 AEGL-3 이
     33,000 ppm 이 아니라 10 ppm 으로 저장됐다 (3,300배 오류).
     부탄은 [10, 77000, 30, 60, 4] 로 10 ppm 이 됐다 (7,700배).
  D2 10분 값이 표 칸 대신 각주에 적힌 물질 9종의 값을 전혀 수집하지 못했다.
     메탄올 40,000 ppm, 아세톤 16,000 ppm, 벤젠 9,700 ppm, n-헥산 12,000 ppm,
     자일렌 7,200 ppm, MEK 10,000 ppm, 프로판 33,000 ppm, 부탄 77,000 ppm,
     MTBE 13,000 ppm.
  D3 긴 지속시간이 NR 인 행을 전량 폐기했다. 'AEGL 3  0.52  0.17  0.086  NR  NR'
     은 숫자가 3개뿐이라 None 이 됐고, 멀쩡한 10분 값을 버렸다. 이 때문에
     Ethyldichloroarsine, Phenyl dichloroarsine, BZ, Ricin 등 군용작용제가
     누락됐다.
  D4 pdfplumber 가 행을 두 줄로 쪼개는 경우를 처리하지 못했다 (Dichlorvos 등).
  D5 숫자 개수에만 의존했으므로 칼럼이 밀려도 알 수 없었다.

재설계 원칙.
  - 표 행과 각주 행을 문법으로 분리한다. 'AEGL 3' 뒤에 하이픈이 오면 각주다.
  - 표 행은 숫자를 긁지 않고 셀을 왼쪽부터 순서대로 식별한다. 셀은 숫자,
    NR, 'see below', 단독 각주표시(*) 중 하나다. 셀이 정확히 5개가 아니면
    값을 쓰지 않고 그 행을 보고한다. 칼럼이 밀리면 반드시 걸린다.
  - 각주의 10분 값은 별도 경로로 읽는다. 등호 뒤 첫 숫자만 취하고,
    등호 앞의 지속시간 라벨은 구조적으로 닿을 수 없다.
"""
import re

# 셀 식별자. 순서가 중요하다. see below -> 단독 각주표시 -> 숫자.
CELL = re.compile(r"""
    (?P<nr>\b(?:NR|ND)\b)
  | (?P<below>\*{0,3}\s*see\s+below\b)
  | (?P<mark>(?<![\d])\*{1,3}(?=\s|$))
  | (?P<num>\*{0,3}\s*(?:\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?))
""", re.I | re.X)

# 표 행: 'AEGL 3  ...' — 하이픈이 뒤따르면 각주이므로 제외한다.
ROW = re.compile(r'^AEGL\s*(?P<lvl>[123])\s+(?!\s*[-–—])(?P<body>.*)$')

# 각주 행: 'AEGL 3 - 10 min = ** 33,000 ppm'
FOOT = re.compile(
    r'^AEGL\s*[-–]?\s*(?P<lvl>[123])\s*[-–—]\s*(?P<spec>[^=]*)=\s*'
    r'(?P<rest>.*)$', re.I)

# 산문형 각주: '** The 10-min AEGL-3 value of 12,000 ppm (42,000 mg/m3) is ...'
PROSE = re.compile(
    r'(?P<dur>10)[-\s]*min\s+AEGL[-\s]*(?P<lvl>[123])\s+value\s+of\s+'
    r'(?P<val>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)\s*(?P<unit>ppm|mg/m3)',
    re.I)

# 각주 spec 이 10분을 포함하는가. '10 min', '10 min/30 min', '10 min/.../8 hr'
TEN_MIN = re.compile(r'(?<![\d])10\s*min', re.I)

_NUMONLY = re.compile(r'^(?:\s*(?:NR|\*{0,3}\s*\d[\d,.]*))+\s*$', re.I)


def _num(tok):
    return float(tok.replace('*', '').replace(',', '').strip())


def cells(body, strict=True):
    """표 행 본문을 셀 리스트로 만든다. 숫자 셀은 float, 그 외는 None.

    strict 면 첫 셀이 본문 선두에서 시작하지 않는 줄을 산문으로 보고 거부한다.
    'AEGL 2 is higher than 1/10 of the Lower Exposure Limit (LEL = 62,000 ...'
    같은 설명문이 표 행으로 인식되어 값을 덮어쓰는 것을 막는다.
    """
    ms = list(CELL.finditer(body))
    if not ms:
        return []
    if strict and body[:ms[0].start()].strip():
        return []
    out = []
    for m in ms:
        if m.group('num'):
            out.append(_num(m.group('num')))
        else:
            out.append(None)
    # 선두의 단독 각주표시는 첫 값에 붙은 표시일 수 있다.
    # 'AEGL 3 ** 10,000 5,200* 3,700* 1,800* 1,400*' 는 표시 1개 + 숫자 5개 = 6셀,
    # 이때 표시는 자리지킴이 아니라 10,000 에 붙은 것이므로 떼어낸다.
    # 반면 'AEGL 3 ** *7500 *5300 *2700 *1900' 는 표시 1개 + 숫자 4개 = 5셀이고,
    # 표시가 10분 칸 자체를 대신하므로 남긴다 (값은 각주에 있다).
    if len(out) == 6 and out[0] is None and ms[0].group('mark') and \
            all(v is not None for v in out[1:]):
        out = out[1:]
    return out


def numeric_only(line):
    """숫자와 NR 만으로 이루어진 줄인가 (pdfplumber 가 쪼갠 행의 나머지)."""
    return bool(line.strip()) and bool(_NUMONLY.match(line))


def parse_row(line, prev_orphan=None):
    """표 행을 파싱한다.

    셀이 5개면 그대로 쓴다. 5개 미만이고 바로 앞에 숫자만 있는 고아 줄이
    있으면, pdfplumber 가 행을 쪼갠 경우로 보고 합친다. 이때 10분 값은
    'AEGL N' 줄의 첫 셀이다 (쪼개진 나머지는 뒤쪽 지속시간이다).
    반환: (level, values5, status)
    """
    m = ROW.match(line)
    if not m:
        return None
    lvl, body = int(m.group('lvl')), m.group('body')
    # 기관이 병기하는 mg/m3 등가값은 대괄호와 소괄호 두 형식으로 나타난다.
    #   'AEGL 3  0.053 [0.38]  0.027 [0.19] ...'            (작용제)
    #   'AEGL 3  17 (64)  13 (47)  10 (38) ...'             (사산화이질소)
    # 둘 다 제거해야 셀이 5개가 된다. 환산 검증용 값은 _bnums 가 따로 읽는다.
    body_nb = re.sub(r'\[[^\]]*\]|\([^)]*\)', ' ', body)
    c = cells(body_nb)
    if len(c) == 5:
        return lvl, c, 'ok'
    if len(c) < 5 and prev_orphan is not None:
        extra = cells(prev_orphan)
        if len(c) + len(extra) == 5:
            return lvl, c + extra, 'merged'
    return lvl, None, f'cells={len(c)}'


def parse_footnote(line):
    """각주에서 10분 값을 읽는다. 반환: (level, value, unit_hint) 또는 None.

    등호 뒤 첫 숫자만 취한다. 세미콜론 뒤는 다른 지속시간이므로 자른다.
    등호 앞의 '10 min/30 min/60 min/4 hr/8 hr' 라벨에는 구조적으로 닿지 않는다.
    """
    m = FOOT.match(line)
    if m:
        if not TEN_MIN.search(m.group('spec')):
            return None
        rest = m.group('rest').split(';')[0]
        g = re.search(r'(\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)', rest)
        if not g:
            return None
        u = re.search(r'(ppm|mg/m3)', rest, re.I)
        return int(m.group('lvl')), _num(g.group(1)), (u.group(1).lower() if u else None)
    p = PROSE.search(line)
    if p:
        return int(p.group('lvl')), _num(p.group('val')), p.group('unit').lower()
    return None


# --------------------------------------------------------------- 자체 시험
CASES_ROW = [
    # (본문, 기대 셀)
    ('12,000          8,500             6,900         4,900   2,100',
     [12000.0, 8500.0, 6900.0, 4900.0, 2100.0]),
    ('0.52 0.17 0.086 NR NR', [0.52, 0.17, 0.086, None, None]),
    ('***27,000 ***27,000 ***22,000 **14,000 *6,800',
     [27000.0, 27000.0, 22000.0, 14000.0, 6800.0]),
    ('** *7500 *5300 *2700 *1900', [None, 7500.0, 5300.0, 2700.0, 1900.0]),
    ('**see below 8,600* 5,700* 2,500 1,700',
     [None, 8600.0, 5700.0, 2500.0, 1700.0]),
    ('** See below 14,000* 7,200* 2,400 1,600',
     [None, 14000.0, 7200.0, 2400.0, 1600.0]),
    ('***see below ***see below ***see below ***see below ***see below',
     [None, None, None, None, None]),
    ('NR NR NR NR NR', [None, None, None, None, None]),
    ('*4,000 *2,900 *2,900 *2,900 *2,900',
     [4000.0, 2900.0, 2900.0, 2900.0, 2900.0]),
    ('1300 920 730 460 300', [1300.0, 920.0, 730.0, 460.0, 300.0]),
    ('64               45                35               8.9             4.4',
     [64.0, 45.0, 35.0, 8.9, 4.4]),
    # 선두 각주표시 + 숫자 5개 = 표시를 떼어낸다 (톨루엔)
    ('** 10,000 5,200* 3,700* 1,800* 1,400*',
     [10000.0, 5200.0, 3700.0, 1800.0, 1400.0]),
    ('ND ND ND ND ND', [None, None, None, None, None]),
]

CASES_ROW_PAREN = [
    # 소괄호 mg/m3 병기 (사산화이질소)
    ('AEGL 3 17 (64) 13 (47) 10 (38) 7.0 (26) 5.7 (21)',
     [17.0, 13.0, 10.0, 7.0, 5.7]),
    ('AEGL 3 0.053 [0.38] 0.027 [0.19] 0.018 [0.13] 0.0098 [0.070] 0.0071 [0.051]',
     [0.053, 0.027, 0.018, 0.0098, 0.0071]),
]

CASES_PROSE_ROW = [
    # 산문은 표 행으로 쓰지 않는다
    'AEGL 2 is higher than 1/10 of the Lower Exposure Limit (LEL = 62,000 mg',
]

CASES_FOOT = [
    ('AEGL 3 - 10 min/30 min/60 min/4 hr/8 hr = *** 33,000 ppm', (3, 33000.0)),
    ('AEGL 3 - 10 min/30 min = ** 10,000 ppm', (3, 10000.0)),
    ('AEGL 3 - 10 min = *** 77,000 ppm; 30 min/ 60 min/ 4hr/ 8hr = *** 53,000 ppm',
     (3, 77000.0)),
    ('AEGL 3 - 10 min = ** 40,000 ppm', (3, 40000.0)),
    ('AEGL 3 - 10 min = ** 7,200 ppm', (3, 7200.0)),
    ('AEGL 3 - 10 min = ** 13,000 ppm', (3, 13000.0)),
    ('AEGL 2 - 10 min = *** 24,000 ppm; 30 min/ 60 min/ 4hr/ 8hr = ** 17,000 ppm',
     (2, 24000.0)),
    ('AEGL 1 - 10 min = ** 10,000 ppm', (1, 10000.0)),
    ('** The 10-min AEGL-3 value of 12,000 ppm (42,000 mg/m3) is higher than the '
     'lower explosive limit of n-', (3, 12000.0)),
    # 10분을 포함하지 않는 각주는 읽지 않는다
    ('AEGL 3 - 30 min/ 60 min/ 4hr/ 8hr = *** 53,000 ppm', None),
]

# 표 행으로 오인하면 안 되는 줄
CASES_NOTROW = [
    'AEGL 3 - 10 min = ** 40,000 ppm',
    'AEGL 2 - 10 min/30 min/60 min/4 hr/8 hr = ** 17,000 ppm',
]

if __name__ == '__main__':
    bad = 0
    for body, want in CASES_ROW:
        got = cells(body)
        if got != want:
            bad += 1
            print(f'  FAIL row  {body[:48]!r}\n        기대 {want}\n        실제 {got}')
    for line, want in CASES_FOOT:
        r = parse_footnote(line)
        got = (r[0], r[1]) if r else None
        if got != want:
            bad += 1
            print(f'  FAIL foot {line[:56]!r}\n        기대 {want}  실제 {got}')
    for line in CASES_NOTROW:
        if ROW.match(line):
            bad += 1
            print(f'  FAIL 각주를 표 행으로 인식 {line[:56]!r}')
    for line, want in CASES_ROW_PAREN:
        r = parse_row(line)
        got = r[1] if r else None
        if got != want:
            bad += 1
            print(f'  FAIL paren {line[:48]!r}\n        기대 {want}  실제 {got}')
    for line in CASES_PROSE_ROW:
        r = parse_row(line)
        if r is not None and r[1] is not None:
            bad += 1
            print(f'  FAIL 산문을 표 행으로 사용 {line[:56]!r} -> {r[1]}')
    total = (len(CASES_ROW) + len(CASES_FOOT) + len(CASES_NOTROW)
             + len(CASES_ROW_PAREN) + len(CASES_PROSE_ROW))
    print(f'[aegl_table] 자체시험 {total}건, 실패 {bad}건')
    raise SystemExit(1 if bad else 0)
