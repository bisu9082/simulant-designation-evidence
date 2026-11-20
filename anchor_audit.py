"""본문 수치 중 어느 앵커의 캡처 그룹에도 들어가지 않는 자리를 센다.

왜 필요한가 (2026-11-18 신설).

verify_consistency.py 의 (B) 검사는 '키 단위' 로 회계한다. 한 키가 어딘가
한 번이라도 맞으면 그 양은 '검사되고 있다' 로 집계되고, 같은 양을 다른
어법으로 적은 자리(초록·결론·도판 캡션)는 조용히 빠진다. 실제로 0.745 는
초록·본문·Conclusions·SI 네 곳에서 말해지는데 사실 키조차 없었고, 초록의
0.028 · 0.360 · 0.0206 은 주입해도 값 불일치가 0 건이었다.

이 스크립트는 회계 단위를 키에서 '자리' 로 바꿔, main.tex 본문의 수치
토큰마다 그것을 덮는 캡처 그룹이 있는지 묻는다. 덮이지 않은 자리는
(a) 앵커를 달거나 (b) DECLARED_CONSTANTS 에 사유와 함께 올린다.

verify_consistency 의 body 구성(주석 제거, \\begin{abstract} 이후, 공백
정규화본 병행)을 그대로 재현한다. 정규화본의 오프셋은 색인으로 원문에
되돌려 매핑하므로, 줄바꿈을 넘는 패턴도 정확히 덮인 것으로 센다.
"""
import re
import sys

import verify_consistency as V

# 본문이 정당하게 담는 상수. 값이 아니라 '그 자리' 를 면제하는 것이 아니라
# 토큰 값을 면제한다. 사유 없이 늘리지 말 것.
DECLARED_CONSTANTS = {
    '0.05': '명목 유의수준',
    '0.95': '신뢰수준 / 상관필터 역치',
    '2.5': '부트스트랩 하위 백분위 라벨',
    '97.5': '부트스트랩 상위 백분위 라벨',
    '298.15': '온도 (K)',
    '25': '온도 (섭씨)',
    '0.79': '\\textwidth 배수 (조판)',
    '1': '단위·서수·분자', '2': '단위·서수', '3': '단위·서수',
    '0': '영',
    }

# \, 는 LaTeX 천단위 구분자다. 토큰 안에 넣어 10\,000 을 하나로 읽는다
# (치환으로 지우면 오프셋이 어긋나 덮임 판정이 틀어진다).
NUMTOK = re.compile(r'(?<![\d.A-Za-z\\])\d(?:[\d,]|\\,)*(?:\.\d+)?')

# 수치가 아니라 이름인 자리. 토큰화 전에 같은 길이의 'X' 로 가린다
# (오프셋이 어긋나면 덮임 판정이 틀어지므로 길이를 보존한다).
IDENT = [
    r'TOP\s*\d+-\d+-\d+[A-Z0-9]*',        # 시험절차 번호
    r'acs-submission-\d+',                  # 예치 태그 이름
    r'quotations\\?_\d+\.csv',               # 인용 대장 파일명
    r'CIEDE\d+',                            # 색차식 이름
    r'_\{\d+\}',                            # 수식 첨자 (\Delta E_{00})
    r'(?<=\d),(?=\d-)',                     # 화학 위치번호 2,4-
    r'(?<=[A-Za-z])~\d+',                    # lewisite~1
    r'AD\d+',                               # DTIC 접수번호
    r'US\s*\d[\d,]*',                       # 특허번호
    # 2026-11-18. 달력 날짜는 측정값이 아니라 이름이다. 앵커를 달 대상이
    # 아니므로 토큰화 전에 가린다 (연도까지 한 덩어리로).
    r'\d{1,2}(?:,\s*\d{1,2})*(?:\s+and\s+(?:again\s+on\s+)?\d{1,2})?'
    r'\s+(?:January|February|March|April|May|June|July|August|September'
    r'|October|November|December)(?:\s+\d{4})?',
    r'(?:January|February|March|April|May|June|July|August|September'
    r'|October|November|December)\s+\d{4}',
]


def mask_identifiers(t):
    for pat in IDENT:
        t = re.sub(pat, lambda m: 'X' * len(m.group(0)), t)
    return t


def flatten(t):
    """공백을 한 칸으로 줄이고, 각 출력 문자의 원문 오프셋을 함께 돌려준다."""
    out, idx, i, n = [], [], 0, len(t)
    while i < n:
        if t[i].isspace():
            j = i
            while j < n and t[j].isspace():
                j += 1
            out.append(' ')
            idx.append(i)
            i = j
        else:
            out.append(t[i])
            idx.append(i)
            i += 1
    return ''.join(out), idx


def main_body():
    tex = V.strip_comments(open(V.TEX).read()) \
        if hasattr(V, 'strip_comments') else None
    if tex is None:                     # strip_comments 가 내부 함수인 판본
        src = open(V.TEX).read()
        keep = []
        for line in src.split('\n'):
            if line.lstrip().startswith('%'):
                continue
            keep.append(re.sub(r'(?<!\\)%.*$', '', line))
        tex = '\n'.join(keep)
    return tex[tex.index('\\begin{abstract}'):]


def covered_spans(body):
    flat, idx = flatten(body)
    spans = []
    for key, pats, _tol in V.CHECKS:
        for pat in pats:
            try:
                for m in re.finditer(pat, body):
                    spans.append((m.start(1), m.end(1), key))
                for m in re.finditer(pat, flat):
                    a, b = m.span(1)
                    spans.append((idx[a], idx[b - 1] + 1, key))
            except re.error:
                print(f'[audit] 정규식 오류 {key}: {pat!r}', file=sys.stderr)
    return spans


def run():
    body = main_body()
    spans = covered_spans(body)
    masked = mask_identifiers(body)
    assert len(masked) == len(body), '가림이 길이를 바꿨다'
    uncovered = []
    for m in NUMTOK.finditer(masked):
        a, b = m.span()
        if any(s <= a and b <= e for s, e, _ in spans):
            continue
        tok = m.group(0).rstrip(',')
        if not tok:
            continue
        a, b = a, a + len(tok)
        if tok in DECLARED_CONSTANTS:
            continue
        ctx = ' '.join(body[max(0, a - 55):b + 25].split())
        uncovered.append((tok, ctx, body[:a].count('\n') + 1))
    return uncovered


if __name__ == '__main__':
    rows = run()
    print(f'[audit] main.tex 본문 미앵커 수치 {len(rows)}곳')
    for tok, ctx, ln in rows:
        print(f'  b{ln:<5d} {tok:>10}  …{ctx}')
    raise SystemExit(0)
