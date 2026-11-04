#!/usr/bin/env python3
"""EPA AEGL 표 수치 파서 — 천단위 쉼표 안전.

2026-10-11 사고 교정. 종전 정규식 r'\\b(\\d+\\.?\\d*)\\b' 은 단어 경계가
쉼표에서 끊겨 "12,000" 을 12 와 000 으로 쪼갰다. AEGL-3 10분 값에 쉼표가
있는 화합물 22종이 전부 1000배 낮게 기록됐고, 한 행의 5개 노출시간이
모두 어긋났다 ("12,000 8,500 6,900 4,900 2,100" -> [12, 0, 8, 500, 6]).

여기서는 쉼표 묶음을 하나의 토큰으로 먼저 잡는다.
"""
import re

# 1,234 / 1,234,567 / 12.5 / 0.059 / 1.2e-3 / 1,400* 를 모두 1개 토큰으로.
# 앞뒤에 숫자·쉼표·소수점이 붙어 있으면 토큰 경계가 아니다.
NUM = re.compile(
    r'(?<![\d,.])'
    r'(\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)'
    r'(?:[eE][-+]?\d+)?'
    r'(?![\d,]|\.\d)'
)


def nums(text):
    """문자열에서 수치 토큰을 순서대로 추출해 float 리스트로 반환."""
    out = []
    for m in NUM.finditer(text):
        tok = m.group(0).replace(',', '')
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def bracket_nums(text):
    """대괄호 안의 환산값 [12,000] 추출 (쉼표 허용)."""
    out = []
    for b in re.findall(r'\[([^\]]+)\]', text):
        out += nums(b)
    return out


# --------------------------------------------------------------- 자체 검증
# 종전 버그를 재현하는 사례를 반드시 포함한다. 쉼표 없는 작은 값만으로
# 검증하면 이 실패 모드를 구조적으로 볼 수 없다 (2026-10-11 사고 원인).
CASES = [
    # (입력, 기대 출력)
    ('12,000          8,500             6,900         4,900   2,100',
     [12000.0, 8500.0, 6900.0, 4900.0, 2100.0]),          # methylene chloride AEGL-3
    ('** 10,000        5,200*                   3,700*              1,800*     1,400*',
     [10000.0, 5200.0, 3700.0, 1800.0, 1400.0]),           # toluene AEGL-3
    ('4,000             4,000              3,200     2,000   1,600',
     [4000.0, 4000.0, 3200.0, 2000.0, 1600.0]),            # chloroform AEGL-3
    ('0.064            0.032             0.022     0.012    0.0087',
     [0.064, 0.032, 0.022, 0.012, 0.0087]),                # 쉼표 없는 소수
    ('1,400*              760                560                 310       250',
     [1400.0, 760.0, 560.0, 310.0, 250.0]),                # toluene AEGL-2
    ('27,000           27,000            27,000    27,000  27,000',
     [27000.0] * 5),                                       # HFC 134A
    ('1,100            370               190       48      24',
     [1100.0, 370.0, 190.0, 48.0, 24.0]),                  # 혼합
    ('NR               NR                NR        NR      NR', []),
    ('290             230               200           NR      NR',
     [290.0, 230.0, 200.0]),
    ('Lower Explosive Limit (LEL) = 14,000 ppm', [14000.0]),
]

BRACKET_CASES = [
    ('0.064 [0.39]', [0.39]),
    ('12 [12,000]', [12000.0]),
]


def self_test(verbose=True):
    bad = []
    for src, want in CASES:
        got = nums(src)
        if got != want:
            bad.append(('nums', src[:46], want, got))
    for src, want in BRACKET_CASES:
        got = bracket_nums(src)
        if got != want:
            bad.append(('bracket', src[:46], want, got))
    if verbose:
        print(f'[aegl_num] 자체 검증 {len(CASES)+len(BRACKET_CASES)}건 중 '
              f'{len(bad)}건 실패')
        for kind, src, want, got in bad:
            print(f'   FAIL [{kind}] {src!r}\n        기대 {want}\n        실제 {got}')
    return not bad


if __name__ == '__main__':
    import sys
    sys.exit(0 if self_test() else 1)
