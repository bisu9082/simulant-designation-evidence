#!/usr/bin/env python3
"""검증기가 잡은 불일치를 원고에 역으로 반영한다 (2026-10-24 신설).

왜 필요한가. 정확종속 규칙을 넣고 파이프라인을 다시 돌리자 28종 패널 계열에서
55건이 어긋났다. 손으로 고치면 반드시 빠뜨린다. 검증기는 이미 각 수치가 어느
정규식으로 본문 어디에 걸리는지 알고 있으므로, 같은 정규식으로 찾아 포획군만
새 값으로 바꾼다.

원칙 두 가지.
 1. 본문에 적힌 자릿수를 그대로 유지한다. 0.028 자리에는 0.031 을 넣지
    0.0308824 를 넣지 않는다.
 2. 영어 수사(three, four)는 숫자와 섞지 않는다. 사전으로만 바꾼다.

--apply 없이 부르면 바꿀 내용만 출력한다.
"""
import re
import sys

import verify_consistency as V

FILES = ['paper/main.tex', 'paper/supplementary.tex',
         'paper/cover_letter.tex', 'paper/highlights.tex']

WORD = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five',
        6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten',
        11: 'eleven', 12: 'twelve', 13: 'thirteen', 14: 'fourteen',
        15: 'fifteen', 16: 'sixteen', 17: 'seventeen', 18: 'eighteen',
        19: 'nineteen', 20: 'twenty'}
WORD_INV = {v: k for k, v in WORD.items()}


def fmt_like(old, new):
    """old 의 표기 관습을 따라 new 를 찍는다."""
    o = old.strip()
    if o.lower() in WORD_INV:
        iv = int(round(new))
        return WORD[iv] if iv in WORD else str(iv)
    if o.lower() in ('zero',):
        return WORD.get(int(round(new)), str(int(round(new))))
    neg = o.startswith(('-', '+'))
    sign = o[0] if neg else ''
    body = o[1:] if neg else o
    body = body.replace(',', '').replace('\\,', '')
    if '.' in body:
        dec = len(body.split('.')[1])
        s = f'{abs(new):.{dec}f}'
    else:
        s = f'{abs(new):.0f}'
    if sign == '+':
        return '+' + s
    if sign == '-':
        return '-' + s
    return s


def main():
    apply = '--apply' in sys.argv
    f = V.load_facts()
    texts = {p: open(p, encoding='utf8').read() for p in FILES}
    changes, skipped = [], []

    for key, pats, tol in V.CHECKS:
        if key not in f:
            continue
        want = float(f[key])
        for path in FILES:
            for pat in pats:
                for m in list(re.finditer(pat, texts[path])):
                    old = m.group(1)
                    try:
                        got = V.tonum(old)
                    except Exception:
                        continue
                    if abs(got - want) <= tol:
                        continue
                    new = fmt_like(old, want)
                    if V.tonum(new) == got:
                        skipped.append((path, key, old, want))
                        continue
                    a, b = m.span(1)
                    texts[path] = texts[path][:a] + new + texts[path][b:]
                    ctx = ' '.join(m.group(0).split())[:58]
                    changes.append((path, key, old, new, ctx))

    for path, key, old, new, ctx in changes:
        print(f'{path.split("/")[-1]:20s} {key:22s} {old:>10s} -> {new:<10s} [{ctx}]')
    for path, key, old, want in skipped:
        print(f'[자릿수 부족] {path.split("/")[-1]} {key}: 본문 {old}, 산출물 {want:g} '
              f'— 표기 자릿수로는 구분되지 않는다. 손으로 볼 것.')
    print(f'\n바꿀 곳 {len(changes)}건, 자릿수 문제 {len(skipped)}건')
    if apply:
        for p, t in texts.items():
            open(p, 'w', encoding='utf8').write(t)
        print('[applied]')
    else:
        print('[미적용] --apply 를 붙이면 실제로 쓴다.')


if __name__ == '__main__':
    main()
