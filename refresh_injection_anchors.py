#!/usr/bin/env python3
"""주입시험의 앵커 문자열을 현행 원고 수치로 맞춘다 (2026-10-24 신설).

주입시험은 (이름, 원문조각, 망친조각) 세 쪽으로 되어 있고, 원문조각이 본문에서
사라지면 그 시험은 조용히 SKIP 된다. 분석을 다시 돌려 수치가 바뀔 때마다 수십
건이 한꺼번에 SKIP 되는데, 그 상태는 '통과'로 보이지 그 시험이 꺼졌다고 보이지
않는다. 이 세션에서 그 경로로 세 번 빠뜨렸다.

여기서는 원문조각 안의 수를 현행 값으로 바꾸고, 망친조각 안의 수는 그대로 둔다.
망친 값이 우연히 현행 값과 같아지면 시험이 무의미해지므로 그 경우만 경고한다.
실제로 잡히는지는 test_verify_consistency.py 가 판정한다.

--apply 없이 부르면 바꿀 내용만 출력한다.
"""
import re
import sys

NUM = re.compile(r'(?<![\d.])([-+]?\d+(?:\.\d+)?)(?![\d.])')


def numbers(s):
    return [(m.group(1), m.span(1)) for m in NUM.finditer(s)]


def main():
    apply = '--apply' in sys.argv
    import verify_consistency as V
    import test_verify_consistency as T
    f = V.load_facts()
    src = open('test_verify_consistency.py', encoding='utf8').read()

    body = ''
    for p in ('paper/main.tex', 'paper/supplementary.tex',
              'paper/cover_letter.tex', 'paper/highlights.tex'):
        body += open(p, encoding='utf8').read() + '\n'

    cases = []
    for name in ('CASES', 'SUPP_CASES', 'MAIN_CASES'):
        cases += list(getattr(T, name, []))

    done, warn = [], []
    for case in cases:
        if len(case) < 3:
            continue
        label, old, new = case[0], case[-2], case[-1]
        if old in body:
            continue                       # 아직 살아 있는 앵커
        on, nn = numbers(old), numbers(new)
        if not on or len(on) != len(nn):
            warn.append((label, '원문조각과 망친조각의 수 개수가 다르다'))
            continue
        # 어느 자리가 망가진 자리인가
        bad = [i for i, (o, n) in enumerate(zip(on, nn)) if o[0] != n[0]]
        if len(bad) != 1:
            warn.append((label, '망가진 자리가 하나가 아니다'))
            continue
        bi = bad[0]
        # 수만 비워 둔 뼈대로 현행 본문에서 같은 문맥을 찾는다
        skel, last = '', 0
        for tok, (a_, b_) in on:
            skel += re.escape(old[last:a_]) + r'([-+]?[\d.]+)'
            last = b_
        skel += re.escape(old[last:])
        m = re.search(skel, body)
        if not m:
            warn.append((label, '본문에서 대응 문맥을 못 찾음'))
            continue
        new_old = m.group(0)
        cur_tokens = [m.group(i + 1) for i in range(len(on))]
        # 망친조각은 현행 문맥에 같은 자리를 같은 망친 값으로 바꿔 만든다
        cn = numbers(new_old)
        a_, b_ = cn[bi][1]
        corrupt = nn[bi][0]
        if corrupt == cur_tokens[bi]:
            warn.append((label, f'망친 값이 현행 값과 같아짐 ({corrupt})'))
            continue
        new_new = new_old[:a_] + corrupt + new_old[b_:]
        if old not in src or new not in src:
            warn.append((label, '시험 파일에서 조각을 못 찾음'))
            continue
        src = src.replace(old, new_old, 1).replace(new, new_new, 1)
        done.append((label, on[bi][0], cur_tokens[bi]))

    for label, o, n in done:
        print(f'  {o:>10s} -> {n:<10s}  {label}')
    for label, why in warn:
        print(f'  [수동] {label}: {why}')
    print(f'\n갱신 {len(done)}건, 수동 확인 {len(warn)}건')
    if apply:
        open('test_verify_consistency.py', 'w', encoding='utf8').write(src)
        print('[applied]')
    else:
        print('[미적용] --apply 를 붙이면 실제로 쓴다.')


if __name__ == '__main__':
    main()
