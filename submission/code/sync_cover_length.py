#!/usr/bin/env python3
"""커버레터의 분량 숫자를 length_audit.json 에서 되쓴다 (2026-11-14 신설).

왜 필요한가. 이 세션에서 다섯 번 같은 결함이 났다. 커버레터는 분량을 네 개의
숫자로 적고 검증기가 그 넷을 outputs/length_audit.json 과 대조한다. 그런데
그 JSON 은 빌드된 PDF 를 읽어야 만들어지므로 build_paper.sh 의 맨 끝에서
생성되고, 커버레터는 그보다 먼저 [6/6] 에서 조판된다. 즉 커버레터는 구조적으로
항상 한 빌드 뒤처진다. 본문을 한 단어라도 고치면 커버레터가 구판이 되고,
손으로 고치면 다음 수정에서 또 틀린다.

이 프로젝트의 반복 결함은 계산 오류가 아니라 전파 실패이고, 이것이 그
전형이다. 숫자를 손으로 들고 있지 않고 빌드에 묶는다.

커버레터 본문은 표현을 유지하고 숫자만 바꾼다. 바뀐 것이 있으면 종료코드 10
을 돌려주며, build_paper.sh 가 그때만 커버레터를 다시 조판한다.

본문 단어수는 main.pdf 에서 세므로 커버레터를 고쳐도 감사값이 변하지 않는다.
그래서 반복이 필요 없다 (고정점이 한 번에 잡힌다).
"""
import json
import re
import sys

TEX = 'paper/cover_letter.tex'
AUDIT = 'outputs/length_audit.json'

# (사실 키, 정규식, 치환될 그룹의 서식)
RULES = [
    ('text_words',
     r'(manuscript stands at\s*\n?)([\d,]+)( words of text)', 'comma'),
    ('ref_words',
     r'(of which\s*\n?)([\d,]+)( words are references)', 'comma'),
    # 2026-11-16. 커버레터가 두 번째 셈법(맨정수 포함)과 참고문헌 제외값도
    # 적는다. 규칙을 걸지 않으면 그 둘만 뒤처진다.
    ('text_words_tokens',
     r'(instead gives\s*\n?)([\d,]+)', 'comma'),
    ('text_words_norefs',
     r'(manuscript is\s*\n?)([\d,]+)( words)', 'comma'),
    # 2026-11-15. 8면 기준 3개 규칙(text_plus_fig / acs_pages / acs_overage)을
    # 제거했다. 편집부(Jie Fu, Managing Editor, 2026-11-15)가 지배 지침을
    # researcher-resources 판 10,000단어로 지정했고, 그 판에는 도판 환산이
    # 없다. 커버레터가 8면 환산치를 더는 적지 않으므로 되쓸 대상도 없다.
    # length_audit.json 은 계속 acs_pages 를 계산하므로, 기준이 다시 바뀌면
    # 규칙만 되살리면 된다.
]


def main():
    a = json.load(open(AUDIT))
    a['text_plus_fig'] = a['text_words'] + a['figure_words']
    # 8면 한도 초과분. 커버레터가 숫자로 적으므로 여기서 함께 되쓴다.
    a['acs_overage'] = a['acs_pages'] - 8.0
    src = open(TEX, encoding='utf-8').read()
    out, changed = src, []
    for key, pat, fmt in RULES:
        want = a[key]
        text = f'{want:,}' if fmt == 'comma' else f'{want:.2f}'
        m = re.search(pat, out)
        if m is None:
            print(f'[cover] 어법 못 찾음: {key} — 정규식을 고칠 것',
                  file=sys.stderr)
            return 1
        if m.group(2) != text:
            changed.append(f'{key}: {m.group(2)} -> {text}')
            out = out[:m.start(2)] + text + out[m.end(2):]
    if not changed:
        print('[cover] 분량 숫자 최신')
        return 0
    open(TEX, 'w', encoding='utf-8').write(out)
    for c in changed:
        print(f'[cover] {c}')
    return 10


if __name__ == '__main__':
    sys.exit(main())
