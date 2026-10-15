#!/usr/bin/env python3
"""AI 문체 스캐너 (2026-10-13 신설).

Ku 지시(2026-08-15 firenose2, 2026-10-01 em-dash 절대금지로 강화)의 금지 목록을
원고·커버레터·하이라이트·도판 주석 전부에 적용한다. 종전에는 본문만 grep 으로
훑었고, 그 바람에 도판 스크립트 안의 em-dash 네 건이 살아남았다.

검사 대상
  1) 금지 어휘
  2) 금지 구조 (문단 첫머리 First/Second/Furthermore 기계적 전개, 삼항 나열,
     문단 끝 요약문 반복)
  3) em-dash (U+2014, U+2013 사이 공백형, LaTeX ---)
LaTeX 주석행(%)과 코드 주석은 출력물이 아니므로 제외한다.
"""
import glob
import re
import sys

VOCAB = [
    r'\bdelve\b', r'\bnotably\b', r'\bcrucially\b', r'\bremarkably\b',
    r'\bintriguingly\b', r'it is worth noting', r'\bworth noting\b',
    r'plays? an? (?:crucial|pivotal|key|vital) role', r'paves? the way',
    r'sheds? light on', r'\bholistic\b', r'\bleverag(?:e|es|ed|ing)\b',
    r'\blandscape\b', r'\bcomprehensive\b', r'\brobust\b',
    r'\bmoreover\b', r'\bfurthermore\b', r'\bunderscores?\b',
    r'\bhighlights? the\b', r'\bdelves? into\b', r'\bin conclusion\b',
    r'\bit should be noted\b', r'\bnoteworthy\b',
]
DASH = [('em-dash U+2014', '—'), ('en-dash 공백형', ' – '),
        ('LaTeX ---', '---')]
STRUCT = [
    (r'^(First|Second|Third|Finally|Furthermore|Moreover|Additionally)\b',
     '문단 첫머리 기계적 전개어'),
]


def strip_comments(text, kind):
    out = []
    for ln in text.split('\n'):
        if kind == 'tex':
            # LaTeX 주석 제거. \% 는 주석이 아니다.
            ln = re.sub(r'(?<!\\)%.*$', '', ln)
        out.append(ln)
    return '\n'.join(out)


def scan(path, kind):
    raw = open(path, encoding='utf-8').read()
    body = strip_comments(raw, kind)
    hits = []
    for pat in VOCAB:
        for m in re.finditer(pat, body, re.I):
            ln = body[:m.start()].count('\n') + 1
            ctx = ' '.join(body[max(0, m.start() - 45):m.end() + 35].split())
            hits.append(('어휘', pat, ln, ctx))
    for lab, ch in DASH:
        for m in re.finditer(re.escape(ch), body):
            ln = body[:m.start()].count('\n') + 1
            ctx = ' '.join(body[max(0, m.start() - 45):m.end() + 35].split())
            hits.append(('대시', lab, ln, ctx))
    for pat, lab in STRUCT:
        for m in re.finditer(pat, body, re.M):
            ln = body[:m.start()].count('\n') + 1
            ctx = ' '.join(body[m.start():m.start() + 70].split())
            hits.append(('구조', lab, ln, ctx))
    return hits


TARGETS = [('paper/main.tex', 'tex'), ('paper/supplementary.tex', 'tex'),
           ('paper/cover_letter.tex', 'tex'),
           ('paper/highlights.tex', 'tex'),
           ('outputs/figure_annotations.json', 'json')]

if __name__ == '__main__':
    total = 0
    for path, kind in TARGETS:
        try:
            hits = scan(path, kind)
        except FileNotFoundError:
            print(f'[skip] {path} 없음')
            continue
        print(f'\n=== {path}  적발 {len(hits)}건')
        for cat, pat, ln, ctx in sorted(hits, key=lambda h: h[2]):
            print(f'  {cat} {pat:<34} {ln:>5}행  {ctx[:88]}')
        total += len(hits)
    # 도판 스크립트의 렌더링 문자열도 본다 (주석 제외)
    fig = 0
    for f in sorted(glob.glob('make_fig_*.py') + glob.glob('make_ga.py')):
        src = open(f, encoding='utf-8').read()
        src = re.sub(r'""".*?"""', '', src, flags=re.S)
        src = re.sub(r'^\s*#.*$', '', src, flags=re.M)
        for lab, ch in DASH:
            for m in re.finditer(re.escape(ch), src):
                ln = src[:m.start()].count('\n') + 1
                print(f'  [도판] {f}:{ln} {lab}  '
                      f'{" ".join(src[max(0,m.start()-40):m.end()+30].split())[:80]}')
                fig += 1
    print(f'\n[ai-style] 원고·부속 {total}건, 도판 스크립트 {fig}건')
    sys.exit(1 if (total or fig) else 0)
