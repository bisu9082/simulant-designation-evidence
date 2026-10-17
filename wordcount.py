#!/usr/bin/env python3
"""투고본 분량 측정 정본 (2026-10-17 신설).

왜 필요한가. 분량을 그때그때 손으로 세다가 기준이 바뀌었다. 2026-10-16 에
'본문 8,025 단어'로 기록했는데 같은 원고를 다시 세니 8,448 이 나왔다. 차이는
원고가 늘어난 것이 아니라 lineno 패키지가 넣은 줄번호를 한 번은 세고 한 번은
빼서 생긴 것이다. 기준을 파일로 고정한다.

규칙.
  - PDF 렌더 결과를 센다. .tex 를 세면 명령어와 주석이 섞인다.
  - 순수 정수 토큰은 제외한다 (lineno 줄번호).
  - 본문은 '1. Introduction' 부터 'References' 직전까지.
  - 초록은 별도로 센다.
  - 표·도판 캡션을 포함한 값과 제외한 값을 함께 보고한다. 저널이 어느
    기준을 쓰는지 공식 문서로 확인하지 못했으므로(robots 차단) 둘 다 둔다.

사용: python3 wordcount.py
"""
import re
import subprocess
import sys

TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9'\-\./$]*")


def words(s):
    return [t for t in TOKEN.findall(s) if not re.fullmatch(r'\d+', t)]


def pdftext(path):
    r = subprocess.run(['pdftotext', '-layout', path, '-'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f'pdftotext 실패: {path}')
    return r.stdout


def sections(body):
    """번호 붙은 절 제목으로 쪼갠다.

    lineno 가 넣은 줄번호가 각 줄 앞에 붙어 절번호처럼 보이므로 먼저
    제거한다. 이걸 빼먹으면 '72 DMMP, dichloromethane...' 같은 가짜 절이
    잡히고 프로파일이 쓸모없어진다.
    """
    body = re.sub(r'^\s*\d+\s+', '', body, flags=re.M)
    heads = list(re.finditer(r'^(\d+(?:\.\d+)?)\.?\s+([A-Z][^\n]{3,70})$',
                             body, re.M))
    out = []
    for i, m in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(body)
        out.append((f'{m.group(1)} {m.group(2).strip()}',
                    len(words(body[m.end():end]))))
    return out


def main():
    txt = pdftext('paper/main.pdf')
    a = txt.index('Chemical warfare agent detection is developed')
    b = txt.index('Keywords')
    n_abs = len(words(txt[a:b]))

    s = txt.index('1. Introduction')
    e = txt.index('References')
    body = txt[s:e]
    n_body = len(words(body))

    # 캡션 제외분. 캡션은 'Table N.' / 'Figure N.' 로 시작해 빈 줄까지.
    cap = sum(len(words(m.group(0))) for m in
              re.finditer(r'^\s*(?:Table|Figure|Fig\.)\s*\d+[.:][\s\S]*?\n\s*\n',
                          body, re.M))

    print(f'초록            {n_abs:>6} 단어')
    print(f'본문 (전체)     {n_body:>6} 단어')
    print(f'본문 (캡션 제외){n_body - cap:>6} 단어   [캡션 {cap}]')
    print()
    for name, n in sections(body):
        print(f'  {name[:52]:<54} {n:>5}')

    sup = pdftext('paper/supplementary.pdf')
    print(f'\nSupplementary   {len(words(sup)):>6} 단어 (전체)')


if __name__ == '__main__':
    main()
