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
    # 2026-11-09. 줄바꿈 하이픈을 결합한다. 'mea-\nsurements' 가 TOKEN 에
    # 두 토큰으로 잡혀 본문이 과대 집계되고 있었다 (본문 51건, 초록 4건,
    # 참고문헌 21건). Word 가 세는 방식에 맞춘다.
    return re.sub(r'(\w)-\n\s*(\w)', r'\1\2', r.stdout)


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


def acs_sections():
    """원고의 \\section 제목을 정본으로 쓴다 (2026-11-06 재작성).

    종전 구현은 JHM 양식을 전제로 'Keywords' 와 '1. Introduction' 을 찾았다.
    achemso 본에는 둘 다 없어서 ValueError 로 죽었고, 9차 채점 때까지 분량
    예산 도구가 현행 원고에서 한 번도 돌지 않았다. 본문이 자체 예산 5,000
    단어를 53% 넘긴 것을 아무도 못 본 이유다. 양식에 의존하지 않도록
    제목 목록을 main.tex 에서 직접 읽는다.
    """
    src = open('paper/main.tex').read()
    return [m.group(1) for m in re.finditer(r'^\\section\*?\{([^}]+)\}', src, re.M)]


def main():
    txt = pdftext('paper/main.pdf')
    heads = acs_sections()
    if not heads:
        raise SystemExit('[wordcount] main.tex 에서 \\section 을 찾지 못했다.')

    first = txt.index(heads[0])
    a = txt.index('Abstract')
    n_abs = len(words(txt[a + len('Abstract'):first]))

    e = txt.index('References')
    body = txt[first:e]
    n_body = len(words(body))

    cap = sum(len(words(m.group(0))) for m in
              re.finditer(r'^\s*(?:Table|Figure|Fig\.)\s*\d+[.:][\s\S]*?\n\s*\n',
                          body, re.M))

    print(f'초록            {n_abs:>6} 단어')
    print(f'본문 (전체)     {n_body:>6} 단어')
    print(f'본문 (캡션 제외){n_body - cap:>6} 단어   [캡션 {cap}]')
    print()
    # 2026-11-08. 종전에는 body.index(h) 로 각 제목의 첫 출현을 썼다. 절
    # 제목과 같은 문구가 산문에 먼저 나오면 그쪽이 잡힌다. 실제로
    # 'Data and code availability' 가 Results 앞의 산문에서 먼저 걸려
    # 마지막 절에 8,258 단어가 붙었고, 섹션 합이 본문 총계를 넘었다.
    # 앞 절 뒤에서부터 순차로 찾아 순서를 강제한다.
    pos, _at = [], 0
    for h in heads:
        j = body.find(h, _at)
        if j < 0:
            continue
        pos.append((h, j))
        _at = j + len(h)
    for i, (h, j) in enumerate(pos):
        k = pos[i + 1][1] if i + 1 < len(pos) else len(body)
        print(f'  {h[:52]:<54} {len(words(body[j + len(h):k])):>5}')

    # ACS Article 환산 (2026-11-07 정정). 저자안내는 총 분량에 본문·그림·
    # 표·참고문헌을 전부 포함하고 TOC 그래픽만 제외한다. 종전 구현은
    # 거꾸로였다 — TOC 를 도판으로 세고 참고문헌을 빼서 면수를 과소·과대
    # 양쪽으로 틀리게 냈다. 1,000 단어 = 1면, 양단 도판 1장 = 500 단어.
    src = open('paper/main.tex').read()
    figs = [m.group(1) for m in
            re.finditer(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', src)]
    # 2026-11-09. 저자안내는 분량에 "all text, figures, tables and
    # references" 를 넣고 TOC 그래픽만 뺀다. 종전 구현은 표를 아예 세지
    # 않아 약 750 단어를 과소 집계했다. tabular 전문에 \textwidth 가
    # 있으면 양단(500), 없으면 단단(250)으로 본다.
    # 열 정의에 중첩 중괄호(@{} 나 p{0.27\textwidth})가 들어가므로
    # [^}]* 로는 첫 } 에서 끊긴다. 줄 끝까지 본다.
    tabs = []
    for _m in re.finditer(r'\\begin\{tabular\}(.*)', src):
        tabs.append(500 if 'textwidth' in _m.group(1) else 250)
    # TOC 그래픽은 분량 환산에서 제외된다. 파일명이 graphical_abstract
    # 로 바뀌어 'toc' 필터를 빠져나가 도판으로 세어지고 있었다.
    figs = [f for f in figs
            if 'toc' not in f.lower()
            and 'graphical_abstract' not in f.lower()]
    refs = len(words(txt[txt.index('References'):]))
    pages = ((n_body - cap) / 1000.0 + len(figs) * 0.5
             + sum(tabs) / 1000.0 + refs / 1000.0)
    print(f'\n  참고문헌      {refs:>6} 단어')
    print(f'  ACS 환산 {pages:.1f}면 = 본문 {(n_body-cap)/1000:.1f}'
          f' + 도판 {len(figs)}장 {len(figs)*0.5:.1f}'
          f' + 표 {len(tabs)}개 {sum(tabs)/1000:.2f}'
          f' + 참고문헌 {refs/1000:.1f}   (TOC 제외, Article 한도 8면)')

    # 2026-11-07. 커버레터가 분량을 숫자로 적으므로 그 숫자가 이 측정의
    # 결과여야 한다. verify_consistency.py 가 이 파일을 읽어 대조한다.
    import json as _j
    import os as _o
    _o.makedirs('outputs', exist_ok=True)
    _j.dump({'body_words': n_body - cap, 'caption_words': cap,
             'ref_words': refs, 'n_figures': len(figs),
             'n_tables': len(tabs), 'table_words': sum(tabs),
             'acs_pages': round(pages, 2)},
            open('outputs/length_audit.json', 'w'), indent=1)
    print('  -> outputs/length_audit.json')

    sup = pdftext('paper/supplementary.pdf')
    print(f'\nSupplementary   {len(words(sup)):>6} 단어 (전체)')


if __name__ == '__main__':
    main()
