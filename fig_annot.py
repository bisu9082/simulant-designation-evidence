#!/usr/bin/env python3
"""도판 주석 문자열을 기록해 검증 대상에 넣는다 (2026-10-11 신설).

필요해진 이유. verify_consistency.py 는 main.tex 만 읽었고, 도판 안의 글자는
PNG 픽셀이라 읽히지 않았다. 그래서 Figure 2 의 주석이 'the least acutely toxic
of the thirteen' 이라는 구판 문자열을 들고 있는데도 본문 검사와 빌드 검사가
모두 통과했다. 본문은 'not the safest' 라고 쓰고 같은 패널의 막대는 메탄올과
n-헥산을 위에 그리는 상태였다. build_paper.sh 는 outputs/ 와 paper/ 의 해시만
비교하므로 그림이 최신이라는 것만 보장했지, 그림의 '내용'이 본문과 맞는지는
보지 않았다.

이제 도판 스크립트가 자기 주석을 여기에 등록하고, verify_consistency.py 가
본문과 같은 규칙으로 그 문자열을 검사한다.
"""
import json
import os

PATH = 'outputs/figure_annotations.json'

# 2026-10-20. 등재부가 누적만 하고 지우지 않아, 폐기된 도판(Fig_panel,
# Fig_qsar_failure)의 주석이 파일에 남아 검증 대상에 계속 들어갔다. 사라진
# 패널의 구판 문자열이 본문 검사를 대신 통과시킬 수 있는 구멍이므로, 한
# 프로세스에서 어떤 도판을 처음 기록할 때 그 도판의 기존 항목을 비운다.
_STARTED = set()


def record(figure, key, text):
    """도판 하나의 주석 한 덩어리를 기록한다."""
    d = {}
    if os.path.exists(PATH):
        try:
            d = json.load(open(PATH))
        except Exception:
            d = {}
    if figure not in _STARTED:
        _STARTED.add(figure)
        d[figure] = {}
    d.setdefault(figure, {})[key] = text
    json.dump(d, open(PATH, 'w'), indent=1, ensure_ascii=False)
    return text


def load():
    if not os.path.exists(PATH):
        return {}
    return json.load(open(PATH))


def all_text():
    """등록된 모든 주석을 한 덩어리 문자열로 돌려준다."""
    d = load()
    out = []
    for fig, items in sorted(d.items()):
        for k, v in sorted(items.items()):
            out.append(f'[{fig}.{k}] {v}')
    return '\n'.join(out)
