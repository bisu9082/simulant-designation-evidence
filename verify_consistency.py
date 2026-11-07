#!/usr/bin/env python3
"""본문 수치를 산출물과 대조한다 (2026-10-11 전면 재작성, 방향을 뒤집었다).

종전 설계의 결함. 핵심 수치 목록을 손으로 적어 두고 그 문자열이 main.tex 에
'존재하는지'만 확인했다. 존재 확인은 모순을 잡지 못한다. 그래서
  - 본문이 같은 양을 208 과 218 로 동시에 적어도 통과했고,
  - 같은 환산값을 41,683 과 41,685 로 적어도 통과했으며,
  - 19,531 과 19,530 도 통과했다.
GATE 8 3차에서 전부 지적됐다.

재설계. 산출물(JSON/CSV)에서 사실표를 만든 뒤, 본문에서 그 양을 가리키는
문맥을 모두 찾아 값을 추출하고 사실표와 비교한다. 두 가지를 본다.
  (A) 문맥 대조: 한 양을 가리키는 모든 등장이 산출물 값과 같은가.
      같은 양을 두 값으로 적으면 반드시 걸린다.
  (B) 근사 불일치: 본문의 수치 중 사실표의 어떤 값과 0.5% 안쪽으로
      가까우면서 같지는 않은 것. 반올림 분기(84.93 대 84.933)를 잡는다.
"""
import glob
import json
import os
import re
import numpy as _np
import sys

import pandas as pd

# 2026-11-04 (결정 J-C). 1차 본문을 ACS본으로 옮긴다.
# 종전에는 paper/main.tex 가 TEX 였고 paper_acs/ 가 추가본이었다. JHM
# 평행본은 archive/jhm_v1/ 로 동결했으므로 축을 뒤집는다. 동결 사유와
# 되살리는 조건은 archive/jhm_v1/FROZEN.md 에 있다.
import glob as _glob

TEX = 'paper_acs/main.tex'
# 투고 패키지에 들어가는 모든 글. 본문과 같은 산출물에 대조한다.
# 2026-11-02. paper_acs/supporting.tex 를 넣는다. 내용이
# paper/supplementary.tex 와 거의 같은 평행 사본이라, 한쪽만 고치면 조용히
# 갈라진다. 둘 다 같은 산출물에 대조하면 어느 쪽에 구판 값이 남아도 걸린다.
# 생성된 표 파일도 대조 대상에 넣는다. 넣지 않으면 그 표의 셀을 가리키는
# 앵커들이 '문맥 미발견' 으로 조용히 빠진다.
# 2026-11-07 (11차 6단계). 종전에는 여덟 개를 손으로 적어 두었고, 뒤에
# 생긴 tab_registers / tab_designations / tab_dyepanel 셋은 이 목록과 아래
# 신선도 검사 양쪽에서 빠져 있었다. 본문에서 레지스터 수치를 SI 표로
# 내리는 순간 그 수치가 무검증 상태가 되는 구조였으므로 글롭으로 바꾼다.
_TAB_TEX = sorted(_glob.glob('paper_acs/tab_*.tex'))
# 2026-11-07 (11차). README 도 대조 대상에 넣는다. README 의 재현 약속
# 단락이 '최소 p 0.1512 / BH 0.159' 같은 구판 수치를 들고 있었고(실측은
# 0.1917 / 0.1632), 그것을 보고 돌린 사람은 코드가 틀렸다고 판단하게 된다.
# 원고만 검사하고 재현 지시서를 검사하지 않으면 같은 종류의 거짓이 남는다.
EXTRA_TEX = ('paper_acs/supporting.tex', 'paper_acs/cover_letter.tex',
             'README.md', *_TAB_TEX)

# 사실표가 읽는 산출물 전부 (2026-10-17 신설).
#
# 왜 필요한가. 저장소를 새로 클론해 검증기를 돌리자 outputs/mde.json 하나가
# 없다는 FileNotFoundError 로 즉시 죽었고, 그 뒤에 또 무엇이 없는지는 알 수
# 없었다. 실제로는 산출물 4건과 도판 1건이 빠져 있었다. 한 건씩 터지는 대신
# 빠진 것을 한 번에 전부 보고한다. 이 목록은 '저장소가 담아야 하는 산출물'의
# 계약서 구실도 한다.
REQUIRED_OUTPUTS = [
    'outputs/aegl_parse_validation.json',
    'outputs/blockbalance.json',
    'outputs/collapse_sensitivity.json',
    'outputs/length_audit.json',
    'outputs/reversal_audit.json',
    'outputs/dye_panel_audit.json',
    'outputs/descriptor_gate_scope.json',
    'outputs/dim_percentiles.json',
    'outputs/distance_table_normalized.csv',
    'outputs/experiment_summary.json',
    'outputs/figure_annotations.json',
    'outputs/h1_candidates.csv',
    'outputs/h1_summary.json',
    'outputs/h1_screen_sensitivity.json',
    'outputs/pca_facts.json',
    'outputs/match_by_panel.json',
    'outputs/mde.json',
    'outputs/minority_regression.json',
    'outputs/multiplicity_c.json',
    'outputs/permutation_sidedness.json',
    'outputs/pesticide_panel.json',
    'outputs/cite_facts.json',
    'outputs/pesticide_family.json',
    'outputs/qsar_heldout.csv',
    'outputs/qsar_tox_summary.json',
    'outputs/rollback_c.json',
    'outputs/stat_supplement.json',
    'outputs/stat_supplement2.json',
    'outputs/toxicity_aegl.csv',
    'outputs/toxicity_predicted.csv',
]

# 도판 스크립트가 추가로 읽는 것. 검증기 자체는 쓰지 않으나 이것이 없으면
# build_paper.sh --figs 가 깨지고, 그러면 원고에 구판 그림이 남는다.
REQUIRED_FOR_FIGURES = [
    'outputs/_electronic_with_gap.csv',
    # 2026-11-03. 휘발성 축 보조 분석 (volatility_space.py).
    'outputs/volatility_candidates.json',
    'outputs/volatility_space.json',
    # 2026-11-03 (작업 A).
    'outputs/shape_reduction.json',
    'outputs/premise_null.json',
    # 2026-11-04 (E2).
    'outputs/lewisite_isomer.json',
    'outputs/_pca_scores.npy',
    'outputs/_null_A_rdkit.npy',
    'outputs/_null_B_electronic.npy',
    'outputs/_null_C_merged.npy',
    'outputs/electronic_distance_table.csv',
    'outputs/electronic_space_summary.json',
    'outputs/h1_gate_tiers.json',
    'outputs/sanity_checks.json',
    'outputs/fig_label_counts.json',
    'outputs/nearest_agent_bootstrap.json',
    'outputs/null_match_expectation.json',
]

# 2026-11-04 (J-C). ACS본의 도판 집합. 종전 목록은 JHM 평행본의 것이어서
# Fig_spaces·Fig_toxicity 를 요구했는데, ACS본은 그 둘을 Fig_talk_response·
# Fig_talk_tox 로 바꿨다. build_paper_acs.sh 가 복사하는 목록과 맞춘다.
FIGURES = ['Fig_main.png', 'Fig_plocal.png', 'Fig_h1.png', 'Fig_mde.png',
           'Fig_talk_response.png', 'Fig_talk_tox.png', 'toc_acs.png']

# 도판 신선도 검사의 mtime 허용오차(초). 2026-10-22. 종전 1초 기준은 갓 클론한
# 저장소에서 거짓으로 걸린다. git clone 은 모든 파일의 mtime 을 체크아웃 시각으로
# 새로 쓰는데, 13 MB 짜리 outputs/ 를 풀다 보면 초 경계를 넘어가기 때문이다.
# 논문이 인용하는 v1.0 을 리뷰어가 클론해서 바로 이 검증기를 돌린다는 점에서
# 거짓 실패는 치명적이다. 실제 구판은 날 단위로 벌어지므로(그래픽 초록 적발
# 당시 11일) 120초면 둘을 가른다.
MTIME_TOL = 120


def _dropped_periods():
    """되쓰기가 문장 끝 마침표를 먹은 자리를 찾는다 (2026-11-03 신설).

    왜 필요한가. --fix 가 수치 157건을 되쓸 때, 캡처가 ([\\d.]+) 꼴이어서
    문장 끝 마침표까지 먹은 자리 네 곳에서 마침표가 사라졌다. tonum 이
    그 마침표를 떼어내고 비교하므로 검사 쪽은 통과했고, 결함은 조판된
    PDF 를 읽어야 보이는 상태로 남았다. --fix 쪽은 고쳤지만, 같은 꼴이
    손 편집으로도 생기므로 가드를 둔다.

    본문(주석 제거)에서 숫자 바로 뒤에 공백과 대문자로 시작하는 단어,
    빈 줄, 또는 환경·절 명령이 오는 자리를 찾는다. 날짜 뒤 줄바꿈 같은
    정상 꼴은 면제한다.
    """
    import os as _o
    import re as _r
    pats = [(_r.compile(r'\d\n[ \t]*\n'), '숫자 뒤 빈 줄'),
            (_r.compile(r'\d\s+\\(?:begin|section|subsection|bibliography)'),
             '숫자 뒤 환경/절 명령'),
            (_r.compile(r'\d\s+[A-Z][a-z]'), '숫자 뒤 대문자 단어')]
    # 날짜 블록과 저자 블록은 숫자로 끝나고 다음 줄이 대문자로 시작한다.
    exempt = _r.compile(r'(?:January|February|March|April|May|June|July|'
                        r'August|September|October|November|December)\s+\d{4}|'
                        r'vspace|\\author|\\affiliation|\\email|'
                        r'Republic of Korea')
    out = []
    # 2026-11-07. 이 가드는 산문용이다. 생성된 표(tab_*.tex)의 셀과 각주는
    # 숫자 뒤에 대문자가 정상적으로 오므로(예: '1 Anthracene') 제외한다.
    # 표의 수치 검증은 앵커와 신선도 검사가 따로 맡는다.
    for f in (TEX, 'paper_acs/supporting.tex', 'paper_acs/cover_letter.tex'):
        if not _o.path.exists(f):
            continue
        raw = open(f).read()
        lines = []
        for ln in raw.split('\n'):
            if ln.lstrip().startswith('%'):
                lines.append('')
            else:
                lines.append(_r.sub(r'(?<!\\)%.*$', '', ln))
        body = '\n'.join(lines)
        for pat, why in pats:
            for m in pat.finditer(body):
                win = body[max(0, m.start() - 120):m.end() + 60]
                if exempt.search(win):
                    continue
                ctx = ' '.join(body[max(0, m.start() - 60):m.end() + 30].split())
                out.append(f'  마침표 누락 의심 [{f}] {why}: {ctx}')
    return sorted(set(out))


def _ordinal_typos():
    """소수 뒤 서수 접미사 오류를 찾는다 (2026-11-03 신설).

    왜. 백분위를 되쓰면서 '24.3rd' 의 접미사를 그대로 두고 숫자만 바꿔
    '22.4rd' 가 되었다. 소수로 끝나는 수의 서수는 영어에서 항상 th 다.
    GATE 8 8차 리뷰어가 네 자리에서 이걸 찾아냈다.
    """
    import os as _o
    import re as _r
    bad = []
    for f in (TEX,) + EXTRA_TEX:
        if not _o.path.exists(f):
            continue
        for m in _r.finditer(r'\d\.\d(?:rd|st|nd)\b', open(f).read()):
            bad.append(f'  서수 접미사 오류 [{f}]: {m.group(0)} -> '
                       f'{m.group(0)[:-2]}th')
    return sorted(set(bad))


def _handmade_tables():
    r"""산출물에서 생성되지 않는 수치 표를 찾는다 (2026-11-03 신설).

    왜 필요한가. SI 의 tab:ablation 은 10열 24수치를 손으로 유지하고
    있었고, 어떤 스크립트도 만들지 않았다. 그래서 T2-1(상관필터 교체)과
    T2-8(형태이성질체 재생성)을 전부 통과해 차원 열이 40/42 구판으로
    살아남았고, 검증기 458건도 통과시켰다. G-D 가드는 인용된 예치 '파일'
    만 보고 표는 보지 않았다.

    여기서는 방향을 뒤집는다. 수치가 여러 줄 들어가는 table 환경은
    \input 으로 들어와야 하며, 본문 안에 직접 적혀 있으면 생성기가 없는
    것으로 보고 보고한다. 판정 기준은 '표 본문에 소수가 12개 이상' 이다.
    작은 대조표나 문자 표는 걸리지 않는다. 면제가 필요한 표는 캡션에
    'hand-tabulated' 를 적어 사유를 남긴다.
    """
    import os as _o
    import re as _r
    out = []
    for f in (TEX,) + EXTRA_TEX:
        if not _o.path.exists(f):
            continue
        s_ = open(f).read()
        # 2026-11-03. 생성된 표 파일 자체가 이 가드에 걸렸다. 생성기가
        # 찍는 머리글이 있으면 '산출물에서 만들어진 표'이므로 건너뛴다.
        if 'make_si_tables.py' in s_.split('\n')[0]:
            continue
        for m in _r.finditer(r'\\begin\{table\}(.*?)\\end\{table\}', s_, _r.S):
            blk = m.group(1)
            lab = _r.search(r'\\label\{([^}]*)\}', blk)
            body = _r.sub(r'\\caption\{.*?\n\n', '', blk, flags=_r.S)
            nums = _r.findall(r'(?<![\w.])\d+\.\d+', body)
            if len(nums) >= 12 and 'hand-tabulated' not in blk:
                out.append(f'  생성기 없는 수치 표 [{f}] '
                           f'{lab.group(1) if lab else "(라벨 없음)"}: '
                           f'소수 {len(nums)}개가 본문에 직접 적혀 있다. '
                           'make_si_tables.py 로 생성하거나 캡션에 '
                           "'hand-tabulated' 와 사유를 적을 것")
    return sorted(set(out))


def _generated_table_sync():
    r"""생성된 표의 원고 사본이 outputs/ 와 어긋나지 않는지 본다.

    왜 필요한가 (2026-11-03 신설). 손표를 생성으로 돌리는 것만으로는
    전파 실패가 사라지지 않는다. 생성기를 다시 돌려도 paper_acs/ 와
    paper/ 의 사본이 구판이면 원고에 박히는 수치는 그대로다. 도판에서
    같은 일이 이미 두 번 있었고(GATE 8 3차, 7차), 표에서도 같은 구멍이
    난다. 빌드 스크립트가 복사하지만, 빌드를 거치지 않고 검증기만 돌린
    경우에 조용히 통과하는 것을 막는다. 사본이 다르면 어느 파일이
    어긋났는지 적어 보고한다.
    """
    import os as _o
    bad = []
    for name in [_o.path.basename(t) for t in _TAB_TEX]:
        src = f'outputs/{name}'
        if not _o.path.exists(src):
            bad.append(f'  생성 표 없음 [{src}]: make_si_tables.py 를 돌릴 것')
            continue
        ref = open(src).read()
        # 2026-11-04 (J-C). paper/ 는 archive/jhm_v1/ 로 동결했다.
        for d in ('paper_acs',):
            dst = f'{d}/{name}'
            if not _o.path.exists(dst):
                bad.append(f'  생성 표 사본 없음 [{dst}]: '
                           'sh build_paper_acs.sh 를 돌릴 것')
            elif open(dst).read() != ref:
                bad.append(f'  생성 표 사본이 구판 [{dst}]: {src} 와 다르다. '
                           'sh build_paper_acs.sh 로 동기화할 것')
    return sorted(set(bad))


def _filter_dependent_staleness():
    r"""상관필터에 의존하는 예치물이 정본 실행보다 오래됐는지 본다.

    왜 필요한가 (2026-11-04 신설, 작업 D). 11-02 08:17 의 구판 실행이 만든
    예치물 일곱 개가 재생성되지 않고 저장소에 남아 있었다. 전부 읽는 코드도
    생성기도 본문 인용도 없는 고아였고, 그래서 기존 고아 검사(본문이 인용한
    예치물에 생성기가 있는가)를 빠져나갔다. 그 안에는 전자공간 좌표를 4개로
    (정본 3), 28종 구조공간을 38개로 (정본 35) 적는 값이 들어 있었고,
    sync_submission.sh 가 outputs/*.json 을 전량 복사하므로 리뷰어가 받는
    패키지에 그대로 들어갔다. 네 명이 같은 자리를 지적했다.

    더 나쁜 것은 distance_table_normalized*.csv 네 벌이 검증기 자신의 거리
    화이트리스트를 넓히고 있었다는 점이다. 구판 거리가 '허용된 값' 이 되면
    본문의 구판 인용이 통과한다.

    여기서는 상관필터와 공간 정의에 의존하는 이름만 본다. 웹에서 받은
    출처(pubchem, volatility_candidates), DFT 실행, 캐시는 정본 실행보다
    오래되어도 정상이므로 대상이 아니다.
    """
    import os as _o
    import glob as _g
    FAMILIES = ('distance_table', 'experiment_summary', 'match_by_panel',
                'panel35_', 'subspace_null', 'permutation_',
                # 2026-11-04 (작업 D, 2차). deltaE_extended.json 이 구판
                # Mantel p 로 Fisher 를 합성해 본문이 0.204 를 싣고 있었다.
                # 정본 p 로 다시 합성하면 0.2373 이다. R3·R4 중복 지적.
                'deltaE_', 'rollback_', 'multiplicity', 'mantel_')
    ref = 'outputs/experiment_summary.json'
    if not _o.path.exists(ref):
        return []
    t0 = _o.path.getmtime(ref)
    bad = []
    for f in sorted(_g.glob('outputs/*.json') + _g.glob('outputs/*.csv')):
        b = _o.path.basename(f)
        if not b.startswith(FAMILIES):
            continue
        if _o.path.getmtime(f) < t0 - MTIME_TOL:
            bad.append(f'  구판 예치물 [{f}]: 정본 실행('
                       f'{ref})보다 오래됐다. 파이프라인에 넣어 재생성하거나 '
                       f'삭제할 것 — 패키지에 들어가면 리뷰어가 본문과 다른 '
                       f'수치를 받는다')
    return bad


def _si_section_labels():
    r"""SI 의 모든 절에 라벨이 있고, 그 라벨이 어딘가에서 인용되는지 본다.

    왜 필요한가 (2026-11-04 신설, G3). SI 26개 절 중 세 개에 라벨이 없어
    아예 가리킬 수 없었고, 세 개는 라벨은 있는데 본문·SI 어디서도 인용되지
    않았다. 그중 하나(res:mde)는 바로 위에 같은 제목의 빈 subsection 이
    붙어 있어 목차에 두 번 찍혔다. 조판은 전부 통과한다.
    """
    import os as _o
    import re as _re
    f_si, f_mn = 'paper_acs/supporting.tex', 'paper_acs/main.tex'
    if not (_o.path.exists(f_si) and _o.path.exists(f_mn)):
        return []
    strip = lambda t: '\n'.join(_re.sub(r'(?<!\\)%.*$', '', l)
                                for l in t.split('\n'))
    si = strip(open(f_si, encoding='utf8').read())
    mn = strip(open(f_mn, encoding='utf8').read())
    secs = [(m.start(), ' '.join(m.group(1).split()))
            for m in _re.finditer(r'\\section\*?\{(.+?)\}', si, _re.S)]
    labs = {m.start(): m.group(1)
            for m in _re.finditer(r'\\label\{([^}]+)\}', si)}
    bad = []
    titles = [t for _, t in secs]
    for i, (pos, title) in enumerate(secs):
        nxt = secs[i + 1][0] if i + 1 < len(secs) else len(si)
        lab = next((l for p, l in sorted(labs.items()) if pos < p < nxt), None)
        if lab is None:
            bad.append(f'  SI 절에 라벨이 없다: "{title[:52]}" — 가리킬 수 '
                       f'없으므로 항해 불가')
            continue
        cited = (_re.findall(r'ref\{' + _re.escape(lab) + r'\}', si)
                 or _re.findall(r'ref\{S-' + _re.escape(lab) + r'\}', mn))
        if not cited:
            bad.append(f'  SI 절이 어디서도 인용되지 않는다: {lab} '
                       f'("{title[:40]}")')
    for t in sorted({t for t in titles if titles.count(t) > 1}):
        bad.append(f'  SI 절 제목 중복: "{t[:52]}"')
    return bad


def _figure_panel_refs():
    r"""인용된 도판 패널 문자가 그 도판 캡션에 실제로 있는지 본다.

    왜 필요한가 (2026-11-04 신설, 작업 D). ACS본은 JHM본의 4패널 독성 도판을
    2패널 Fig_talk_tox 로 바꿨는데 SI 의 인용은 'Fig.~\ref{M-fig:tox}c and d'
    로 남아 있었다. 없는 패널을 가리키므로 조판은 통과하고 독자는 찾지
    못한다. 수치 검사로는 잡히지 않는 종류다.
    """
    import os as _o
    import re as _re
    caps = {}
    for f in ('paper_acs/main.tex',):
        if not _o.path.exists(f):
            continue
        txt = open(f, encoding='utf8').read()
        for blk in _re.findall(r'\\begin\{figure\*?\}(.*?)\\end\{figure\*?\}',
                               txt, _re.S):
            lab = _re.search(r'\\label\{(fig:[^}]+)\}', blk)
            if not lab:
                continue
            caps[lab.group(1)] = set(_re.findall(r'\((\w)\)', blk))
    bad = []
    for f in ('paper_acs/main.tex', 'paper_acs/supporting.tex'):
        if not _o.path.exists(f):
            continue
        txt = open(f, encoding='utf8').read()
        for m in _re.finditer(r'ref\{(?:M-)?(fig:[^}]+)\}([a-z])\b', txt):
            lab, panel = m.group(1), m.group(2)
            if lab in caps and panel not in caps[lab]:
                bad.append(f'  없는 패널 인용 [{f}]: {lab} 에 ({panel}) 가 '
                           f'없다. 캡션의 패널은 '
                           f'{sorted(caps[lab]) or "없음"}')
    return sorted(set(bad))


def _title_sync():
    r"""제목이 본문·SI·커버레터에서 같은지 본다 (2026-11-03 신설, 작업 B).

    왜 필요한가. 제목의 'Sixteen' 을 빼면서 main.tex 만 고쳤고 SI 와
    커버레터는 구판 제목을 들고 있었다. 금지항목이 그 한 문자열을 잡아
    주었을 뿐이고, 다음 제목 교체에서는 금지항목이 없으므로 조용히
    통과한다. 세 파일이 같은지를 직접 본다.
    """
    import re as _re
    import os as _o
    want, bad = None, []
    for path, pat in [
            ('paper_acs/main.tex', r'\\title\{(.+?)\}'),
            ('paper_acs/supporting.tex', r'\\title\{(.+?)\}'),
            ('paper_acs/cover_letter.tex',
             r'We submit \\textit\{(.+?)\} for consideration'),
            # 2026-11-07 (12차 R5-B3). README 머리말도 제목을 싣는다.
            # 세 tex 만 보던 탓에 README 가 구판 제목·저자 2인·JHM 투고를
            # 들고 공개돼 있었고, 예치물 링크를 따라간 심사자가 가장 먼저
            # 보는 자리였다.
            ('README.md', r'^> \*\*(.+?)\*\*$')]:
        if not _o.path.exists(path):
            continue
        m = _re.search(pat, open(path, encoding='utf8').read(),
                       _re.S | _re.M)
        if not m:
            bad.append(f'  제목을 찾지 못함 [{path}]')
            continue
        # README 는 인용블록이라 줄이 이어지면 '> ' 가 끼어든다.
        got = ' '.join(t for t in m.group(1).split() if t != '>')
        if want is None:
            want = got
        elif got != want:
            bad.append(f'  제목이 본문과 다르다 [{path}]:\n'
                       f'    본문 {want}\n    이쪽 {got}')
    return bad


def _label_scope_guards():
    r"""수치는 맞고 그 수치를 설명하는 '말'이 틀린 경우 (2026-11-03 신설).

    왜 필요한가. 리뷰어 1 이 결함 10건을 들고 왔고 전부 사실이었는데,
    검증기 505건은 모두 통과했다. 이유는 하나다. -0.189 도, 0.3539 도,
    59.3 도 산출물에 실제로 있는 값이다. 틀린 것은 그 값을 설명하는
    단어였다: 'four' G agents(실제 3종), 'eighteen' tests(실제 72검정),
    'tier one'(출처가 Extrapolated 라고 적은 값), 'exhaustive'(실제 표집).
    종전 가드는 '표 셀이 구판인가'를 보았고 이것은 다른 계열이다.

    기계적으로 잡히는 세 가지만 여기서 본다. 분모 오류(3종을 4종이라
    적는 류)는 문장마다 앵커를 붙이는 수밖에 없고, 그렇게 했다.
    """
    import json as _j
    import os as _o
    import re as _r
    bad = []

    # (A) 전수열거 주장 대 표집 산출물. #6 이 이 경로였다.
    CLAIM = _r.compile(r'exhaustiv\w+', _r.I)
    sampled = []
    for fn in sorted(_o.listdir('outputs')):
        if not fn.endswith('.json'):
            continue
        try:
            txt = open(f'outputs/{fn}').read()
        except OSError:
            continue
        if '"enumerated": false' in txt.replace(' ', '').replace(
                '"enumerated":false', '"enumerated": false'):
            sampled.append(fn)
        elif _r.search(r'"(n_draws|n_rep|n_replicates)":\s*\d+', txt):
            sampled.append(fn)
    for f in (TEX,) + EXTRA_TEX:
        if not _o.path.exists(f):
            continue
        raw = open(f).read()
        body = '\n'.join('' if ln.lstrip().startswith('%')
                          else _r.sub(r'(?<!\\\\)%.*$', '', ln)
                          for ln in raw.split('\n'))
        for m in CLAIM.finditer(body):
            win = ' '.join(body[max(0, m.start() - 320):m.end() + 320].split())
            # 범위를 한정한 문장은 통과시킨다. 'sampled' 를 함께 적고
            # 표집 횟수를 밝힌 문장이 그것이다.
            if _r.search(r'are sampled|sampled controls|sampled at|draws',
                         win, _r.I):
                continue
            # 'exhaustive representations of all compounds' 처럼 순열검정과
            # 무관한 뜻으로 쓰인 자리는 대상이 아니다.
            if not _r.search(r'permutation|enumerat|relabell', win, _r.I):
                continue
            bad.append(f'  전수열거 주장에 범위가 없다 [{f}]: '
                       f'"...{win[260:400]}..." — 표집 산출물 '
                       f'{len(sampled)}건 존재 ({", ".join(sampled[:3])} 등). '
                       '전수인 검정과 표집인 검정을 문장에서 구분할 것')

    # (B0) 내부 문서가 공개 예치물에 들어가면 멈춘다 (2026-11-03 신설).
    #      리뷰어 4 가 outputs/HO-7.json 을 읽었다. 그 파일에는 타깃 저널이
    #      JHM 으로, 자기인용 비율 축소 여부가 'Ku 판단' 으로, 교신저자
    #      이메일 문제와 A계열 동정 미검증이 적혀 있었다. 원고가 저장소를
    #      예치물로 인용하므로 심사자가 저자의 내부 심의를 읽게 된다.
    INTERNAL = ('handoff_id', 'docx_report', '확인 필요', '자기인용',
                'elsarticle', 'Hazardous Materials 투고')
    for fn in sorted(_o.listdir('outputs')):
        if not fn.endswith('.json'):
            continue
        txt = open(f'outputs/{fn}').read()
        hit = [w for w in INTERNAL if w in txt]
        # 인용 메타데이터는 저널명을 정당하게 담는다. 핸드오프 카드만 잡는다.
        if 'handoff_id' in hit or 'docx_report' in hit or (
                '자기인용' in hit and 'journal' not in txt[:200]):
            bad.append(f'  내부 문서가 예치물에 있다 [outputs/{fn}]: '
                       f'{hit} — 공개 저장소에서 제거할 것')

    # (B0b) 내부 심의 문서와 폐기본이 추적되고 있는지 (2026-11-07 신설,
    #       12차 R5-A5/A6).
    #
    # 왜 필요한가. 위 B0 는 outputs/*.json 만 훑는다. 그래서 sources/ 의
    # GATE8 채점 메모 7건과 outputs/_superseded/ 18건이 공개 저장소에
    # 그대로 들어가 있었고 어느 검사도 묻지 않았다. 가드의 범위가 위험
    # 보다 좁으면 가드가 아니다. 경로 패턴으로 직접 센다.
    import subprocess as _sp2
    _INTERNAL_PAT = ('sources/GATE8_', 'sources/T1_', 'sources/T2_',
                     'sources/NOTE_', 'outputs/_superseded/')
    try:
        _tracked = _sp2.run(['git', '-c', 'core.quotepath=false',
                             'ls-files'], capture_output=True, text=True,
                            timeout=60)
        _files = _tracked.stdout.splitlines() if _tracked.returncode == 0 \
            else None
    except Exception:
        _files = None
    if _files is None:
        print('[verify] 검사 B0b: git 을 쓸 수 없어 건너뛴다.')
    else:
        _leak = [f for f in _files
                 if any(f.startswith(p_) for p_ in _INTERNAL_PAT)]
        if _leak:
            bad.append(f'  내부 문서·폐기본이 예치물에 추적되고 있다 '
                       f'(검사 B0b): {len(_leak)}건, 예: '
                       f'{", ".join(_leak[:3])} — .gitignore 에 넣고 '
                       'git rm --cached 할 것')

    # (B) tier 주장 대 출처 문자열. #8 이 이 경로였다.
    vp = 'outputs/volatility_space.json'
    if _o.path.exists(vp):
        d = _j.load(open(vp))['vapour_pressure_mmHg']
        BAD_WORD = _r.compile(r'extrapolat|estimat|predict|calculat', _r.I)
        for n, e in d.items():
            if e['tier'] in ('T1', 'T2') and BAD_WORD.search(e['source']):
                bad.append(f"  측정값 tier 에 추정값 [{n}]: tier {e['tier']} "
                           f"인데 출처 문자열이 '{e['source'][-30:]}' 다. "
                           'TX 로 내리거나 다른 온도의 측정값을 쓸 것')

    # (C) 생성 표 캡션이 말하는 가족 크기 대 실제로 읽은 산출물. #4 가 이 경로였다.
    #     tab_ablation 의 마지막 열은 '해당 패널 구성의 18검정' 이어야 한다.
    tf = 'outputs/tab_ablation.tex'
    mc = 'outputs/multiplicity_c.json'
    if _o.path.exists(tf) and _o.path.exists(mc):
        want = {round(r[4], 4) for r in _j.load(open(mc))}
        txt = open(tf).read()
        if 'eighteen tests' in txt:
            got = {round(float(x), 4) for x in
                   _r.findall(r'&\s*([\d.]+)\s*\\\\', txt)}
            if not (want & got):
                bad.append('  생성 표의 캡션과 데이터가 다른 가족 '
                           f'[{tf}]: 캡션은 18검정 BH 라고 적는데 '
                           f'{sorted(want)} 가 표에 없다. 72검정 파일을 '
                           '붙였는지 확인할 것')
    return sorted(set(bad))


def preflight():
    """산출물·도판 존재를 먼저 전수 확인한다. 없으면 전부 나열하고 멈춘다."""
    import os
    missing = [p for p in REQUIRED_OUTPUTS if not os.path.exists(p)]
    figless = [p for p in REQUIRED_FOR_FIGURES if not os.path.exists(p)]
    nofig = [f for f in FIGURES if not os.path.exists(f'paper_acs/{f}')]
    # 그래픽 초록은 build_paper.sh 바깥에 있어서 분석을 다시 돌려도 갱신되지
    # 않았다. 2026-10-22 투고 직전에 p=0.297, Spearman r=-0.022, n=214 를
    # 그대로 들고 있는 10-11 판이 발견됐다. 산출물보다 오래되면 멈춘다.
    if os.path.exists('outputs/graphical_abstract.png'):
        _ga = os.path.getmtime('outputs/graphical_abstract.png')
        _src = max(os.path.getmtime(p) for p in
                   ('outputs/sanity_checks.json', 'outputs/h1_summary.json',
                    'outputs/h1_candidates.csv', 'outputs/stat_supplement.json')
                   if os.path.exists(p))
        if _ga < _src - MTIME_TOL:
            # 2026-11-06. 메시지를 고친다. 이 항목은 outputs/ 의 파일인데
            # 호출부가 'paper_acs/' 를 붙여 출력하고 있었다. paper_acs 에는
            # 그런 파일이 없으므로(투고 시스템에 따로 올리는 첨부다) 읽는
            # 사람이 없는 경로를 찾게 된다. 9차에서 세 번 헛돌았다.
            nofig.append('[outputs/] graphical_abstract.png 가 산출물보다 '
                         '오래됐다 — python3 make_ga.py 로 재생성할 것')
    # 2026-10-20 신설. figure_annotations.json 은 누적만 하므로, 폐기된 도판의
    # 주석이 남아 본문 검사를 대신 통과시킬 수 있다. 등재부에 현행 도판이 아닌
    # 항목이 하나라도 있으면 멈춘다.
    stale = []
    if os.path.exists('outputs/figure_annotations.json'):
        import fig_annot as _fa
        # 2026-11-07 (11차 냉클론 시험). graphical_abstract 는 구판
        # make_ga.py (JHM 판) 의 첨부였고 ACS 빌드는 그 스크립트를 부르지
        # 않는다. 그래서 빈 저장소를 클론하면 그 등재가 아예 없고, 그
        # 등재에 걸려 있던 앵커(ga_q28)가 통째로 꺼졌다. make_ga.py 를
        # archive/jhm_v1/ 로 옮기고 등재를 지웠다. ACS 의 TOC 도판은
        # make_toc_acs.py 가 자기 글자를 'toc_acs' 로 등재한다.
        # 2026-11-04 (J-C). 종전에는 FIGURES(현 원고가 싣는 도판)와
        # 비교했다. JHM 평행본을 동결하고 FIGURES 를 ACS 집합으로 줄이자
        # Fig_spaces·Fig_toxicity 가 폐기로 걸렸는데, 그 둘은 생성기가
        # 여전히 만들고 있다(make_fig_main.py, make_fig_toxicity.py).
        # 이 검사의 목적은 '만드는 곳이 없어진 주석' 을 잡는 것이므로
        # outputs/ 에 실제로 생성된 도판 전부와 비교한다.
        live = ({f[:-4] for f in FIGURES} | {'toc_acs'}
                | {os.path.basename(x)[:-4]
                   for x in glob.glob('outputs/Fig_*.png')})
        stale = sorted(set(_fa.load()) - live)
    # 2026-11-02 신설. 원고가 '반복측정' 을 주장했는데 원자료에는 조합당
    # 한 행뿐이었다. 주장을 자료에 걸어 둔다. 어느 피벗에든 (dye, agent,
    # conc) 가 두 번 나오면 그때는 반복측정이 실제로 있는 것이므로, 위의
    # SUPERSEDED 금지항목을 함께 풀어야 한다는 뜻이다.
    # 2026-11-02 신설 (GATE 8 7차). 인용된 예치물의 생성자 검사.
    #
    # 왜. outputs/multiplicity_across_panels.json 이 저장소에 예치되어 있고
    # SI 가 이름을 대서 인용하는데, 그 파일을 만드는 코드가 없었다. 검증기만
    # 읽고 있었다. 생성 코드 없는 예치물의 수치는 출처 불명이고, 실제로 그
    # 값이 현행 산출물과 어긋나 있었다 (0.156 대 0.169).
    #
    # 범위. outputs/*.json 전부를 정적으로 훑어 '생성 vs 읽기' 를 가르려 했으나
    # 파일명을 변수로 조립하는 스크립트가 많아 위양성이 쏟아졌다. 그래서
    # 원고가 \texttt{} 로 이름을 대서 인용한 것만 본다. 독자가 찾아갈 수
    # 있다고 원고가 약속한 파일이므로, 그 약속의 범위가 맞다.
    #
    # 새 파일을 인용하면 여기에 생성자를 등록해야 통과한다. 등록을 강제하는
    # 것이 이 검사의 요점이다.
    GENERATED_BY = {
        'conformer_sensitivity_ensemble.json': 'conformer_sensitivity.py',
        'conformer_sensitivity_single.json': 'conformer_sensitivity.py',
        'multiplicity_across_panels.json': 'multiplicity_panels.py',
        # 2026-11-03. 휘발성 보조 분석.
        'volatility_candidates.json': 'volatility_fetch.py',
        'volatility_space.json': 'volatility_space.py',
        'shape_reduction.json': 'shape_reduction.py',
        'collapse_sensitivity.json': 'collapse_sensitivity.py',
        'length_audit.json': 'wordcount.py',
        'reversal_audit.json': 'reversal_audit.py',
        'dye_panel_audit.json': 'dye_panel_audit.py',
        'premise_null.json': 'premise_null.py',
        'stat_supplement.json': 'stat_supplement.py',
        'lewisite_isomer.json': 'lewisite_isomer.py',
        'subspace_null.json': 'subspace_null.py',
        'dft_dipole.csv': 'dft_omega.py',
        'dft_electronic_space.json': 'dft_electronic_space.py',
        'dft_omega.csv': 'dft_omega.py',
        # 2026-11-02 (T2). 새 예치물 넷.
        'order_invariance.json': 'order_invariance.py',
        'match_rule_audit.json': 'match_rule_audit.py',
        'deltaE_noise.json': 'deltaE_noise.py',
        'purity_sensitivity.json': 'purity_sensitivity.py',
        'channel_classification.json': 'channel_classification.py',
        # 2026-11-03 (T1-14). 희소 지시자 인코딩 민감도.
        'sparse_indicators.json': 'sparse_indicators.py',
    }
    orphan = []
    _cited = set()
    for _doc in ('paper_acs/main.tex', 'paper_acs/supporting.tex'):
        if not os.path.exists(_doc):
            continue
        _t = '\n'.join(l for l in open(_doc).read().split('\n')
                        if not l.lstrip().startswith('%'))
        for _m in re.finditer(r'\\texttt\{([^}]*\.(?:json|csv))\}', _t):
            _cited.add(os.path.basename(_m.group(1).replace('\\_', '_')))
    for _name in sorted(_cited):
        _maker = GENERATED_BY.get(_name)
        if _maker is None:
            orphan.append(f'{_name} 를 원고가 인용하는데 GENERATED_BY 에 '
                          f'생성 스크립트가 등록되지 않았다')
        elif not os.path.exists(_maker):
            orphan.append(f'{_name} 의 생성 스크립트 {_maker} 가 없다')
        elif not glob.glob(f'outputs/{_name}'):
            orphan.append(f'{_name} 를 원고가 인용하는데 outputs/ 에 없다')
    if orphan:
        print('[verify] 사전 점검 실패 — 인용된 예치물.')
        for _m in orphan:
            print(f'   {_m}')
        return False

    # 2026-11-02 신설 (GATE 8 7차). 작용제 수 서술 검사.
    #
    # 왜. 검증기가 417건 PASS 인 상태에서 제목이 작용제를 17종이라고 적고
    # 있었다. 실제는 16종이고 17은 모사체 DMMP 를 포함한 화합물 수다.
    # (B) 검사는 등록된 어법의 '숫자'만 보므로 영어 수사는 보이지 않는다.
    #
    # 범위를 좁힌 이유. 처음에는 모든 수사를 피벗 실측과 대조했는데, 원고에는
    # 정당한 다른 수가 여럿 있다 ('eleven agents' = 서술자 패널의 작용제 수,
    # 'seven-agent run' = 2차 런). 그래서 실제 실패 양식 하나만 잡는다:
    # 측정 화합물 수(17)를 작용제 수로 적는 것. 'non-agents' 는 C(17,3) 의
    # 비작용제 17종이므로 면제한다.
    cnt = []
    agents = set()
    for _p in sorted(glob.glob('deltaE/d?.uv_deltaE_pivot.csv')):
        try:
            agents |= set(pd.read_csv(_p)['Agent'].astype(str))
        except Exception as e:
            cnt.append(f'{_p}: 읽기 실패 {e}')
    SIMULANTS = {'DMMP', 'DPGME', 'DCM'}
    n_cpd, n_agent = len(agents), len(agents - SIMULANTS)
    if agents and n_cpd == 17 and n_agent == 16:
        for _doc in ('paper_acs/main.tex', 'paper_acs/supporting.tex',
                     'paper_acs/cover_letter.tex'):
            if not os.path.exists(_doc):
                continue
            body = '\n'.join(l for l in open(_doc).read().split('\n')
                              if not l.lstrip().startswith('%'))
            # 제목은 'Seventeen Authentic Agents' 처럼 대문자다. re.I 필수
            # (2026-11-02 주입시험에서 이것 때문에 처음에 안 잡혔다).
            for m in re.finditer(r'seventeen\s+(?:\w+\s+){0,2}?agents?\b',
                                 body, re.I):
                if 'non-agent' in m.group(0):
                    continue
                cnt.append(f'{_doc}: "{" ".join(m.group(0).split())}" '
                           f'— 측정 화합물은 17종이지만 작용제는 {n_agent}종이다 '
                           f'(DMMP 는 모사체). sixteen 으로 적을 것')
    if cnt:
        print('[verify] 사전 점검 실패 — 작용제 수 서술.')
        for _m in cnt:
            print(f'   {_m}')
        return False

    # 2026-11-02 신설. ACS 판본 도판 점검. 검증기는 paper_acs/main.tex 를
    # 글로는 읽지만 그 도판이 실재하는지, 최신인지는 보지 않았다. TOC 도판은
    # outputs/deltaE_premise.json 에서 수치를 읽어 그리므로 그 산출물이
    # 바뀌면 조용히 구판이 된다. 또 ACS TOC 는 규격이 고정(3.25 x 1.75 in)
    # 이라 실척으로 만들어 두었고, 픽셀 수가 바뀌면 규격이 깨진다.
    acsbad = []
    _acs_figs = []
    for _doc in ('paper_acs/main.tex', 'paper_acs/supporting.tex'):
        if os.path.exists(_doc):
            _acs_figs += re.findall(
                r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', open(_doc).read())
    for _f in sorted(set(_acs_figs)):
        if not os.path.exists(f'paper_acs/{_f}'):
            acsbad.append(f'paper_acs/{_f} 없음')
    # 투고 디렉터리에 쓰이지 않는 도판이 남아 있으면 알린다. 구판 TOC 가
    # 남아 본문 도판처럼 보이던 일이 있었다 (2026-11-02).
    for _p in sorted(glob.glob('paper_acs/*.png')):
        if os.path.basename(_p) not in set(_acs_figs):
            acsbad.append(f'{_p} 는 어느 문서도 쓰지 않는다 — 지우거나 쓸 것')
    _toc, _TOCPX = 'paper_acs/toc_acs.png', (1950, 1050)
    if os.path.exists(_toc):
        try:
            from PIL import Image as _Im
            _sz = _Im.open(_toc).size
            if _sz != _TOCPX:
                acsbad.append(
                    f'{_toc} 픽셀 {_sz} != 규격 {_TOCPX} '
                    f'(3.25 x 1.75 in, 600 dpi) — make_toc_acs.py 재실행')
        except ImportError:
            pass
        _srcs = [p for p in ('outputs/deltaE_premise.json',) if os.path.exists(p)]
        if _srcs and os.path.getmtime(_toc) < max(
                os.path.getmtime(p) for p in _srcs) - MTIME_TOL:
            acsbad.append(f'{_toc} 가 산출물보다 오래됐다 — '
                          f'python3 make_toc_acs.py 로 재생성할 것')
    if acsbad:
        print('[verify] 사전 점검 실패 — ACS 판본 도판.')
        for _m in acsbad:
            print(f'   {_m}')
        return False

    dup = []
    for _p in sorted(glob.glob('deltaE/*_deltaE_pivot.csv')):
        try:
            _d = pd.read_csv(_p)
        except Exception as e:
            dup.append(f'{_p}: 읽기 실패 {e}')
            continue
        # 피벗은 wide 다. 행은 (Agent, Concentration), 열이 Dye1..Dye29.
        # 반복측정이 있었다면 같은 (작용제, 농도)가 두 행으로 나타난다.
        _keys = [c for c in ('Agent', 'Concentration') if c in _d.columns]
        if len(_keys) != 2:
            dup.append(f'{_p}: 기대한 열(Agent, Concentration)이 없다 '
                       f'— 자료 형식이 바뀌었으므로 이 검사가 무력하다')
        elif _d.duplicated(_keys).any():
            _n = int(_d.duplicated(_keys).sum())
            dup.append(f'{_p}: (Agent, Concentration) 중복 {_n}행')
    if dup:
        print('[verify] 사전 점검 실패 — 원자료에 반복측정이 있다.')
        for _m in dup:
            print(f'   {_m}')
        print('[verify] 원고의 "반복측정 없음" 서술과 SUPERSEDED 금지항목'
              ' "어레이 반복측정 주장" 을 함께 재검토할 것.')
        return False
    if missing or figless or nofig or stale:
        print('[verify] 사전 점검 실패 — 필요한 산출물이 없다.')
        for p in missing:
            print(f'   사실표 필수 누락: {p}')
        for p in figless:
            print(f'   도판 생성 필수 누락: {p}')
        for f in nofig:
            print(f'   원고 도판 누락: paper_acs/{f}')
        for f in stale:
            print(f'   도판 주석 등재부에 폐기 항목: {f} '
                  f'(outputs/figure_annotations.json 을 지우고 도판 재생성)')
        print('[verify] sh run_analysis.sh && sh build_paper_acs.sh --figs 를 먼저 '
              '실행하거나, 저장소에서 빠진 파일을 추가할 것.')
        return False
    return True


def load_facts():
    """산출물에서만 읽는다. 여기에 손으로 적은 수치는 없다."""
    f = {}
    # 2026-11-07. ga_q28 (구판 그래픽 초록 주석) 을 폐기하고 TOC 주석으로
    # 옮긴다. TOC 는 독자가 가장 먼저 보는 자리이고 9차에 철회한 주장을
    # 들고 있었다.
    # 2026-11-07 (12차 Step 5). TOC 패널 (a) 가 DMMP 코사인에서 역전 계수로
    # 바뀌었다. toc_cos_* 다섯은 폐기하고 역전 두 사실을 건다. 값의 출처는
    # reversal_audit.json 이고 TOC 주석이 적는 글자와 같은 수다.
    _tc = json.load(open('outputs/deltaE_premise.json'))
    _trv = json.load(open('outputs/reversal_audit.json'))['totals']
    f['toc_rev_uncentred'] = _trv['reversals_uncentred']
    f['toc_rev_centred'] = _trv['reversals_centred']
    _trho = [r['spearman'] for r in _tc['space_test']]
    f['toc_rho_min'] = min(_trho)
    f['toc_rho_max'] = max(_trho)
    f['toc_cells'] = len(_trho)
    NE = json.load(open('outputs/null_match_expectation.json'))['sets']
    f['nullexp_k3'] = NE['incumbent3']['null_mean']
    f['nullexp_k8'] = NE['incumbent_plus_lit8']['null_mean']
    NA = json.load(open('outputs/nearest_agent_bootstrap.json'))['simulants']
    f['na_dmmp_top'] = NA['DMMP']['top_pct']
    f['na_dmmp_g'] = NA['DMMP']['g_series_pct']
    f['na_dcm_g'] = NA['dichloromethane']['g_series_pct']
    f['na_dpgme_top'] = NA['DPGME']['top_pct']
    f['na_dpgme_g'] = NA['DPGME']['g_series_pct']
    L = json.load(open('outputs/fig_label_counts.json'))
    f['fig3_from_panel'] = L['n_from_panel']
    f['fig3_other_labelled'] = L['n_other_labelled']
    f['fig3_unlabelled'] = L['n_unlabelled']
    q = json.load(open('outputs/qsar_tox_summary.json'))
    f['n_train'] = q['n_train']
    f['n_heldout'] = q['n_heldout_cwa']
    f['n_features'] = q['n_features']
    f['ad_threshold'] = q['ad_threshold']
    for row in q['cv']:
        if row['tag'].startswith('RandomForest'):
            f['rf_r2'] = row['r2_mean']
            f['rf_r2_sd'] = row['r2_sd']
            f['rf_rmse'] = row['rmse_mean']
            f['rf_rmse_lo'], f['rf_rmse_hi'] = row['rmse_ci95']
        if row['tag'].startswith('baseline-1'):
            f['b1_rmse'] = row['rmse_mean']
            f['b1_lo'], f['b1_hi'] = row['rmse_ci95']
        if row['tag'].startswith('baseline-2'):
            f['b2_rmse'] = row['rmse_mean']

    h = json.load(open('outputs/h1_summary.json'))
    f['n_candidates'] = h['n_candidates']
    f['spearman_r'] = h['spearman_r']
    f['spearman_lo'], f['spearman_hi'] = h['spearman_ci']
    f['pearson_r'] = h['pearson_r']
    f['pearson_lo'], f['pearson_hi'] = h['pearson_ci']
    f['dcm_aegl3'] = h['dcm_aegl3_mgm3']
    t = h['safety_tiers']
    f['n_relaxed_10'] = t['relaxed_10']['n']
    f['n_relaxed_100'] = t['relaxed_100']['n']
    f['n_incumbent'] = t['incumbent']['n']
    f['n_more_faithful_10'] = t['relaxed_10']['n_more_faithful_than_dcm']
    f['relaxed_bar'] = t['relaxed_10']['bar_mgm3']
    f['fid_cut'] = h['fidelity_cut']
    f['n_screened_distinct'] = h['n_screened_distinct']

    hs = json.load(open('outputs/h1_screen_sensitivity.json'))
    f['sens_kept_n'] = hs['kept']['n_candidates']
    f['sens_kept_rho'] = hs['kept']['spearman_r']

    m = json.load(open('outputs/minority_regression.json'))
    f['n_ache'] = m['n_ache']
    f['n_ache_op'] = m['n_by_group']['organophosphorus_ester']
    f['n_ache_carb'] = m['n_by_group']['carbamate_ester']
    f['low10_err'] = m['low10_err']
    f['low10_n'] = m['low10_n']
    f['high10_err'] = m['high10_err']
    f['phorate_fold'] = m['most_toxic_ache']['fold_too_safe']
    f['malathion_fold'] = m['least_toxic_ache']['fold_too_toxic']

    v = json.load(open('outputs/aegl_parse_validation.json'))
    f['xreader_pairs'] = v['cross_reader_pairs']
    f['xreader_mismatch'] = v['cross_reader_mismatch']
    f['trainset_checked'] = v['trainset_checked']
    f['panel_checked'] = v['panel_checked']
    f['comma_rows'] = v['comma_rows_source']

    T = pd.read_csv('outputs/toxicity_aegl.csv')
    f['panel_with_aegl3'] = int(T.AEGL3_10min_mgm3.notna().sum())
    f['panel_from_footnote'] = int((T.AEGL3_provenance == 'footnote').sum())
    _H = T[T.AEGL3_10min_mgm3.notna()]
    _dv = float(_H[_H.compound_name == 'dichloromethane'].AEGL3_10min_mgm3.iloc[0])
    f['dcm_panel_rank'] = int((_H.AEGL3_10min_mgm3 > _dv).sum()) + 1
    f['_dcm_is_safest'] = f['dcm_panel_rank'] == 1
    f['_panel_safer_than_dcm'] = sorted(
        _H[_H.AEGL3_10min_mgm3 > _dv].compound_name.tolist())
    f['panel_without'] = int(T.AEGL3_10min_mgm3.isna().sum())

    E = pd.read_csv('outputs/qsar_heldout.csv')
    f['heldout_n'] = int(len(E))
    f['heldout_in_domain'] = int(E.in_domain.sum())
    f['heldout_worst_in_domain'] = float(E[E.in_domain].error_log10.max())
    f['heldout_err_min'] = float(E.error_log10.min())
    f['heldout_err_max'] = float(E.error_log10.max())

    C = pd.read_csv('outputs/h1_candidates.csv')
    _d = C[C.label == 'dichloromethane']
    _t = C[C.label == 'Trichloroethylene']
    if len(_d) and len(_t):
        f['tce_ratio'] = float(_d.fidelity_dist.iloc[0] / _t.fidelity_dist.iloc[0])
        f['dcm_above'] = int((C.y_tox > float(_d.y_tox.iloc[0])).sum())
    f['panel_without'] = int(T.AEGL3_10min_mgm3.isna().sum())
    for nm, key in [('dichloromethane', 'dcm'), ('methanol', 'methanol'),
                    ('n-hexane', 'hexane'), ('chloroform', 'chloroform'),
                    ('toluene', 'toluene'), ('acetone', 'acetone'),
                    ('benzene', 'benzene'), ('p-xylene', 'pxylene'),
                    ('cyclosarin (GF)', 'gf'), ('VX', 'vx'),
                    # 2026-11-04 (J-C). SI 의 AEGL-3 순서 문장에 앵커가
                    # 없었다. 주입시험 하나가 숫자를 바꿔 끼우다 그 문장에
                    # 잘못 붙는 사고로 드러났다. 순서를 적는 다섯 값 전부를
                    # 사실로 올린다.
                    ('soman (GD)', 'gd'), ('sarin (GB)', 'gb'),
                    ('tabun (GA)', 'ga')]:
        r = T[T.compound_name == nm]
        if len(r) and r.AEGL3_10min_mgm3.notna().all():
            f[f'aegl3_{key}'] = float(r.AEGL3_10min_mgm3.iloc[0])

    P = pd.read_csv('outputs/toxicity_predicted.csv')
    d = P[P.compound_name == 'DMMP']
    if len(d):
        f['dmmp_ad'] = float(d.ad_dist.iloc[0])

    s = json.load(open('outputs/permutation_sidedness.json'))
    f['p_struct'] = s['A_rdkit']['one_sided']
    f['p_elec'] = s['B_electronic']['one_sided']
    f['p_merged'] = s['C_merged']['one_sided']

    # 롤백 C 의 헤드라인 수치. 2026-10-13 신설. 이 값들이 본문의 주장을 떠받치므로
    # 전량 등재한다. 종전에는 헤드라인이 검증기에 없어서 전파 실패가 통과했다.
    rb = json.load(open('outputs/rollback_c.json'))
    dc = rb['dimensionality_control']
    f['pc3_p28'] = dc['panel28']['pc3_p']
    f['omega_dev'] = rb['omega_dependence']['panel28'][
        'max_deviation_from_definition']
    f['r_lumo_omega'] = rb['omega_dependence']['panel28']['r_lumo_omega']
    f['omega_distrank'] = rb['omega_dependence']['panel28'][
        'dist_rank_correlation_4d_vs_3d']
    f['omega_pc3_var'] = 100.0 * rb['omega_dependence']['panel28'][
        'pc3_cumulative_variance']
    # 2026-10-25. Mantel 다중성 값이 등록돼 있지 않아, SI 가 한 절 안에서
    # 0.112 와 '0.177 에서 0.124 로' 를 동시에 적고도 통과했다. 12셀은
    # 산출물에서, 8셀(병합 제외 구판 가족)은 같은 p 에서 다시 계산한다.
    _mc = rb['deltaE']['mantel']
    f['mantel_bh12'] = float(rb['deltaE']['mantel_bh']['min_q'])

    def _bh_min(ps):
        _p = _np.asarray(sorted(ps), float)
        _m = len(_p)
        _q = _np.minimum.accumulate(
            (_p * _m / _np.arange(1, _m + 1))[::-1])[::-1]
        return float(min(_np.minimum(_q, 1.0)))

    f['mantel_bh8'] = _bh_min([c['p_mantel'] for c in _mc
                               if c['space'] != 'merged'])

    for m in rb['deltaE']['mantel']:
        if m['dataset'] == 'd2' and m['mode'] == 'uv' and \
                m['space'] == 'structural':
            f['mantel_d2uv_struct'] = m['p_mantel']
            f['mantel_d2uv_struct_rho'] = m['spearman']
            f['naive_d2uv_struct'] = m['p_pairs_as_independent']
        if m['dataset'] == 'd2' and m['mode'] == 'uv' and m['space'] == 'merged':
            f['mantel_d2uv_merged'] = m['p_mantel']
            f['mantel_d2uv_merged_rho'] = m['spearman']
        if m['dataset'] == 'd2' and m['mode'] == 'normal' and \
                m['space'] == 'electronic3':
            f['mantel_d2nor_elec'] = m['p_mantel']
            f['mantel_d2nor_elec_rho'] = m['spearman']
            f['naive_d2nor_elec'] = m['p_pairs_as_independent']
    # 2026-10-24 신설. 중복 채널 축약은 본문이 '23 effective channels' 라고
    # 쓰는 새 수치다. 손으로 적은 숫자가 아니라 산출물에서 온 것임을 묶는다.
    # 화합물 단위 검정력. deltaE_power.py 와 같은 식을 쓴다.
    from scipy import stats as _st
    from scipy.optimize import brentq as _bq
    def _pw(r, n, a=0.05):
        z = _np.arctanh(r) * _np.sqrt(n - 3)
        zc = _st.norm.ppf(1 - a / 2)
        return _st.norm.cdf(z - zc) + _st.norm.cdf(-z - zc)
    # 2026-10-24 2차 교정. Fisher z 근사(0.924)는 실제로 돌린 정확 Mantel 의
    # 검정력을 과소평가한다. 시뮬레이션 값(mantel_power.py)이 정본이다.
    f['mde_rho_compound'] = float(_bq(lambda r: _pw(r, 6) - 0.8, .01, .999))
    import os as _os
    if _os.path.exists('outputs/mantel_power.json'):
        _mp = json.load(open('outputs/mantel_power.json'))
        f['mde_rho_mantel'] = _mp['mde_rho_mantel']
        f['mde_rho_mantel_bh'] = _mp['mde_rho_mantel_bh']
        f['mantel_power_perm'] = _mp['n_perm']
        # 2026-11-06. GATE 8 9차 C3-1. SI 가 'type I rate ... at 0.049' 라고
        # 적고 있었는데 이 스크립트에 귀무 조건이 없었다. 수행하지 않은 결과를
        # 적은 것이다. mantel_power.py 에 run_null() 을 신설해 실제로 돌렸고
        # 그 값을 여기에 묶는다. 원고가 다시 떠다니면 FAIL 로 잡힌다.
        f['mantel_type_i'] = _mp['null']['type_i_rate']
        f['mantel_type_i_nrep'] = _mp['null']['n_rep']
    # 2026-10-26. 다중성 표(tab:family) 28칸을 전부 등록한다. 종전에는 몇 칸만
    # 걸려 있어 BH 0.185, meff 5.00, 사전지정 0.046 이 구판인 채 통과했다.
    _fam = json.load(open('outputs/stat_supplement.json'))['family']
    _wy = json.load(open('outputs/stat_supplement2.json'))['westfall_young']
    _pf = json.load(open('outputs/pesticide_family.json'))
    _cols = ['p35', 'p35lit', 'p28']
    _keys = ['panel35_incumbent3', 'panel35_incumbent_plus_lit8',
             'panel28_incumbent3']
    for _lab, _k in (('minp', 'min_p'), ('bh18', 'bh_all_min_q'),
                     ('bh6', 'bh_by_space_min_q'), ('meff', 'meff'),
                     ('bonf', 'meff_adjusted_min_p'),
                     ('pre', 'prespecified_bh_min_q')):
        for _c, _kk in zip(_cols, _keys):
            f[f'fam_{_lab}_{_c}'] = float(_fam[_kk][_k])
        f[f'fam_{_lab}_pest'] = float(_pf[_k])
    for _c, _kk in zip(_cols, _keys):
        f[f'fam_wy_{_c}'] = float(_wy[_kk]['wy_adjusted_p'])
    f['fam_wy_pest'] = float(_pf['westfall_young_adjusted_p'])

    # 2026-10-27. 패널 축까지 포함한 72검정 통합 BH.
    if os.path.exists('outputs/multiplicity_across_panels.json'):
        _ap = json.load(open('outputs/multiplicity_across_panels.json'))
        f['ap_ntests'] = _ap['n_tests']
        _m = _ap['min_q_by_config']
        f['ap_p28'] = _m['panel28_designated3']['q']
        f['ap_p35'] = _m['panel35_designated3']['q']
        f['ap_p35lit'] = _m['panel35_designated_plus_literature8']['q']
        f['ap_pest'] = _m['pesticide10']['q']

    # 2026-10-27. 부분공간 대조의 귀무 보정. 본문이 36.0% 를 '평범한 draw' 로
    # 읽었던 자리를 뒤집는 수치이므로 전부 건다.
    # 2026-11-04 (G3). SI 의 두 '수' 가 산출물에 없거나 앵커가 없었다.
    # 207종 중 일곱만 더 안전(경계값인 DCM 자신은 분면 8에 포함되므로
    # 둘은 다른 수다)과 한 자릿수 아래 여덟 후보.
    _h1s = json.load(open('outputs/h1_summary.json'))
    f['h1_safer7'] = _h1s['n_strictly_safer_than_incumbent']
    f['h1_relaxed10'] = _h1s['safety_tiers']['relaxed_10']['n']
    f['h1_ncand'] = _h1s['n_candidates']
    if os.path.exists('outputs/subspace_null.json'):
        _sn = json.load(open('outputs/subspace_null.json'))
        # 2026-11-04 (작업 D). 순위 자체에 검사가 없었다. 37좌표 구판의
        # 'fifth' 를 되살리는 주입이 통과했다. 서수를 사실로 올린다.
        f['sn_rank'] = _sn['rank_of_actual']
        f['sn_actual'] = 100 * _sn['actual_designation_frac']
        f['sn_mean'] = 100 * _sn['null_mean_frac']
        f['sn_median'] = 100 * _sn['null_median_frac']
        f['sn_p95'] = 100 * _sn['null_p95_frac']
        f['sn_max'] = 100 * _sn['null_max_frac']
        f['sn_rank'] = _sn['rank_of_actual']
        f['sn_p'] = _sn['p_of_actual']
        f['sn_nsub'] = _sn['n_subspaces']
        f['sn_ndesig'] = _sn['n_designations']

    # 2026-10-26. 형태이성질체 민감도. 3D 서술자 9개가 모든 거리에 들어가므로
    # 임베딩 시드에 대한 폭을 원고가 싣는다. 두 처리(앙상블/단일)를 모두 건다.
    _cs = {}
    for _m in ('ensemble', 'single'):
        _fp = f'outputs/conformer_sensitivity_{_m}.json'
        if not _os.path.exists(_fp):
            continue
        _d = json.load(open(_fp))
        _tab = [r['panel28']['d_norm']['tabun (GA)'] for r in _d['runs']]
        _nm = {r['panel28']['n_match'] for r in _d['runs']}
        _n35 = {r['panel35']['n_match'] for r in _d['runs']}
        _pp = [r['panel28']['perm_p_merged'] for r in _d['runs']]
        _cs[_m] = dict(lo=min(_tab), hi=max(_tab), seeds=len(_d['runs']),
                       nm=_nm, n35=_n35, plo=min(_pp), phi=max(_pp))
    if 'ensemble' in _cs and 'single' in _cs:
        _e, _s1 = _cs['ensemble'], _cs['single']
        # 2026-11-02 (T2-8). 종전 관문은 두 처리 모두 매칭 4건이고 시드에
        # 흔들리지 않는다는 것을 전제로 박아 두었다. 결정적 상관필터로
        # 교체하자 그 전제가 거짓이 되었다. 앙상블은 시드 6개에서 3건과
        # 4건으로 갈리고(타분 1.9943~2.0106), 단일은 6개 모두 3건이다
        # (타분 2.0172~2.0273). 즉 타분은 역치 2 의 양쪽 0.03 안에 있고
        # 매칭 여부가 임베딩 시드의 함수다. 관문이 제 몫을 했다.
        #
        # 고정해야 하는 것은 '4건' 이 아니라 주 패널의 0건이다. 그쪽이
        # 흔들리면 멈추고, 28종의 분포는 사실표에 실어 원고가 실측 폭을
        # 인용하게 한다. 28종 매칭 수가 {3,4} 밖으로 나가면 그것은 다른
        # 변경이므로 역시 멈춘다.
        if _e['n35'] != {0} or _s1['n35'] != {0}:
            raise SystemExit('[verify] 중단: 형태이성질체 민감도에서 주 패널 '
                             f"매칭 수가 0 이 아니다 (앙상블 {sorted(_e['n35'])}, "
                             f"단일 {sorted(_s1['n35'])}). 원고의 헤드라인을 "
                             '먼저 고쳐야 한다.')
        if not (_e['nm'] | _s1['nm']) <= {3, 4}:
            raise SystemExit('[verify] 중단: 형태이성질체 민감도의 28종 매칭 '
                             f"수가 {sorted(_e['nm'] | _s1['nm'])} 로 "
                             '{3,4} 를 벗어났다. 원고의 서술을 먼저 고칠 것.')
        f['conf_seeds'] = _e['seeds']
        f['conf_ens_lo'] = _e['lo']
        f['conf_ens_hi'] = _e['hi']
        f['conf_sgl_lo'] = _s1['lo']
        f['conf_sgl_hi'] = _s1['hi']
        f['conf_p_hi'] = _s1['phi']
        f['conf_ens_match_lo'] = min(_e['nm'])
        f['conf_ens_match_hi'] = max(_e['nm'])
        f['conf_sgl_match'] = max(_s1['nm'])
        f['conf_ens_n4'] = sum(
            1 for r in json.load(
                open('outputs/conformer_sensitivity_ensemble.json'))['runs']
            if r['panel28']['n_match'] == 4)

    # 2026-11-03 (T1-14). 희소 지시자 인코딩. 본문 2.4 와 SI 가 인용한다.
    if _os.path.exists('outputs/sparse_indicators.json'):
        _sp = json.load(open('outputs/sparse_indicators.json'))
        _V = {v['variant']: v for v in _sp['variants']}
        f['sp_n_singleton'] = len(_sp['singleton_indicators'])
        f['sp_pedestal'] = _sp['pedestal_per_singleton']
        f['sp_base_dims'] = _sp['baseline_dims']
        f['sp_tabun_share'] = 100.0 * _sp['share_singletons'][
            'tabun (GA) -> DMMP']['share']
        f['sp_tabun_base'] = _V['baseline']['tabun_d_norm']
        f['sp_tabun_drop1'] = _V['drop_binary_nonzero_le_1']['tabun_d_norm']
        f['sp_match_base'] = _V['baseline']['n_match']
        f['sp_match_drop1'] = _V['drop_binary_nonzero_le_1']['n_match']
        f['sp_match_drop3'] = _V['drop_binary_nonzero_le_3']['n_match']
        f['sp_match_lowcard'] = _V['drop_cardinality_le_5']['n_match']
        f['sp_p_base'] = _V['baseline']['p_one_sided']
        f['sp_p_lowcard'] = _V['drop_cardinality_le_5']['p_one_sided']
        # 2026-11-07 (10차 R3). 종전에는 세 변종만 키였다. sparse_indicators
        # 가 일곱 전부를 돌므로 일곱 전부에 앵커를 붙인다. SI 가 여섯 수를
        # 나란히 적으므로 자리별로 따로 검사한다.
        _ms = {r['variant']: r for r in _sp['mantel_stability']}
        f['sp_mantel_q_base'] = _ms['baseline']['min_q']
        f['sp_mantel_q_drop1'] = _ms['drop_binary_nonzero_le_1']['min_q']
        f['sp_mantel_q_drop2'] = _ms['drop_binary_nonzero_le_2']['min_q']
        f['sp_mantel_q_drop3'] = _ms['drop_binary_nonzero_le_3']['min_q']
        f['sp_mantel_q_drop4'] = _ms['drop_binary_nonzero_le_4']['min_q']
        f['sp_mantel_q_drop5'] = _ms['drop_binary_nonzero_le_5']['min_q']
        f['sp_mantel_q_lowcard'] = _ms['drop_cardinality_le_5']['min_q']
        f['sp_mantel_q_max'] = max(r['min_q'] for r in _sp['mantel_stability'])
        f['sp_mantel_nvar'] = len(_sp['mantel_stability'])

    # 전제 검정의 쌍별 위반. 초록이 '평균' 표현을 쓰는 근거이자 그 한계.
    if _os.path.exists('outputs/deltaE_premise.json'):
        _pm = json.load(open('outputs/deltaE_premise.json'))['discrimination']
        f['prem_pairwise_viol'] = sum(r['n_pairwise_violations'] for r in _pm)
        f['prem_n_compound_cells'] = sum(r['n_compounds'] for r in _pm)
        f['prem_p_lo'] = min(r['p_exact_signflip'] for r in _pm)
        f['prem_p_hi'] = max(r['p_exact_signflip'] for r in _pm)
        f['prem_dz_lo'] = min(r['cohens_dz'] for r in _pm)
        f['prem_dz_hi'] = max(r['cohens_dz'] for r in _pm)

    # 타분 정규화거리. 신판 정본과 구판 대조값 양쪽을 사실표에 둔다.
    if _os.path.exists('outputs/match_by_panel.json'):
        _mbp = json.load(open('outputs/match_by_panel.json'))
        for _a in _mbp['panel28_incumbent3']['agents']:
            if _a['agent'] == 'tabun (GA)':
                f['t3_GA_full'] = _a['d_norm_full']
    # 구판 필터의 '원래 열 순서' 결과. order_invariance 는 분포만 담으므로
    # (최대값은 치환으로 나온 다른 값이다) match_rule_audit 의 legacy 팔에서
    # 읽는다. 이 값은 SI 열순서 절이 대조값으로 인용한다.
    if _os.path.exists('outputs/match_rule_audit.json'):
        _mra = json.load(open('outputs/match_rule_audit.json'))
        for _a in _mra['panel28_incumbent3']['legacy/wide']['agents']:
            if _a['agent'] == 'tabun (GA)':
                f['t3_GA_legacy'] = _a['d_norm']

    # 2026-11-03. 위 앵커들이 읽는 사실. 전부 산출물에서만 뽑는다.
    if _os.path.exists('outputs/rollback_c.json'):
        _rc = json.load(open('outputs/rollback_c.json'))['deltaE']
        f['mantel_bh12_main'] = _rc['mantel_bh']['min_q']
        for _m in _rc['mantel']:
            if _m['dataset'] == 'd2' and _m['mode'] == 'uv' \
                    and _m['space'] == 'structural':
                f['mantel_d2uv_struct_p'] = _m['p_mantel']
            if _m['dataset'] == 'd2' and _m['mode'] == 'normal' \
                    and _m['space'] == 'electronic3':
                f['mantel_d2nor_elec_p'] = _m['p_mantel']
        f['mantel_rho_max_abs'] = max(abs(x['spearman'])
                                      for x in _rc['mantel'])
        _d1 = [x for x in _rc['mantel'] if x['dataset'] == 'd1']
        if _d1:
            f['d1_rho_max_abs'] = max(abs(x['spearman']) for x in _d1)
            f['d1_p_min'] = min(x['p_mantel'] for x in _d1)
    # 농도와 응답크기의 순위상관. 모든 셀에서 같은 값이어야 한다.
    if _os.path.exists('outputs/deltaE_noise.json'):
        # 2026-11-06 (9차 R3 6-4). 종전에는 concentration_vs_magnitude
        # 를 읽었는데 그 표는 DMMP 가 있는 두 셀뿐이다. 본문은 '모든
        # 런-채널 셀' 이라고 적으므로 여섯 셀 전수 표를 읽는다.
        _dnj = json.load(open('outputs/deltaE_noise.json'))
        _dn = _dnj.get('concentration_vs_norm_all_cells') or _dnj[
            'concentration_vs_magnitude']
        _rs = {round(r['rho_conc_vs_norm'], 6) for r in _dn}
        if len(_rs) == 1:
            f['conc_norm_rho'] = _rs.pop()

    # 중심화 바닥 여섯 셀. T1-13 이 본문·SI 양쪽에 싣는다.
    if _os.path.exists('outputs/rollback_c.json'):
        _cf = json.load(open('outputs/rollback_c.json'))['deltaE'][
            'centring_floor']
        f['floor_gap_lo'] = min(x['gap'] for x in _cf)
        f['floor_gap_hi'] = max(x['gap'] for x in _cf)
        _d3n = [x for x in _cf if x['dataset'] == 'd3'
                and x['mode'] == 'normal']
        if _d3n:
            f['floor_d3normal_obs'] = _d3n[0]['observed_mean_cosine']
            f['floor_d3normal_gap'] = _d3n[0]['gap']

    # 2026-10-25. 전자공간을 B3LYP 궤도에너지로 갈아끼운 대조. 심사
    # 시늄레이션이 제기하고 재현으로 확인된 축이다. 기준 재현(GFN2)이 정본
    # 0.275 와 1/11 을 내는 것까지 같은 파일에 들어 있다.
    if _os.path.exists('outputs/dft_electronic_space.json'):
        _es = json.load(open('outputs/dft_electronic_space.json'))
        f['dft_es_gfn2_p'] = float(_es['gfn2']['p_one_sided'])
        f['dft_es_gfn2_match'] = int(_es['gfn2']['n_match'])
        f['dft_es_b3lyp_p'] = float(_es['b3lyp']['p_one_sided'])
        f['dft_es_b3lyp_match'] = int(_es['b3lyp']['n_match'])
        # 2026-10-27. 쌍극자까지 B3LYP 로 바꾼 3좌표 배치. 종전에는 쌍극자를
        # GFN2 로 남겨 두고 "어떤 통상 수준에서도 유의하지 않다"고 적었는데,
        # mf.dip_moment() 는 이미 수렴한 SCF 에서 공짜로 나오므로 남길 이유가
        # 없었다. 세 좌표 전부 바꾸면 p = 0.0368 이고 그 서술은 거짓이다.
        f['dft_es_full_p'] = float(_es['b3lyp_full']['p_one_sided'])
        f['dft_es_full_match'] = int(_es['b3lyp_full']['n_match'])
        # 병합공간 대조. GFN2 쪽이 정본(39좌표, 4/11, 0.0250)을 재현하는 것이
        # 이 비교의 관문이다.
        f['dft_mg_gfn2_p'] = float(_es['merged_gfn2']['p_one_sided'])
        f['dft_mg_gfn2_match'] = int(_es['merged_gfn2']['n_match'])
        f['dft_mg_b3lyp_p'] = float(_es['merged_b3lyp']['p_one_sided'])
        f['dft_mg_b3lyp_match'] = int(_es['merged_b3lyp']['n_match'])
    # 2026-11-02 (GATE 8 7차). SI 가 쌍극자 재계산을 설명하면서 두 방법의
    # 일치도를 수치로 내놓지 않았다. 산출물에서 직접 센다.
    if (os.path.exists('outputs/dft_dipole.csv')
            and os.path.exists('descriptor_cache_xtb.csv')):
        _dd = pd.read_csv('outputs/dft_dipole.csv').set_index('name')
        _xx = pd.read_csv('descriptor_cache_xtb.csv').set_index('compound_name')
        _g = _xx['dipole_moment_D'].reindex(_dd.index)
        _b = _dd['dipole_D']
        _ok = _g.notna() & _b.notna()
        if _ok.sum() > 2:
            f['dip_r'] = float(_g[_ok].corr(_b[_ok]))
            f['dip_mae'] = float((_b[_ok] - _g[_ok]).abs().mean())
            _dv = (_b[_ok] - _g[_ok]).abs()
            f['dip_maxdev'] = float(_dv.max())
            f['dip_dmmp_gfn2'] = float(_g['DMMP'])
            f['dip_dmmp_b3lyp'] = float(_b['DMMP'])
        # 2026-11-02. 종전에는 정본 값 (39좌표, 4/11) 을 상수로 박아 두었다.
        # T2-1 에서 상관필터를 결정적 규칙으로 교체해 정본이 (38, 3) 으로
        # 바뀌자 이 관문이 멈췄다. 관문의 뜻은 'DFT 대조의 GFN2 팔이 정본
        # 병합공간과 같은 공간인가' 이므로, 상수 대신 살아 있는 정본 산출물과
        # 대조한다. 그래야 다음에 정본이 또 바뀌어도 관문이 뜻을 유지한다.
        if _os.path.exists('outputs/match_by_panel.json'):
            _mb = json.load(open('outputs/match_by_panel.json'))
            _c = _mb['panel28_incumbent3']
            if (_es['merged_gfn2']['n_coords'] != _c['n_dims']
                    or _es['merged_gfn2']['n_match'] != _c['n_match_wide']):
                raise SystemExit(
                    '[중단] 병합공간 GFN2 대조가 정본'
                    f"({_c['n_dims']}좌표, {_c['n_match_wide']}/11)을 재현하지 "
                    f"못한다: {_es['merged_gfn2']['n_coords']}좌표, "
                    f"{_es['merged_gfn2']['n_match']}/11")

    # 2026-10-25. 병리지수 제거 대조. 구판은 omega 가 아직 든 4좌표 전자공간에서
    # 돌아간 산출물이었고 본문은 그 값을 '모든 분석을 반복했다'며 인용했다.
    if _os.path.exists('outputs/permutation_nopatho.json'):
        _nop = json.load(open('outputs/permutation_nopatho.json'))
        if _nop['B_electronic']['n_kept'] != 3:
            raise SystemExit('[verify] 중단: permutation_nopatho 의 전자공간이 '
                             f"{_nop['B_electronic']['n_kept']}좌표다. omega 제거 "
                             '이후 재실행되지 않은 산출물이다.')
        f['nopatho_struct_p'] = float(_nop['A_rdkit']['p_one_sided'])
        f['nopatho_elec_p'] = float(_nop['B_electronic']['p_one_sided'])
        f['nopatho_merged_p'] = float(_nop['C_merged']['p_one_sided'])

    # 2026-10-25 신설. DFT 대조 수치가 SI 에 들어갔다. GFN2 의 omega 축에서
    # 끌어낸 순서가 방법 종속임을 보이는 값들이므로 산출물에 묶는다.
    if _os.path.exists('outputs/dft_vs_gfn2.csv'):
        import pandas as _pd
        _dv = _pd.read_csv('outputs/dft_vs_gfn2.csv')
        _ag = _dv[_dv.compound_class == 'real_agent']
        _it = _dv[_dv.compound_class == 'interferent']
        _gb = float(_ag[_ag.name.str.startswith('sarin')].omega_eV.iloc[0])
        _dm = float(_dv[_dv.name == 'DMMP'].omega_eV.iloc[0])
        f['dft_omega_min'] = float(_dv.omega_eV.min())
        f['dft_omega_max'] = float(_dv.omega_eV.max())
        f['gfn2_omega_min'] = float(_dv.gfn2_w.min())
        f['gfn2_omega_max'] = float(_dv.gfn2_w.max())
        f['dft_ps'] = float(_dv[_dv.name.str.startswith('chloropicrin')].omega_eV.iloc[0])
        f['dft_cg'] = float(_dv[_dv.name.str.startswith('phosgene')].omega_eV.iloc[0])
        f['dft_gb'] = _gb
        f['dft_dmmp'] = _dm
        f['dft_itf_above_gb'] = int((_it.omega_eV > _gb).sum())
        f['dft_n_dft'] = int(len(_dv))
        _lo = _ag.sort_values('omega_eV').head(3)
        f['dft_lo1'] = float(_lo.omega_eV.iloc[0])
        f['dft_lo2'] = float(_lo.omega_eV.iloc[1])
        f['dft_lo3'] = float(_lo.omega_eV.iloc[2])
        f['dft_rho_omega'] = float(_dv[['gfn2_w', 'omega_eV']]
                                   .corr(method='spearman').iloc[0, 1])
        f['dft_rho_lumo'] = float(_dv[['E_LUMO_eV', 'lumo_eV']]
                                  .corr(method='spearman').iloc[0, 1])

    # 2026-10-24. 다른 세션이 올린 AsCl3 민감도 표(SI)가 검증기에 등록돼 있지
    # 않았다. 그 표의 28종 행은 이번 정확종속 규칙으로 3/11 에서 4/11 로
    # 움직였고, 등록돼 있었다면 기계가 잡았을 자리다.
    if _os.path.exists('outputs/match_by_panel_ascl3.json'):
        _ac = json.load(open('outputs/match_by_panel_ascl3.json'))
        _AC = {'panel28_incumbent3': 'a28', 'panel35_incumbent3': 'a35',
               'panel35_incumbent_plus_lit8': 'a35l',
               'panel36_incumbent3': 'a36', 'panel36_incumbent_plus_lit': 'a36l'}
        for _k, _t in _AC.items():
            if _k in _ac:
                f[f'ascl3_{_t}_match'] = _ac[_k]['n_match']
                _L = [x for x in _ac[_k]['agents'] if 'lewisite' in x['agent']]
                if _L:
                    f[f'ascl3_{_t}_dist'] = _L[0]['d_norm']

    # 2026-10-24 신설. 하단꼬리와 중첩 조건부 귀무는 초록까지 올라간 새 수치다.
    _nm = json.load(open('outputs/null_match_expectation.json'))
    f['null3_lower'] = _nm['sets']['incumbent3']['p_at_most_observed']
    f['null8_nested_mean'] = _nm['nested_null_8']['null_mean']
    f['null8_nested_p'] = _nm['nested_null_8']['p_at_least_observed']
    _ch = rb['deltaE'].get('channels', {})
    if 'd2.uv' in _ch:
        f['n_eff_d2uv'] = _ch['d2.uv']['n_effective']
        f['n_declared_dyes'] = _ch['d2.uv']['n_declared']
    if 'd3.uv' in _ch:
        f['n_eff_d3uv'] = _ch['d3.uv']['n_effective']
    for fl in rb['deltaE']['centring_floor']:
        if fl['dataset'] == 'd1' and fl['mode'] == 'uv':
            f['floor_d1uv_obs'] = fl['observed_mean_cosine']
            f['floor_d1uv_floor'] = fl['algebraic_floor']
    for rk in rb['dmmp_ranking']:
        if rk['dataset'] == 'd1' and rk['mode'] == 'uv':
            for e in rk['ranking']:
                if e['agent'] == 'GA':
                    f['dmmp_ga_uv'] = e['cosine']
                if e['agent'] == 'GB':
                    f['dmmp_gb_uv'] = e['cosine']

    # 2026-11-01. ACS 판본(paper_acs/main.tex)이 싣는 응답공간 수치를 등록한다.
    # 종전에는 elsarticle 판 문구에만 패턴이 걸려 있어, 같은 양을 다르게 쓴
    # ACS 문장이 무검사로 자랐다. 이 원고에서 반복해 터진 '문구별 검사' 결함이다.
    _pre = json.load(open('outputs/deltaE_premise.json'))
    for _d in _pre['dmmp']:
        if _d['mode'] == 'uv':
            f['dmmp_gd_uv'] = _d['cos_GD']
            f['dmmp_gf_uv'] = _d['cos_GF']
            f['dmmp_vx_uv'] = _d['cos_VX']
            f['dmmp_g4mean_uv'] = _d['g4_series_mean']
            f['dmmp_gmean_uv_3'] = _d['g_series_mean']
            f['gvx_uv'] = _d['g_minus_vx']
            f['gvx_uv_lo'] = _d['ci_lo']
            f['gvx_uv_hi'] = _d['ci_hi']
        else:
            f['dmmp_g4mean_nor'] = _d['g4_series_mean']
            f['dmmp_gmean_nor_3'] = _d['g_series_mean']
            f['gvx_nor'] = _d['g_minus_vx']
            f['gvx_nor_lo'] = _d['ci_lo']
            f['gvx_nor_hi'] = _d['ci_hi']
    _ext = json.load(open('outputs/deltaE_extended.json'))
    for _r in _ext['concentration']:
        if _r['mode'] == 'uv' and _r['agent'] == 'GB':
            f['conc_gb_lo'] = _r['cosines'][0]
            f['conc_gb_hi'] = _r['cosines'][-1]
        if _r['mode'] == 'uv' and _r['agent'] == 'GF':
            f['conc_gf_lo'] = abs(_r['cosines'][0])
            f['conc_gf_hi'] = abs(_r['cosines'][-1])
    # 추세 10건의 BH 최소 q. 리터럴로 두면 재분석 때 조용히 구판이 된다.
    _tp = sorted(x['trend_p_exact'] for x in _ext['concentration']
                 if 'trend_p_exact' in x)
    _m = len(_tp)
    _q = [min(1.0, pv * _m / (i + 1)) for i, pv in enumerate(_tp)]
    for _i in range(_m - 2, -1, -1):
        _q[_i] = min(_q[_i], _q[_i + 1])
    f['trend_bh_q'] = round(min(_q), 3)
    f['trend_monotone'] = sum(1 for x in _ext['concentration']
                              if abs(x.get('trend_rho', 0)) == 1.0)
    _M = _ext['matching']
    f['rank_obs'] = _M['n_rank1_observed']
    f['rank_exp'] = _M['n_rank1_expected']
    f['rank_p'] = _M['p_binom_one_sided']
    for _c in _ext['dyes']['controls']:
        f['dye_null_' + _c['mode']] = _c['null_half_mean']
    f['fisher_minq'] = _ext['combined']['min_q']
    f['fisher_cells'] = _ext['combined']['n_tests']

    f['p_struct_fig'] = f['p_struct']
    f['p_elec_fig'] = f['p_elec']
    f['p_merged_fig'] = f['p_merged']

    f['mantel_minq'] = rb['deltaE']['mantel_bh']['min_q']
    f['mantel_cells'] = rb['deltaE']['mantel_bh']['n_tests']
    f['mantel_nominal'] = sum(1 for _m in rb['deltaE']['mantel']
                              if _m['p_mantel'] < 0.05)

    ss = json.load(open('outputs/stat_supplement.json'))
    w = ss['whitening']['panel28']
    f['rand3_frac28'] = 100.0 * w['raw']['frac_below_05']
    f['rand3_med28'] = w['raw']['quantiles']['q50']
    f['white_frac28'] = 100.0 * w['whitened']['frac_below_05']
    # 2026-11-03 (2차). 쌍독립 가정 p 와 pooled 구조공간 rho 가 사실표에
    # 없어서, 그 자리의 오기를 주입해도 검증기가 통과했다.
    import pandas as _pdx
    _st = _pdx.read_csv('outputs/deltaE_space_test.csv')
    for _m, _sp, _k in [('uv', 'structural', 'pairs_d2uv_struct'),
                        ('uv', 'merged', 'pairs_d2uv_merged'),
                        ('normal', 'electronic', 'pairs_d2nor_elec')]:
        _r = _st[(_st.dataset == 'd2') & (_st['mode'] == _m)
                 & (_st.space == _sp)]
        f[_k] = float(_r.p.iloc[0])
    _pl = _pdx.read_csv('outputs/deltaE_pooled.csv')
    for _m, _tag in [('uv', 'uv'), ('normal', 'nor')]:
        _r = _pl[(_pl['mode'] == _m) & (_pl.space == 'structural')]
        f[f'pooled_struct_{_tag}'] = float(_r.pooled_rho.iloc[0])
        f[f'pooled_struct_{_tag}_lo'] = float(_r.ci_lo.iloc[0])
        f[f'pooled_struct_{_tag}_hi'] = float(_r.ci_hi.iloc[0])
    f['white_med28'] = w['whitened']['quantiles']['q50']
    s2 = json.load(open('outputs/stat_supplement2.json'))
    wy = s2['westfall_young']
    f['wy35'] = wy['panel35_incumbent3']['wy_adjusted_p']
    f['wy35lit'] = wy['panel35_incumbent_plus_lit8']['wy_adjusted_p']
    f['wy28'] = wy['panel28_incumbent3']['wy_adjusted_p']
    bt = s2['distance_bootstrap']
    f['boot_inc_med'] = bt['incumbent3']['median']
    f['boot_inc_lo'] = bt['incumbent3']['ci'][0]
    f['boot_inc_hi'] = bt['incumbent3']['ci'][1]
    f['boot_inc_zero'] = 100.0 * bt['incumbent3']['frac_zero']
    f['boot_inc_le1'] = 100.0 * bt['incumbent3']['frac_le1']
    f['boot_lit_med'] = bt['incumbent_plus_lit8']['median']
    f['boot_diff_med'] = bt['difference']['median']
    f['boot_diff_pos'] = round(bt['difference']['frac_positive'] * bt['n_boot'])
    _r3 = bt['incumbent3']['retention']
    _r8 = bt['incumbent_plus_lit8']['retention']
    f['ret_gb3'] = 100.0 * _r3.get('sarin (GB)', 0)
    f['ret_gd3'] = 100.0 * _r3.get('soman (GD)', 0)
    f['ret_gf3'] = 100.0 * _r3.get('cyclosarin (GF)', 0)
    f['ret_ga8'] = 100.0 * _r8.get('tabun (GA)', 0)

    tf = ss['threshold_free']
    f['unbounded35'] = tf['panel35_incumbent3']['n_match_unbounded']
    f['null_exp35'] = tf['panel35_incumbent3']['null_expectation_outside_closer']
    f['null_exp28'] = tf['panel28_incumbent3']['null_expectation_outside_closer']

    fam = ss['family']['panel28_incumbent3']
    f['q_by_space28'] = fam['bh_by_space_min_q']
    f['q_prespec28'] = fam['prespecified_bh_min_q']
    f['meff28'] = fam['meff']
    f['meff_p28'] = fam['meff_adjusted_min_p']

    md = json.load(open('outputs/mde.json'))
    f['mde_delta35'] = md['panel35_incumbent3']['mde_delta']
    f['mde_gap_median'] = md['panel35_incumbent3']['gap_median']
    f['mde_gap_min'] = md['panel35_incumbent3']['gap_min']
    f['mde_p0_35'] = md['panel35_incumbent3']['p_at_zero']
    f['mde_p0_28'] = md['panel28_incumbent3']['p_at_zero']

    mb = json.load(open('outputs/match_by_panel.json'))
    f['match35_inc'] = mb['panel35_incumbent3']['n_match']
    f['match35_lit'] = mb['panel35_incumbent_plus_lit8']['n_match']
    f['match28_inc'] = mb['panel28_incumbent3']['n_match']
    _lit = {a['agent']: a for a in mb['panel35_incumbent_plus_lit8']['agents']}
    _inc = {a['agent']: a for a in mb['panel35_incumbent3']['agents']}
    f['dfp_gb'] = _lit['sarin (GB)']['d_norm']
    f['dfp_gd'] = _lit['soman (GD)']['d_norm']
    f['dfp_gf'] = _lit['cyclosarin (GF)']['d_norm']
    f['dmmp_gb_35'] = _inc['sarin (GB)']['d_norm']
    f['dmmp_gd_35'] = _inc['soman (GD)']['d_norm']
    f['dmmp_gf_35'] = _inc['cyclosarin (GF)']['d_norm']
    f['cees_hd'] = _lit['sulfur mustard (HD)']['d_norm']

    # 2026-10-17. 표 3 (tab:match) 은 이 산출물을 쓰면서도 검증기에 등록된 적이
    # 없었다. 144건 검사를 통과하는 동안 11행 중 9행이 구판 값이었다. 동료
    # 검토에서 지적돼 전수 대조로 돌린다. 같은 이유로 res:match 산문의
    # 공간별 수치도 전부 사실표에 올린다.
    _p28 = {a['agent']: a for a in mb['panel28_incumbent3']['agents']}
    _T3 = {'GB': 'sarin (GB)', 'GD': 'soman (GD)',
           'CK': 'cyanogen chloride (CK)', 'GF': 'cyclosarin (GF)',
           'AC': 'hydrogen cyanide (AC)', 'CG': 'phosgene (CG)',
           'HD': 'sulfur mustard (HD)', 'GA': 'tabun (GA)',
           'L': 'lewisite 1 (L)', 'PS': 'chloropicrin (PS)', 'VX': 'VX'}
    for _code, _nm in _T3.items():
        f[f't3_{_code}'] = _p28[_nm]['d_norm']

    # 공간별 주성분 분산. 본문 4.1절 수치가 파이프라인 밖 폐기본에서 오던 것을
    # 2026-10-24 에 끊었다. 검증기가 이 값을 들고 있지 않아 못 잡았다.
    _pf = json.load(open('outputs/pca_facts.json'))
    f['pc1_struct28'] = round(100 * _pf['panel28']['A_rdkit']['pc1'], 2)
    f['pc12_elec28'] = round(100 * _pf['panel28']['B_electronic']['pc1_pc2'], 2)
    f['pc1_struct35'] = round(100 * _pf['panel35']['A_rdkit']['pc1'], 2)
    f['pc12_elec35'] = round(100 * _pf['panel35']['B_electronic']['pc1_pc2'], 2)
    f['pc12_merged28'] = round(100 * _pf['panel28']['C_merged']['pc1_pc2'], 2)
    f['pc12_merged35'] = round(100 * _pf['panel35']['C_merged']['pc1_pc2'], 2)

    # GFN2 가상궤도 에너지의 범위. 한계 절이 이 수치로 'Koopmans 로 읽지
    # 말라'를 건다 (2026-10-23 신설).
    _x35 = pd.read_csv('descriptor_cache_xtb_35.csv')
    f['lumo_min'] = float(_x35.E_LUMO_eV.min())
    f['lumo_max'] = float(_x35.E_LUMO_eV.max())
    f['n_lumo_below4'] = int((_x35.E_LUMO_eV < -4).sum())

    # 공간별 정규화 거리와 포스겐의 최근접 경쟁자. 산문이 세 공간을 섞어
    # 적었던 자리이므로 공간을 키에 박아 둔다.
    _dt = pd.read_csv('outputs/distance_table_normalized.csv')
    _SP = {'A_rdkit': 'struct', 'B_electronic': 'elec', 'C_merged': 'merged'}
    for _cond, _g in _dt.groupby('condition'):
        _tag = _SP[_cond]
        _nn = float(_g.nn_median_space.iloc[0])
        for _, _r in _g.iterrows():
            _code = re.search(r'\((\w+)\)$', _r.real_agent)
            _code = _code.group(1) if _code else 'VX'
            f[f'{_tag}_{_code}'] = float(_r.d_sim_norm)
            if _code == 'CG':
                f[f'{_tag}_CG_rival'] = float(_r.d_interferent) / _nn

    # 2026-11-07 (10차 R5-B5). SI meth:conformer 가 관문의 범위를 '53 of the
    # 61 numeric columns' 과 'eight columns' 로 적는다. 손으로 적힌 수가
    # 아니라 관문 실행의 결과여야 한다.
    # 2026-11-07 (10차 R2/R3). 역치가 여섯째 자유도다. 본문과 SI 가 적는
    # 3 -> 5, panel35 0, 귀무기대 9.06/9.63 을 산출물에서 읽는다.
    # 2026-11-07 (10차 R3). 붕괴 문단의 여섯 수치를 collapse_sensitivity.py
    # 실측에서 읽는다. 종전 값은 어느 실행에서도 나오지 않았다.
    # 2026-11-07 (11차). README 재현 약속 6수치. 구판을 들고 있었다.
    _mc = {r[0]: r for r in json.load(open('outputs/multiplicity_c.json'))}
    f['rm_p28'] = _mc['panel28_designated3'][2]
    f['rm_p35'] = _mc['panel35_designated3'][2]
    f['rm_p35lit'] = _mc['panel35_designated_plus_literature8'][2]
    f['rm_q28'] = _mc['panel28_designated3'][4]
    f['rm_q35'] = _mc['panel35_designated3'][4]
    f['rm_q35lit'] = _mc['panel35_designated_plus_literature8'][4]

    # 2026-11-07. 커버레터의 분량 숫자. wordcount.py 실측이다.
    _la = json.load(open('outputs/length_audit.json'))
    f['len_body'] = _la['body_words']
    f['len_refs'] = _la['ref_words']
    f['len_pages'] = _la['acs_pages']
    f['len_figs'] = _la['n_figures']

    _cs = json.load(open('outputs/collapse_sensitivity.json'))
    f['coll_q_off'] = _cs['bh_min_q_uncollapsed']
    f['coll_q_on'] = _cs['bh_min_q_collapsed']
    _cm = {(c['dataset'], c['mode'], c['space']): c for c in _cs['cells']}
    f['coll_merged_off'] = _cm[('d2', 'uv', 'merged')]['rho_uncollapsed']
    f['coll_merged_on'] = _cm[('d2', 'uv', 'merged')]['rho_collapsed']
    f['coll_struct_off'] = _cm[('d2', 'uv', 'structural')]['rho_uncollapsed']
    f['coll_struct_on'] = _cm[('d2', 'uv', 'structural')]['rho_collapsed']

    # 2026-11-07 (11차 주입시험). 자유도 개수와 매칭 범위에 키가 없었다.
    # 'Five analyst choices' 와 'between three and six' 가 주입되어도
    # 통과했다. 개수는 본문이 열거하는 여섯이고, 범위 상한은 인코딩
    # 민감도의 최대 매칭 수(고유값<=5 변종)다.
    # 2026-11-07 (12차 A-2). 비교군이 일곱째다. 실험절이 적는 문자
    # 그대로의 판독(작용제 포함)은 감사된 적이 없었고, 재보니 이
    # 목록에서 가장 크게 움직인다 (28종 3->0, 문헌8 5->2).
    f['n_freedoms'] = 7
    f['match_range_hi'] = _V['drop_cardinality_le_5']['n_match'] \
        if 'drop_cardinality_le_5' in _V else 6

    _ts = json.load(open('outputs/stat_supplement.json'))['threshold_free']
    f['thr_free_28'] = _ts['panel28_incumbent3']['n_match_unbounded']
    f['thr_free_35'] = _ts['panel35_incumbent3']['n_match_unbounded']
    f['thr_null_28'] = _ts['panel28_incumbent3']['null_expectation_outside_closer']
    f['thr_null_35'] = _ts['panel35_incumbent3']['null_expectation_outside_closer']

    _gs = json.load(open('outputs/descriptor_gate_scope.json'))
    # 2026-11-07 (12차 A-4). 관문이 다섯 캐시를 돌므로 수치열 수가
    # 캐시마다 다르다 (53~62). 산문이 그 범위를 적으므로 하한을 건다.
    f['gate_numeric'] = _gs['n_numeric_min']
    f['gate_checked'] = _gs['n_checked']
    f['gate_exempt'] = len(_gs['exempt'])
    # 2026-11-07 (12차 A-4). 관문이 다섯 캐시를 돈다. 캐시 수와 수치열
    # 상한을 사실로 건다. 종전에는 refs[0] 하나만 기록됐고, 그래서
    # '각 캐시' 라는 서술을 어느 검사도 받지 않았다.
    f['gate_caches'] = _gs['n_caches']
    f['gate_numeric_max'] = _gs['n_numeric_max']

    # 염료 패널 분류 (2026-11-07 신설, 12차 A-9). SI meth:deltaE 가 적는
    # 비염료 개수와 채널 비를 산출물에 건다. 종전에는 '다섯' 이 손으로만
    # 적혀 있었고 어느 측정에도 걸려 있지 않았다.
    _dp = json.load(open('outputs/dye_panel_audit.json'))
    f['dye_n_components'] = _dp['n_components']
    f['dye_n_nondye'] = _dp['n_marked_nondye']
    f['dye_n_fluor'] = _dp['n_fluorophores']
    f['dye_top_ratio'] = _dp['highest_ratio']['ratio']
    f['dye_top_uv'] = _dp['highest_ratio']['uv_max']
    f['dye_top_normal'] = _dp['highest_ratio']['normal_max']
    _byname = {c['name']: c for c in _dp['components']}
    f['dye_anthracene_ratio'] = _byname['Anthracene']['ratio_max']
    f['dye_pyrene_ratio'] = _byname['Pyrene']['ratio_max']

    # 양수 가상궤도 (2026-11-07 신설, 12차 A-5).
    _pl = json.load(open('outputs/dft_electronic_space.json'))['positive_lumo']
    f['lumo_n'] = _pl['n_compounds']
    f['lumo_b3lyp_pos'] = _pl['n_b3lyp_positive']
    f['lumo_gfn2_pos'] = _pl['n_gfn2_positive']
    f['lumo_b3lyp_max'] = _pl['b3lyp_max']

    pc = json.load(open('outputs/dim_percentiles.json'))
    f['pct_struct'] = pc['struct_full']
    # 2026-10-27. 28종 구조공간 명목 p 가 어떤 어법으로도 등록돼 있지 않았다.
    # 그래서 'value of 0.029 sits at the 24.3rd' 같은 주입이 통과했다.
    f['p28_struct'] = float(json.load(open('outputs/permutation_statistics.json'))
                            ['panel28_designated']['A_rdkit']['statistics']
                            ['mean_diff']['p_one_sided'])
    f['pct_struct_caption'] = pc['struct_full']
    f['pct_elec'] = pc['electronic']
    f['pct_pc3'] = pc['struct_pc3']

    # 블록 가중 민감도 (2026-10-17 신설, SI res:blockbalance).
    bb = json.load(open('outputs/blockbalance.json'))
    _S6 = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
           'n_closer', 'n_valid']
    _BBTAG = {'panel35_incumbent3': '35',
              'panel35_incumbent_plus_lit8': '35lit',
              'panel28_incumbent3': '28'}
    for _lab, _tag in _BBTAG.items():
        _p = bb[_lab]['block_balanced']['p']
        f[f'bb_minp{_tag}'] = min(_p[s] for s in _S6)
        # 2026-11-07 (12차 A-3). 서술자 균등 쪽 최소값은 사실로 등록돼
        # 있지 않았고, 그래서 SI 산문의 0.151 이 어느 산출물과도 대조되지
        # 않았다. 실측은 0.192 다. 두 값을 모두 등록한다.
        f[f'bb_minpdesc{_tag}'] = min(
            bb[_lab]['equal_per_descriptor']['p'][s] for s in _S6)
        f[f'bb_match{_tag}'] = bb[_lab]['block_balanced']['n_valid']
        f[f'bb_cur_match{_tag}'] = bb[_lab]['equal_per_descriptor']['n_valid']
        f[f'bb_meandiff{_tag}'] = _p['mean_diff']
    f['bb_nstruct'] = bb['panel28_incumbent3']['n_struct']
    f['bb_nelec'] = bb['panel28_incumbent3']['n_elec']

    # 블록균등 병합공간 + 기존 구조·전자 = 18검정 가족의 BH 최소 q.
    def _bh_min(ps):
        o = sorted(ps)
        m = len(o)
        q = [o[i] * m / (i + 1) for i in range(m)]
        for i in range(m - 2, -1, -1):
            q[i] = min(q[i], q[i + 1])
        return q[0]
    _rc = json.load(open('outputs/rollback_c.json'))
    _RCMAP = {'panel28_designated3': ('panel28_incumbent3', '28'),
              'panel35_designated3': ('panel35_incumbent3', '35'),
              'panel35_designated_plus_literature8':
                  ('panel35_incumbent_plus_lit8', '35lit')}
    for _lab, (_bblab, _tag) in _RCMAP.items():
        _sp = _rc['panels'][_lab]
        _ps = [_sp[s]['statistics'][k]['p']
               for s in ('structural', 'electronic3') for k in _S6]
        _ps += [bb[_bblab]['block_balanced']['p'][k] for k in _S6]
        f[f'bbq_{_tag}'] = _bh_min(_ps)
        # 사전지정 통계량 3검정도 가중에 따라 달라질 수 있으므로 같이 둔다.
        f[f'bbpre_{_tag}'] = _bh_min(
            [_sp['structural']['statistics']['mean_diff']['p'],
             _sp['electronic3']['statistics']['mean_diff']['p'],
             bb[_bblab]['block_balanced']['p']['mean_diff']])

    # 네 번째 구성 — 농약 2종을 모사체로 분류 (2026-10-20 신설).
    # 본문 §2.1 이 '쟁점의 여지가 있다'고 적고 돌리지 않았던 라벨 결정을
    # 실제로 돌린 결과다. 이 논문에서 다중검정 보정을 통과하는 유일한 구성이라
    # 수치 하나가 틀리면 주장이 통째로 바뀐다. 표와 산문을 모두 등록한다.
    pp = json.load(open('outputs/pesticide_panel.json'))
    f['pp_min_p'] = pp['min_p']
    f['pp_bh_q'] = pp['bh_min_q']
    f['pp_n_perm'] = pp['spaces']['merged']['n_perm']
    f['pp_n_sim'] = pp['n_simulants']
    for _sp, _tag in [('structural', 'st'), ('electronic3', 'el'),
                      ('merged', 'mg')]:
        _d = pp['spaces'][_sp]
        for _k, _short in [('mean_diff', 'md'), ('median_diff', 'mdn'),
                           ('mean_logratio', 'lr'), ('mean_dsim', 'ds'),
                           ('n_closer', 'nc'), ('n_valid', 'nv')]:
            f[f'pp_{_tag}_{_short}'] = _d['p'][_k]
        f[f'pp_{_tag}_match'] = _d['n_valid']

    # 외부 문헌에서 옮겨 적은 수치 (2026-10-20). 계산값이 아니라 전사값이므로
    # 어떤 산출물과도 대조되지 않는 상태에 있었다. 표 3 사고와 같은 경로다.
    cf = json.load(open('outputs/cite_facts.json'))['papers']
    f['j18_identified'] = cf['james2018']['n_simulants_identified']
    f['j18_unsuitable'] = cf['james2018']['n_unsuitable']
    f['j21_identified'] = cf['james2021']['n_chemicals_identified']
    f['j21_discounted'] = cf['james2021']['n_discounted']
    f['j21_suitable'] = cf['james2021']['n_suitable']
    f['j21_possible'] = cf['james2021']['n_possible']

    # 네 번째 구성의 다중검정 체계 5종 (2026-10-20). 표 2 네 번째 열.
    pf = json.load(open('outputs/pesticide_family.json'))
    f['pf_min_p'] = pf['min_p']
    f['pf_bh18'] = pf['bh_all_min_q']
    f['pf_bh6'] = pf['bh_by_space_min_q']
    f['pf_wy'] = pf['westfall_young_adjusted_p']
    f['pf_meff'] = pf['meff']
    f['pf_bonf_meff'] = pf['meff_adjusted_min_p']
    f['pf_prespec'] = pf['prespecified_bh_min_q']

    # 서술자 감소 정산 (2026-10-20). 본문이 51 과 38 을 따로 적고 잇지 않았다.
    # 2026-10-25. 종전에는 13 이 리터럴이었다. 재분석으로 28종이 36좌표가 된
    # 뒤에도 그 리터럴이 남아 본문의 '두 패널에서 같다'는 서술을 받쳐 주었다.
    # 두 패널을 따로 산출물에서 읽는다.
    # 2026-10-26. 종전에는 35종 값을 experiment_summary_p35all.json 에서
    # 읽었다. 그 파일은 파이프라인이 만들지 않는 고아였고 전자공간에 omega 가
    # 아직 들어 있었다. 두 패널 모두 rollback_c.json 의 panels 에서 읽는다.
    _c28 = json.load(open('outputs/experiment_summary.json'))['conditions'][0]
    _pan = rb['panels']
    f['struct_candidates'] = int(_c28['n_in'])
    f['struct_kept_35'] = int(_pan['panel35_designated3']['structural']['n_descriptors'])
    f['struct_kept_28'] = int(_pan['panel28_designated3']['structural']['n_descriptors'])
    f['struct_removed'] = f['struct_candidates'] - f['struct_kept_35']
    f['struct_removed_28'] = f['struct_candidates'] - f['struct_kept_28']
    if f['struct_removed'] == f['struct_removed_28']:
        raise SystemExit('[verify] 중단: 두 패널의 서술자 감소가 같아졌다. '
                         '본문 서술과 검사 설계를 함께 손봐야 한다.')

    mp = json.load(open('outputs/multiplicity_c.json'))
    for lab, _m, _minp, _bonf, bhq, _k in mp:
        if lab == 'panel28_designated3':
            f['bh_q28'] = bhq
            f['minp28'] = _minp
        if lab == 'panel35_designated3':
            f['bh_q35'] = bhq
            f['minp35'] = _minp
        if lab == 'panel35_designated_plus_literature8':
            f['bh_q35lit'] = bhq
            f['minp35lit'] = _minp
    # 2026-11-03 (2차). 별칭은 사실표 끝에서 만든다. 앞에 두었더니
    # 뒤에 정의되는 사실 3건이 아직 없어 앵커가 조용히 꺼졌다.
    for _a, _b in [('rand3_frac28', 'rand3_frac28_2'),
                   ('rand3_med28', 'rand3_med28_2'),
                   ('white_med28', 'white_med28_2'),
                   ('pct_struct', 'pct_struct_2'),
                   ('fam_bh18_p28', 'fam_bh18_p28_2'),
                   ('mantel_d2nor_elec', 'mantel_d2nor_elec_2'),
                   ('struct_kept_28', 'struct_kept_28_2')]:
        if _a in f:
            f[_b] = f[_a]
    # ---------------------------------------------- 휘발성 축 (2026-11-03)
    # 보조 분석이지만 본문·SI 가 25개 수치를 인용한다. 앵커를 붙이지 않으면
    # 이 프로젝트가 반복해서 빠진 구멍(가드 없는 새 절)이 그대로 생긴다.
    _vol = json.load(open('outputs/volatility_space.json'))
    _vd = _vol['designation'][0]
    f['vol_p_md'] = _vd['statistics']['mean_diff']['p']
    f['vol_rank_md'] = _vd['statistics']['mean_diff']['rank']
    f['vol_obs_md'] = abs(_vd['statistics']['mean_diff']['observed'])
    f['vol_p_nvalid'] = _vd['statistics']['n_valid']['p']
    f['vol_rank_nvalid'] = _vd['statistics']['n_valid']['rank']
    f['vol_nperm'] = _vd['n_permutations']
    f['vol_range'] = _vd['sanity']['log10_range']
    f['vol_t1'] = _vd['sanity']['tiers']['T1']
    f['vol_tx'] = _vd['sanity']['tiers']['TX']
    # 2026-11-07 (10차 C3-2). 철회된 사전지정 통계량의 값과 4셀 단일가족
    # BH 를 원고에 실으면서 둘 다 묶는다. 10차 R3 가 'q=0.0430 이 사전지정
    # 실패 뒤 채택된 가족에 있는데 초록이 prespecified 라 부른다'를 잡았고,
    # 0.0809 는 그때까지 원고 어디에도 없었다.
    f['vol_withdrawn_p'] = _vol['agreement']['superseded_pooled_four_cells']['p_value']
    # BH 최소 q 를 순수 파이썬으로 센다. 여기서 numpy 를 지역 import 하면
    # 위쪽 중첩함수가 닫고 있는 _np 를 가려 NameError 가 난다 (2026-11-07).
    _vps = sorted(c['p_mantel'] for c in _vol['agreement']['cells'])
    _vn = len(_vps)
    _run, _vq = 1.0, []
    for _i in range(_vn - 1, -1, -1):
        _run = min(_run, _vps[_i] * _vn / (_i + 1))
        _vq.append(_run)
    f['vol_four_cell_bh'] = _vq[-1]
    _lw = _vol['lewisite_isomer_sensitivity']
    f['vol_lew_cis_p'] = _lw['cis']['p']
    f['vol_lew_cis_rank'] = _lw['cis']['rank']
    f['vol_t2'] = _vd['sanity']['tiers']['T2']
    f['vol_t3'] = _vd['sanity']['tiers']['T3']
    _bp = _vol['valid_match_detail']['sole_match_breakpoint']
    f['vol_gap'] = _bp['observed_gap_log10']
    f['vol_bp_ratio_gap'] = _bp['gap_breaking_ratio_condition_log10']
    f['vol_bp_ratio_vp'] = _bp['vp_breaking_ratio_condition_mmHg']
    f['vol_bp_closer_gap'] = _bp['gap_breaking_closer_condition_log10']
    f['vol_bp_closer_vp'] = _bp['vp_breaking_closer_condition_mmHg']
    # 2026-11-03 (작업 A·B). 형태좌표 축소 민감도, 전제검정 라벨순열 귀무,
    # MDE 의 BH18 문턱, 재표집 회복률. 본문에 실린 수치만 사실로 올린다.
    _sr = json.load(open('outputs/shape_reduction.json'))
    _srf = _sr['variants']['full']['panel28_incumbent3']
    _srr = _sr['variants']['reduced']['panel28_incumbent3']
    _srf35 = _sr['variants']['full']['panel35_incumbent3']
    f['sr_nshape'] = _srf['block_share_merged']['n_shape_kept']
    f['sr_dims28'] = _srf['n_dims']['merged']
    f['sr_dims35'] = _srf35['n_dims']['merged']
    f['sr_maxr28'] = _srf['shape_corr']['max_abs_r']
    f['sr_maxr35'] = _srf35['shape_corr']['max_abs_r']
    f['sr_share_full'] = _srf['block_share_merged']['shape_pct']
    f['sr_share_red'] = _srr['block_share_merged']['shape_pct']
    f['sr_nderived'] = len(_sr['shape_derived_dropped'])
    _acn_f = sorted(_srf['block_share_merged']['near_linear_pairs_shape_pct']
                    .values(), reverse=True)
    _acn_r = sorted(_srr['block_share_merged']['near_linear_pairs_shape_pct']
                    .values(), reverse=True)
    f['sr_acn_hi_full'], f['sr_acn_lo_full'] = _acn_f[0], _acn_f[1]
    f['sr_acn_hi_red'], f['sr_acn_lo_red'] = _acn_r[0], _acn_r[1]
    f['sr_match_full'] = -_srf['merged']['statistics']['n_valid']['observed']
    f['sr_match_red'] = -_srr['merged']['statistics']['n_valid']['observed']
    f['sr_mantel_q_red'] = _srr['mantel']['min_q']
    f['sr_id_lo'] = min(v for k, v in _srf['identity_audit'].items()
                        if k != 'note') * 100.0
    f['sr_id_hi'] = max(v for k, v in _srf['identity_audit'].items()
                        if k != 'note') * 100.0

    # 2026-11-04 (E2). 루이사이트 기하이성질체 민감도.
    _li = json.load(open('outputs/lewisite_isomer.json'))
    f['lew_dnorm_E'] = _li['E_trans']['d_norm']
    f['lew_dnorm_Z'] = _li['Z_cis']['d_norm']
    _pn = json.load(open('outputs/premise_null.json'))
    f['pn_cells'] = _pn['bh']['n_cells']
    f['pn_min_q'] = _pn['bh']['min_q']
    f['pn_n_sig'] = _pn['bh']['n_sig_q']
    # 2026-11-03 (작업 C). 바닥 논증을 대체한 문단이 쓰는 범위 전부.
    _pc = _pn['cells']
    f['pn_null_min'] = min(c['null_mean'] for c in _pc)
    f['pn_null_max'] = max(c['null_mean'] for c in _pc)
    f['pn_obs_min'] = min(c['observed'] for c in _pc)
    f['pn_obs_max'] = max(c['observed'] for c in _pc)
    f['pn_within_min'] = min(c['mean_within'] for c in _pc)
    f['pn_within_max'] = max(c['mean_within'] for c in _pc)
    f['pn_between_hi'] = max(c['mean_between'] for c in _pc)
    f['pn_between_lo'] = min(c['mean_between'] for c in _pc)
    f['pn_exact_n'] = max(c['n_null'] for c in _pc if c['enumerated'])
    _zs = [v for v in _pn['zero_sum_residual'].values() if 'error' not in v]
    # 2026-11-06 (9차 B1). 아래 min 이 'after_renorm_max' 를 읽고 있었다.
    # 셀별 최대값들의 최소를 구한 셈이라, 전역 최소(0.0445, d3_normal)가
    # 아니라 d3_normal 의 최대값(0.1292)이 나왔다. 원고가 그 0.1292 를
    # 하한으로 싣고 있었고 검증기가 같은 실수를 공유했으므로 461건 PASS 가
    # 이 양을 덮지 못했다. 필드를 바로잡는다 — 실제 하한은 2.9배 더 작고,
    # 즉 영합 전제가 깨졌다는 논거에 **불리한** 쪽이 아니라 유리한 쪽이다.
    f['pn_zs_min'] = min(v['after_renorm_min'] for v in _zs)
    f['pn_zs_max'] = max(v['after_renorm_max'] for v in _zs)
    f['pn_zs_before'] = max(v['before_renorm_max'] for v in _zs)

    # 2026-11-06 (9차 R3 6-4). 본문이 '모든 런-채널 셀에서 정확히 +1.000'
    # 이라고 적는데 예치물에는 두 셀뿐이었다(DMMP 없는 런을 건너뛴다).
    # deltaE_noise 에 여섯 셀 전수 산출을 신설하고 여기에 묶는다.
    _dn = json.load(open('outputs/deltaE_noise.json'))
    if 'concentration_vs_norm_all_cells' in _dn:
        f['conc_norm_cells'] = len(_dn['concentration_vs_norm_all_cells'])

    # 2026-11-07 (10차 C3-1). SI 가 37/25/15 를 적고 있었는데 예치물은
    # 34/24/16 이다. 9차에서 --n 40 을 호출부에 고정하면서 그 고정이
    # 맞춰야 할 본문 수치를 대조하지 않은 구멍이다. 세 수치를 묶는다.
    if _os.path.exists('outputs/order_invariance.json'):
        _oi = json.load(open('outputs/order_invariance.json'))
        _lg = _oi['cases']['legacy/panel28_incumbent3']
        f['oi_legacy28_sets'] = _lg['n_distinct_kept_sets']
        f['oi_legacy28_m3'] = _lg['n_match']['3']
        f['oi_legacy28_m4'] = _lg['n_match']['4']
        f['oi_legacy35_sets'] = (
            _oi['cases']['legacy/panel35_incumbent3']['n_distinct_kept_sets'])

    _md = json.load(open('outputs/mde.json'))
    f['mde_bh18_35'] = _md['panel35_incumbent3']['thresholds']['bh18']['mde_delta']
    f['mde_bh18_28'] = _md['panel28_incumbent3']['thresholds']['bh18']['mde_delta']
    f['mde_ratio_bh35'] = (_md['panel35_incumbent3']['thresholds']['bh18']
                           ['ratio_median_gap_to_mde'])
    f['mde_ratio_bh28'] = (_md['panel28_incumbent3']['thresholds']['bh18']
                           ['ratio_median_gap_to_mde'])
    f['mde_gap_median_28'] = _md['panel28_incumbent3']['gap_median']

    _nb = json.load(open('outputs/nearest_agent_bootstrap.json'))
    f['nb_nboot'] = _nb['n_boot']
    f['nb_dmmp'] = _nb['simulants']['DMMP']['top_pct']
    f['nb_dpgme'] = _nb['simulants']['DPGME']['top_pct']
    f['nb_dcm'] = _nb['simulants']['dichloromethane']['top_pct']

    _ag = _vol['agreement']
    # 2026-11-03 (작업 A-3). 4셀 일괄 Fisher 를 철회하고 채널내 런간 합성
    # + 채널 BH 로 교체했다. 같은 플레이트 두 채널은 독립이 아니다.
    _vp = _ag['prespecified']
    f['vol_fisher_uv'] = _vp['by_channel']['uv']['p_value']
    f['vol_fisher_nor'] = _vp['by_channel']['normal']['p_value']
    f['vol_fisher_minq'] = _vp['min_q']
    f['vol_channels'] = _vp['n_channels']
    for _c in _ag['cells']:
        f[f"vol_rho_{_c['cell']}"] = _c['spearman']
    # 2026-11-07 (11차 주입시험). 같은 수치가 본문과 SI 양쪽에 있을 때
    # 키 하나에 앵커 하나면 둘 중 먼저 걸리는 자리만 지켜지고 나머지는
    # 무방비가 된다. 주입시험이 본문 자리 두 곳에서 그것을 입증했다
    # (0.0430 -> 0.1430, 0.4647 -> 0.0647 가 통과했다). 본문 전용 키를 둔다.
    f['vol_fisher_minq_main'] = _vp['min_q']
    f['vol_p_nvalid_main'] = _vd['statistics']['n_valid']['p']

    _pt = {r['agent']: r for r in _vol['pairing_table']['rows']}
    f['vol_gap_cg'] = _pt['phosgene (CG)']['gap_to_nearest_simulant']
    f['vol_gap_ac'] = _pt['hydrogen cyanide (AC)']['gap_to_nearest_simulant']
    f['vol_gap_gb'] = _pt['sarin (GB)']['gap_to_nearest_simulant']
    f['vol_gap_gb_pool'] = _pt['sarin (GB)']['gap_to_nearest_pool']
    f['vol_gap_hd'] = _pt['sulfur mustard (HD)']['gap_to_nearest_simulant']
    f['vol_gap_vx'] = _pt['VX']['gap_to_nearest_simulant']
    f['vol_n_int_nearest'] = _vol['pairing_table'][
        'n_agents_whose_nearest_is_interferent']

    return f

def nums_in(s):
    return [float(x.replace(',', '').replace('\\,', ''))
            for x in re.findall(r'(?<![\d.])(\d[\d,]*(?:\\,\d{3})*(?:\.\d+)?)', s)]


# (사실표 키, 본문에서 그 양을 가리키는 문맥 정규식, 허용 오차)
# 문맥 정규식의 첫 캡처 그룹이 검사할 수치다.
# 2026-10-24 신설. 표 S1 의 칸 검사는 지금까지 열 간격을 글자 그대로 박아
# 두었다 ('& 38 &', 'structural  &'). 표를 산출물에서 다시 생성하자 간격이
# 달라져 12건이 한꺼번에 '문맥 미발견'으로 꺼졌다. 꺼진 검사는 통과처럼
# 보이므로 가장 위험하다. 행 이름·공간·열 번호로 패턴을 만든다.
def row_pat(label, space, col):
    """표 S1 의 (행, 공간, col 번째 통계) 칸을 잡는 정규식."""
    skip = r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*'
    cap = r'(?:\\textbf\{)?([\d.]+)\}?'
    return (re.escape(label) + r'\s*&\s*' + space
            + r'\s*&\s*\d+\s*&\s*' + skip * col + cap)


# 본문에 인용하지 않기로 한 사실. 키마다 이유를 적는다. 비워 두는 것이 정상이고,
# 여기 들어오는 항목은 '검사를 끈다'는 선언이므로 심사 대상이다.
UNQUOTED = {
    # 2026-11-03. 처음에는 pp_el_* / pp_mg_* 여섯 건도 여기 올렸다. 그것은
    # 내 오판이었다. 행이 삭제된 것이 아니라, make_si_tables 의 ROWS 가
    # 블록마다 이름 두 개만 두고 농약 구성의 세 줄짜리 이름을 한 칸에
    # 넣었던 탓에 이름이 공간과 한 칸씩 어긋나 찍혀 앵커가 빗나간 것이다.
    # 생성기의 배치를 손표대로 되돌려 여섯 건 모두 되살렸다. 검사를 끄는
    # 쪽으로 기울면 이런 오판이 통과한다는 사례로 적어 둔다.
    #
    # C(24,10) 전수열거 횟수는 농약 구성 절을 줄이면서 본문에서 빠졌다.
    # 표본이 아니라 전수열거라는 사실 자체는 2.1 이 C(24,3) 으로 이미
    # 말하고 있다. 이 한 건만 면제로 남긴다.
    'pp_n_perm': '농약 구성의 전수열거 횟수를 본문에 싣지 않기로 함',

    # ------------------------------------------------------------------
    # 2026-11-04 (결정 J-C). JHM 평행본을 archive/jhm_v1/ 로 동결하고
    # 1차 본문을 paper_acs/ 로 옮겼다. 그러자 아래 사실들이 '문맥
    # 미발견' 으로 떠올랐다. 전부 JHM 본문·SI 에만 적혀 있던 수치이고,
    # ACS본은 분량을 줄이면서 그 서술을 싣지 않는다. 주석을 벗기고
    # 값별로 대조해 확인했다.
    #
    # 중요한 점: 이 항목들이 틀렸다는 뜻이 아니다. 투고되는 문서가
    # 그 수치를 말하지 않으므로 대조할 문장이 없다는 뜻이다. ACS본이
    # 실어야 할 서술이 이 안에 있는지는 편집 판단이고, 검증기가
    # 답할 문제가 아니다. 동결본을 되살릴 때(ACS 데스크리젝 시)
    # 이 블록을 통째로 되돌리면 그쪽 검사가 다시 살아난다.
    # ------------------------------------------------------------------
    'pc1_struct28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pc12_elec28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pc1_struct35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pc12_elec35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pc12_merged28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pc12_merged35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'lumo_min': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'lumo_max': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'n_lumo_below4': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'nullexp_k3': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'nullexp_k8': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'na_dmmp_top': 'ACS본 미수록 (J-C). nb_dmmp 가 같은 양을 ACS 어법으로 검사한다',
    'na_dmmp_g': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'na_dcm_g': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'na_dpgme_top': 'ACS본 미수록 (J-C). nb_dpgme 가 같은 양을 ACS 어법으로 검사한다',
    'na_dpgme_g': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'fig3_from_panel': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'fig3_unlabelled': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'r_lumo_omega': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'omega_distrank': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'omega_pc3_var': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'nopatho_struct_p': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'nopatho_elec_p': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'nopatho_merged_p': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'null3_lower': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'null8_nested_mean': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'null8_nested_p': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'match35_lit': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'match28_inc': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dfp_gb': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dfp_gd': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dfp_gf': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dmmp_gb_35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dmmp_gd_35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'dmmp_gf_35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'cees_hd': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'mde_gap_min': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'null_exp35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'null_exp28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'boot_diff_pos': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'ret_gb3': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'ret_gd3': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'ret_gf3': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'ret_ga8': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_GD': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_CK': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_GF': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_AC': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_CG': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_HD': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_GA': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_L': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_PS': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    't3_VX': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_GD': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_GF': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_VX': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'merged_VX': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'elec_AC': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_CG': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'elec_CG': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_CG_rival': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'merged_CG_rival': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'elec_CG_rival': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pp_min_p': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_removed': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_removed_28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'struct_kept_28': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    # 순도 민감도의 '38 to 37' 은 28종 패널의 병합 차원이고 이 사실의
    # 35종 구조좌표 수와 값만 같다. 앵커로 쓰려다 검증기가 41 을 잡아
    # 냈다. 우연한 동일값에 앵커를 걸지 않는다.
    'struct_kept_35': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'pp_bh_q': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j18_identified': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j18_unsuitable': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j21_identified': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j21_discounted': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j21_suitable': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
    'j21_possible': 'ACS본 미수록 (J-C) — JHM 동결본 전용 수치',
}


# 2026-11-02 신설. 문헌 인용 수치.
#
# load_facts() 는 산출물에서만 읽는다는 원칙을 지킨다. 그런데 ACS 2.1 은
# 전제를 우리 계산이 아니라 세 원논문 자신의 분류 성능으로 세우므로, 그
# 수치는 우리 산출물에 존재하지 않는다. 손으로 적되, 근거를 함께 적는다.
#
# 규칙: 값마다 (1) 어느 논문인지 (2) 그 논문의 원문 한 문장을 함께 둔다.
# 원문이 없는 값은 여기 들어올 수 없다. 원문은 2026-11-02 에 Ku 가 올린
# 세 PDF 에서 직접 추출했다.
LITERATURE = {
    # d1 = hazref, J. Hazard. Mater. 501 (2026) 140703
    'lit_d1_acc': (100.0,
                   'hazref',
                   'our linear discriminant analysis classifier achieved '
                   '100% classification accuracy'),
    'lit_d1_vis': (87.5,
                   'hazref',
                   'elevating the classification performance from 87.5% to 100%'),
    # d2 = multiplexref, Discover Sensors 1(1) (2025) 29
    'lit_d2_acc': (0.964,
                   'multiplexref',
                   'The overall classification accuracy was 0.964, with a '
                   'macro F1-score of 0.964'),
    # d3 = group3ref, Sens. Actuators B 457 (2026) 139706
    'lit_d3_acc': (1.000,
                   'group3ref',
                   'demonstrated perfect accuracy (accuracy = 1.000, '
                   'macro F1 = 1.000) across all classifications'),
    'lit_d3_lod': (10.0,
                   'group3ref',
                   'The analytical working range was 10-500 uM with an '
                   'calculated limit of detection of ca. 10 uM'),
}

CHECKS = [
    # 2026-11-02. ACS 2.1 전제를 떠받치는 문헌 수치. 값의 출처는 LITERATURE.
    ('lit_d1_acc', [r'dimethyl methylphosphonate at (\d+)\\%\s*\n?accuracy'], 0),
    ('lit_d1_vis', [r'against (\d+\.\d)\\% under visible illumination alone'], 0.001),
    ('lit_d2_acc', [r'all seven agents at an accuracy of\s*\n?(\d\.\d{3})'], 0.0005),
    ('lit_d3_acc', [r'four Novichok agents at an accuracy of (\d\.\d{3})'], 0.0005),
    ('lit_d3_lod', [r'with a detection limit near (\d+)~\$\\mu\$M'], 0),
    ('n_train', [r'random forest trained on (\d+) AEGL chemicals',
     r'(\d[\d,]*) training (?:chemicals|compounds)',
                 r'leaving (\d[\d,]*) for\s*\n?training',
                 r'trained on (\d[\d,]*) AEGL'], 0),
    ('n_heldout', [r'(?:the|all) (twenty|seventeen|nineteen|eighteen) '
                   r'(?:military )?agents held out',
                   r'On the (twenty|seventeen) agents held out'], 0),
    ('n_candidates', [r'(?:the|Across the|across the) (\d[\d,]*) candidates',
                      r'excluding the authentic agents\s*\n?themselves \(\$n=(\d+)\$\)',
                      r'(\d[\d,]*)-compound\s*\n?AEGL pool',
                      r'all but \w+ of the (\d[\d,]*) compounds in the',
                      r'association across the (\d[\d,]*)\s*\n?compounds',
                      r'for the (\d[\d,]*) candidates\s*\n?carrying both axes',
                      r'of the (\d[\d,]*), none combines'], 0),
    # 2026-11-04 (J-C). SI 의 AEGL-3 순서 문장. 다섯 값 전부를 가린다.
    ('aegl3_vx', [r'VX \((0\.\d+)~mg~m\$\^\{-3\}\$\),'], 0.00005),
    ('aegl3_gd', [r'soman \((0\.\d+)\), sarin'], 0.0005),
    ('aegl3_gb', [r'sarin \((0\.\d+)\), cyclosarin'], 0.0005),
    ('aegl3_gf', [r'cyclosarin \((0\.\d+)\), tabun'], 0.0005),
    ('aegl3_ga', [r'tabun \((0\.\d+)\), sulfur'], 0.0005),
    ('dcm_aegl3', [r'dichloromethane[,\s]+(\d[\d,\\]*)~mg~m',
                   r'at\s*\n?(\d[\d,\\]*)~mg~m\$\^\{-3\}\$ it sits at the',
                   r"dichloromethane's (\d[\d,\\]*)~mg~m"], 1.0),
    ('aegl3_chloroform', [r'chloroform \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_toluene', [r'toluene \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_methanol', [r'methanol \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_hexane', [r'n-hexane\s*\n?\((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_acetone', [r'acetone \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_benzene', [r'benzene \((\d[\d,\\]*)\)'], 1.0),
    ('aegl3_pxylene', [r'p-xylene \((\d[\d,\\]*)\)'], 1.0),
    ('rf_r2', [r'reaches \$R\^\{2\}=([\d.]+)\$ in cross-validation',
     r'R\^\{2\}=([\d.]+) \\pm'], 0.001),
    ('rf_r2_sd', [r'R\^\{2\}=[\d.]+ \\pm ([\d.]+)'], 0.001),
    ('rf_rmse', [r'RMSE \$?([\d.]+)\$? in'], 0.001),
    ('b1_rmse', [r'training-mean baseline of\s*\n?([\d.]+) \('], 0.001),
    ('b2_rmse', [r'linear baseline of ([\d.]+)\.'], 0.001),
    ('dmmp_ad', [r"DMMP's domain distance, ([\d.]+),"], 0.005),
    ('n_ache', [r'(Thirteen|Eleven|Twelve) of the\s*\n?\d+ training compounds carry'], 0),
    ('n_relaxed_10', [r'(Eight|Seven|Nine) compounds then combine',
                      r'there are (eight|seven|nine) candidates'], 0),
    ('n_relaxed_100', [r'the count rises to (fifteen|seventeen|sixteen)'], 0),
    ('n_incumbent', [r'Faithful and safe at that bar & (\d+) '], 0),
    ('pc1_struct28', [r'size axis carrying ([\d.]+)\\% of the variance'], 0.06),
    ('pc12_elec28', [r'the first two\s*\n?carrying ([\d.]+)\\%'], 0.06),
    ('pc1_struct35', [r'The primary panel behaves the same way, at ([\d.]+)\\%'], 0.06),
    ('pc12_elec35', [r'the same way, at [\d.]+\\% and ([\d.]+)\\%'], 0.06),
    ('pc12_merged28', [r'([\d.]+)\\% on the smaller panel'], 0.06),
    ('pc12_merged35', [r'and ([\d.]+)\\% on the primary one'], 0.06),
    ('lumo_min', [r'they run from \$([-+][\d.]+)\$~eV for chloropicrin'], 0.005),
    ('lumo_max', [r'to \$([-+][\d.]+)\$~eV for water'], 0.005),
    ('n_lumo_below4', [r'and (twenty-four|twenty-three|twenty-five) of the\s*\n?'
                       r'thirty-five fall below'], 0),
    ('n_more_faithful_10', [r'more faithful than dichloromethane itself & (\d+) ',
                            r'and all (eight|seven|nine) are more faithful than'], 0),
    ('fid_cut', [r'own fidelity of\s*\n?0\.78 falling outside the ([\d.]+) quartile',
                 r'methanol at 2\.00 against the ([\d.]+) boundary'], 0.005),
    ('n_screened_distinct', [r'(Eight|Seven|Nine) further records were set aside',
                             r'Retaining all (eight|seven|nine)\s*\n?records'], 0),
    ('sens_kept_n', [r'raises the candidate count to (\d+) and'], 0),
    ('sens_kept_rho', [r'moves the Spearman coefficient\s*\n?from \$\+0\.062\$ to '
                       r'\$([-+][\d.]+)\$'], 0.0005),
    ('panel_with_aegl3', [r'(Nineteen|Thirteen|Fourteen|Eighteen) of the '
                          r'twenty-eight compounds have an AEGL-3',
                          r'the (nineteen|thirteen|fourteen) panel compounds',
                          r'(nineteen|thirteen|fourteen) compounds that have one'], 0),
    ('xreader_pairs', [r'agree on all (\d[\d,]*) entries for which both'], 0),
    ('trainset_checked', [r'on all (\d[\d,]*) of\s*\n?those that resolved to a structure'], 0),
    ('comma_rows', [r'(\w+[-\w]*) values reach 1000 or more'], 0),
    ('panel_from_footnote',
     [r'(Five|Nine|Four) of the nineteen values come from the',
      r'(Five|Nine|Four) of the nineteen panel values'], 0),
    ('nullexp_k3',
     [r'designating three\s*\n?non-agents at random matches (\d+\.\d+) agents',
      r'random matches 2\.68 agents on average against (\d+\.\d+) for three'], 0.05),
    ('nullexp_k8', [r'designating eight non-agents at\s*\n?random matches (\d+\.\d+) agents'], 0.05),
    ('na_dmmp_top', [r'DMMP is nearest sarin in\s*\n?(\d+\.\d+)\\% of replicates'], 0.05),
    ('na_dmmp_g',
     [r'nearest one of the four G-series agents in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dcm_g',
     [r'Dichloromethane is nearest a G-series agent in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dpgme_top',
     [r'DPGME divides evenly, sulfur mustard in\s*\n?(\d+\.\d+)\\%'], 0.05),
    ('na_dpgme_g', [r'and the\s*\n?G series in (\d+\.\d+)\\%'], 0.05),
    # 2026-11-07 (12차 Step 5). TOC 패널 (a) 가 DMMP 코사인에서 역전 계수로
    # 바뀌었으므로 toc_cos_* 다섯 사실을 폐기하고 역전 두 사실로 건다.
    # 폐기한 사실은 facts 쪽에서도 지웠다 (등재부가 더 이상 그 글자를 적지
    # 않으므로 앵커만 남기면 영구 미발견이 된다).
    ('toc_rev_uncentred',
     [r'different compound in (\d+) of 34 cases before centring'], 0),
    ('toc_rev_centred',
     [r'before centring and (\d+) of 34 after'], 0),
    ('toc_rho_min', [r'12 cells: ([-+][\d.]+) to'], 0.0005),
    ('toc_rho_max', [r'12 cells: [-+][\d.]+ to ([-+][\d.]+)'], 0.0005),
    ('toc_cells', [r'rank correlation, (\d+) cells'], 0),
    ('fig3_from_panel', [r"name the (ten|nine|eleven) agents of this paper's panel"], 0),
    ('fig3_other_labelled',
     [r'the\s*\n?(three|two|four) nitrogen mustards'], 0),
    ('fig3_unlabelled',
     [r'the remaining (seven|six|eight) held-out agents are plotted'], 0),
    ('panel_without', [r'\(b\) The (nine|fifteen|ten)\s+compounds with no\s+regulator-grade',
                       r'(Nine|Fifteen|Ten) compounds have no value'], 0),
    ('heldout_worst_in_domain',
     [r'worst of those fourteen is wrong by ([\d.]+) log',
     r'worst of those fourteen being\s*\n?wrong by ([\d.]+) log units',
      r'is still wrong\s*\n?by ([\d.]+) log units'], 0.005),
    ('heldout_in_domain',
     [r'certifies (\w+) of the twenty\s+as inside the domain',
     r'(Fourteen|Thirteen|Fifteen) of the twenty lie inside',
      r'(fourteen|thirteen|fifteen) of them pass the domain check'], 0),
    ('tce_ratio', [r'trichloroethylene[, ]+is ([\d.]+) times closer',
                   r'trichloroethylene ([\d.]+) times closer'], 0.05),
    ('dcm_above', [r'with (seven|eight|two) compounds\s*\n?above it',
                   r'all but (seven|eight|two) of the 214 compounds'], 0),
    ('panel_without', [r'inhalation\s*\n?\s*value exists for any of these (\d+)'], 0),
    ('dcm_panel_rank', [r'no designating\s*\n?document, ranks (\d+) of \d+'], 0),
    ('panel_with_aegl3', [r'no designating\s*\n?document, ranks \d+ of (\d+)'], 0),
    ('heldout_in_domain', [r'(\d+) of the 20 agents pass the domain check'], 0),
    ('heldout_n', [r'Applied to (\w+) military agents held out',
     r'and all (\d+) are wrong by'], 0),
    ('heldout_worst_in_domain',
     [r'worst inside the domain is wrong by ([\d.]+) log'], 0.005),
    ('spearman_r', [r'Spearman \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?',
                    r'Spearman \$r=([-+][\d.]+)\$'], 0.0006),
    ('spearman_lo', [r'Spearman.{0,60}?95\\% CI \$?([-+][\d.]+)\$? to'], 0.001),
    ('spearman_hi', [r'Spearman.{0,60}?95\\% CI .{0,20}?to\s*\n?\$?([-+][\d.]+)\$?'], 0.001),
    ('pearson_r', [r'Pearson \$?r\$? ?[&=] ?\$?([-+][\d.]+)\$?'], 0.0006),
    # 2026-10-13: 본문이 재집필되어 종전 정규식의 문맥이 사라졌다. kappa 민감도
    # 문장과 Table 1 의 '28종 패널' 블록을 동시에 대조한다. 표에는 패널 3구성이
    # 들어 있으므로 '28 cpd.' 행에 앵커를 걸어야 한다. 앵커 없이 쓰면 35종 행의
    # 값까지 물어와 거짓 불일치가 난다.
    ('p_struct', [r'with the structural\s+space at ([\d.]+)',
     r'condition \(structural \$p=([\d.]+)\$',
                  row_pat('28 cpd., incumbent 3', 'structural', 0)], 0.0005),
    # 2026-11-03. 'electronic $p=...$' 어법은 본문에서 사라졌고, 삭제판
    # 행 이름도 (omits the five) -> (omits the seven) 로 바뀌어 두 패턴이
    # 함께 꺼져 있었다. 살아 있는 두 자리로 다시 가리킨다.
    ('p_elec', [r'one-sided \$p\$ of ([\d.]+) to five matches',
                r'a one-sided \$p\$ of ([\d.]+)\.',
                row_pat('(omits the seven)', 'electronic', 0)], 0.0005),
    # 각 블록의 세 번째 행은 행 이름 칸이 비어 있다. 앞 행을 닻으로 쓴다.
    ('p_merged', [r'merged space returns a nominal one-sided \$p\$ of\s+([\d.]+)',
     r'merged \$p=([\d.]+)\$',
                  r'\(omits the five\)[^\n]*\n\s*&\s*merged\s*&\s*\d+\s*&\s*'
                  r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0005),
    ('pc3_p28', [r'from \$p=[\d.]+\$ to\s*\n?\$p=([\d.]+)\$'], 0.0005),
    ('rand3_frac28', [r'and ([\d.]+)\\% of the subspaces reach'], 0.05),
    ('rand3_med28', [r'panel the median \$p\$ is ([\d.]+)'], 0.0005),
    ('r_lumo_omega', [r'r\(E_\{\\mathrm\{LUMO\}\}, \\omega\) =\s*\n?([-+][\d.]+)'],
     0.0005),
    ('omega_distrank', [r'0\.934 on the primary panel, ([\d.]+) on the smaller one'], 0.0005),
    ('omega_pc3_var', [r'account for ([\d.]+)\\% of its variance'], 0.05),
    ('bh_q28', [r'over the eighteen tests gives ([\d.]+),\s+correcting within',
     r'([\d.]+) for Benjamini--Hochberg over the\s*\n?eighteen tests',
     r'eighteen-test family\s*\n?gives a smallest adjusted value '
                r'of ([\d.]+), correcting',
                r'Benjamini-Hochberg \$q =\s*\n?([\d.]+)\$'], 0.0005),
    ('bh_q35', [r'Benjamini--Hochberg adjusted value is\s+([\d.]+)',
     r'Benjamini-Hochberg over 18\s+& (\d+\.\d+) &'], 0.0005),
    ('bh_q35lit', [r'Benjamini-Hochberg over 18\s+& [\d.]+ & (\d+\.\d+) &'], 0.0005),
    # ACS 판본 Conclusions 가 2.4 의 같은 값을 다시 센다. 자리별 전용 앵커.
    ('bh_q35_concl',
     [r'at a Benjamini--Hochberg adjusted value of\s+(\d+\.\d+)'], 0.0005),
    # ACS 커버레터 전용 앵커 (2026-11-02). 본문과 같은 산출물에 건다.
    # 2026-11-03 (T1-10). 표지서신에서 분류 정확도 주장을 뺐다. 그 수치는
    # 본문 2.1 에 남아 있고 lit_d2_acc / lit_d3_acc 가 검사한다. 자리별
    # 별칭을 은퇴시킨다. 표지서신이 그 주장을 다시 실으면 아래를 되살릴 것.
    #   ('cl_d2_acc', [r'cross-validated accuracy of (\d\.\d{3}) on'], 0.0005),
    #   ('cl_d3_acc', [r'seven-agent run and (\d\.\d{3}) on the Novichok'], 0.0005),
    ('cl_dmmp_ga', [r'correlates at\s*\n?\$\+(\d\.\d{2})\$ with tabun'], 0.005),
    # 2026-11-07 (12차 Step 5). 커버레터 문장에 'in the ultraviolet
    # channel' 이 끼었다.
    ('cl_dmmp_gb', [r'\$\+(\d\.\d{2})\$ with\s*\n?sarin[^.]{0,45}?but at'],
     0.005),
    ('cl_dmmp_gd', [r'but at\s*\n?\$(-\d\.\d{2})\$ with soman'], 0.005),
    ('cl_dmmp_gf', [r'\$(-\d\.\d{2})\$ with cyclosarin'], 0.005),
    ('cl_dmmp_vx', [r'\$(-\d\.\d{2})\$ with VX'], 0.005),
    # 2026-11-03 (T1-9). 표지서신에서 농도 추세 주장을 뺐다. 정확 p 가
    # 0.083 밑으로 갈 수 없고 SNR 교란을 분리할 수 없으므로 표지서신에
    # 실을 수 없다. 본문 쪽은 conc_gb_lo / conc_gb_hi 가 검사한다.
    #   ('cl_conc_hi', [...], 0.005),
    #   ('cl_conc_lo', [...], 0.005),
    # 2026-11-01 신설. ACS 판본 Results 2.1~2.3 의 수치.
    ('dmmp_gd_uv', [r'anti-correlates with soman at\s+\$([-][\d.]+)\$'], 0.0015),
    # 2026-11-07 (12차 Step 5). 종전 정규식은 'cyclosarin at $-x$' 를 문맥
    # 없이 잡았고, 신설한 기저선 문장의 '사린 대 사이클로사린' 값(-0.820)
    # 에도 걸렸다. 그 쌍은 DMMP 와 무관하다. 'against' 앞머리를 배제한다.
    ('dmmp_gf_uv', [r'(?<!against )cyclosarin\s+at \$([-][\d.]+)\$'], 0.0015),
    ('dmmp_vx_uv', [r'and VX at \$([-][\d.]+)\$'], 0.0015),
    # 2026-11-03 (2차) 폐기. 'negative in both channels' 문장이 분모
    # 오류여서 교체됐다. 같은 산출물 값(3종 평균)은 dmmp_gmean_uv_3 /
    # dmmp_gmean_nor_3 가 새 어법 자리에서 검사한다.
    # 2026-11-06 (9차 R3-3). 문장을 재집필했다. 재표집 단위가 설계농도
    # 4개라는 것과 G 평균이 타분을 뺀 3종이라는 것을 본문에 적고, 구간을
    # 검정이 아니라 서술 범위로 읽도록 바꿨다. 앵커도 함께 옮긴다.
    ('gvx_uv', [r'excluding tabun, is \$\+([\d.]+)\$ in the\s*\n?ultraviolet channel'], 0.0015),
    ('gvx_uv_lo', [r'with them, \$\+([\d.]+)\$ to'], 0.0015),
    ('gvx_uv_hi', [r'with them, \$\+[\d.]+\$ to \$\+([\d.]+)\$ and'], 0.0015),
    ('gvx_nor', [r'ultraviolet channel and \$\+([\d.]+)\$ in the reflectance',
                 r'ultraviolet\s*\n?channel and \$\+([\d.]+)\$ in reflectance'], 0.0015),
    ('gvx_nor_lo', [r'and \$([-][\d.]+)\$ to \$\+[\d.]+\$, come from'], 0.0015),
    ('gvx_nor_hi',
     [r'and \$[-][\d.]+\$ to \$\+([\d.]+)\$, come from'], 0.0015),
    # 2026-11-02. 원논문 Methods 를 직접 읽어 농도 단위를 ppm 에서 µM 으로
    # 바로잡았다. 단위를 앵커에서 빼면 같은 실수를 다시 못 잡으므로 단위째
    # 교정하고 남겨 둔다.
    ('conc_gb_lo', [r'rises\s+monotonically from \$\+([\d.]+)\$ at 10~\$\\mu\$M'], 0.0015),
    ('conc_gb_hi', [r'to \$\+([\d.]+)\$ at 500~\$\\mu\$M'], 0.0015),
    ('conc_gf_lo', [r'falls monotonically from \$-([\d.]+)\$'], 0.0015),
    ('conc_gf_hi', [r'falls monotonically from \$-[\d.]+\$ to \$-([\d.]+)\$'], 0.0015),
    ('trend_monotone',
     [r'Monotonicity is perfect in (\w+) of the ten',
      r'monotonicity perfect in\s*\n?(\w+) of the ten'], 0),
    ('trend_bh_q',
     [r'over the ten leaves a smallest adjusted\s+value of ([\d.]+)'], 0.0015),
    # 2026-11-06 (9차 분량). 세 앵커가 접속사까지 글자로 박고 있어 문장을
    # 줄이자 끊겼다. 'against four' 와 같은 취약성이므로 수치 주변만 남긴다.
    ('rank_obs', [r'the observed number is (\w+)[,.]'], 0),
    ('rank_exp',
     [r'over ten\s+agent-channel cells\s+is ([\d.]+)[,.]'], 0.0015),
    ('rank_p', [r'one-sided binomial \$p\$ (?:of |is )([\d.]+)'], 0.005),
    ('dye_null_uv', [r'gives a mean of ([\d.]+) and [\d.]+\.'], 0.0015),
    ('dye_null_normal', [r'gives a mean of [\d.]+ and ([\d.]+)\.'], 0.0015),
    ('fisher_minq',
     [r'smallest adjusted value over the six is ([\d.]+)'], 0.0015),
    ('mantel_cells', [r'Of (twelve|eight|ten) panel-channel-space cells',
                      r'(twelve|eight|ten) cells are one\s*\n?family'], 0),
    ('mantel_nominal', [r'(three|two|four|five) of the twelve panel-channel-space cells reach'], 0),
    # rho 를 패턴에 박아 두면 rho 가 바뀔 때마다 검사가 조용히 빠진다
    # (2026-10-24 에 두 번 겪었다). 셀을 가리키는 문구로 건다.
    ('mantel_d2uv_struct',
     [r'that same set and channel \(\$\\rho = [-+][\d.]+\$, Mantel\s*\n?'
      r'\$p = ([\d.]+)\$'], 0.0005),
    ('ap_ntests', [r'one family of (\d+) tests'], 0),
    ('ap_p35', [r'Benjamini-Hochberg over all 72\s*& (\d+\.\d+) &'], 0.0006),
    ('ap_p35lit', [r'Benjamini-Hochberg over all 72\s*& [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('ap_p28', [r'Benjamini-Hochberg over all 72\s*& [\d.]+ & [\d.]+ & (\d+\.\d+) &',
                r'panel rises from 0\.165 to (\d+\.\d+)'], 0.0006),
    ('ap_pest', [r'configuration stays significant at (\d+\.\d+)'], 0.0006),
    ('sn_mean', [r'where the mean is ([\d.]+)\\% and',
     r'averages (\d+\.\d+)\\% across the 680'], 0.06),
    ('sn_median', [r'and the median ([\d.]+)\\%',
     r'with a median of (\d+\.\d+)\\%'], 0.06),
    ('sn_p95', [r'95th percentile of (\d+\.\d+)\\%'], 0.06),
    ('sn_max', [r'a maximum\s*\n?of (\d+\.\d+)\\%'], 0.06),
    # 2026-11-04 (작업 D). 순위를 어법에 박아 두면 순위가 바뀔 때 앵커가
    # 죽는다. 실제로 37좌표 -> 35좌표 교정에서 5위가 7위로 바뀌며 죽었다.
    # 서수를 자리표시로 바꾼다.
    ('sn_p', [r'ranks [a-z]+ of \d+\s*\n?\(\$p=([\d.]+)\$\)',
     r'ranks [a-z]+ of 680, \$p=(\d+\.\d+)\$'], 0.0006),
    # 2026-11-04 (작업 D). 본문이 '0.0824' 를 싣고 있었고 어느 산출물에도
    # 없는 값이었다. 생성표 tab_family 는 같은 양을 0.077 로 적는다.
    # R2·R5 가 같은 자리를 지적했다. 앵커가 없어서 통과했다.
    ('q_by_space28', [r'correcting within each space gives (0\.\d+)',
                      r'tests, (0\.\d+) within each space'],
     0.00005),
    ('h1_safer7', [r'only ([a-z]+) of the \d+ candidate compounds are'], 0),
    ('h1_ncand', [r'only [a-z]+ of the (\d+) candidate compounds are',
                  r'Across the (\d+) candidates carrying both'], 0),
    ('h1_relaxed10', [r'there are ([a-z]+) compounds, all of'], 0),
    ('sn_nsub', [r'Across (\d+) random three-descriptor subspaces',
     r'Over (\d+) subspaces'], 0),
    # 2026-11-01 신설. 사실로는 적재되는데 검사 항목이 없던 셋.
    ('sn_ndesig', [r'ranks [a-z]+ of (\d+)'], 0),
    ('sn_rank', [r'designation ranks ([a-z]+) of 680',
                 r'ranks ([a-z]+) of 680, \$p='], 0),
    ('heldout_err_min',
     [r'without exception, by ([\d.]+) to [\d.]+\s+log units'], 0.0006),
    ('heldout_err_max',
     [r'without exception, by [\d.]+ to ([\d.]+)\s+log units'], 0.0006),
    ('sn_actual', [r'reaches \$p<0\.05\$ in ([\d.]+)\\% of them',
     r"incumbent designation's (\d+\.\d+)\\% ranks fifth"], 0.06),
    ('fam_minp_p35', [r'Smallest nominal \$p\$ \(18 tests\)\s*& (\d+\.\d+) &'], 0.0006),
    ('fam_minp_p35lit', [r'Smallest nominal \$p\$ \(18 tests\)\s*&\s*'
                         r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0006),
    ('fam_minp_p28', [r'Smallest nominal \$p\$ \(18 tests\)\s*&\s*'
                      r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0006),
    ('fam_bh18_p35', [r'Benjamini-Hochberg over 18\s*& (\d+\.\d+) &'], 0.0006),
    ('fam_bh18_p35lit', [r'Benjamini-Hochberg over 18\s*& [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_bh18_p28', [r'Benjamini-Hochberg over 18\s*& [\d.]+ & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_bh6_p35', [r'Benjamini-Hochberg over 6, per space & (\d+\.\d+) &'], 0.0006),
    ('fam_bh6_p35lit', [r'Benjamini-Hochberg over 6, per space & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_bh6_p28', [r'Benjamini-Hochberg over 6, per space & [\d.]+ & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_wy_p35', [r'Westfall-Young min-\$P\$ over 18\s*& (\d+\.\d+) &'], 0.0006),
    ('fam_wy_p35lit', [r'Westfall-Young min-\$P\$ over 18\s*& [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_wy_p28', [r'Westfall-Young min-\$P\$ over 18\s*& [\d.]+ & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_meff_p35', [r'of 18 \(Li-Ji\)\s*& ([\d.]+)\s*&'], 0.006),
    ('fam_meff_p35lit', [r'of 18 \(Li-Ji\)\s*& [\d.]+\s*& ([\d.]+)\s*&'], 0.006),
    ('fam_meff_p28', [r'of 18 \(Li-Ji\)\s*& [\d.]+\s*& [\d.]+\s*& ([\d.]+)\s*&'], 0.006),
    ('fam_bonf_p35', [r'Bonferroni at \$m_\{\\mathrm\{eff\}\}\$\s*& (\d+\.\d+) &'], 0.0006),
    ('fam_bonf_p35lit', [r'Bonferroni at \$m_\{\\mathrm\{eff\}\}\$\s*& [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_bonf_p28', [r'Bonferroni at \$m_\{\\mathrm\{eff\}\}\$\s*& [\d.]+ & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_pre_p35', [r'Statistic specified in advance \(3 tests\) & (\d+\.\d+) &'], 0.0006),
    ('fam_pre_p35lit', [r'Statistic specified in advance \(3 tests\) & [\d.]+ & (\d+\.\d+) &'], 0.0006),
    ('fam_pre_p28', [r'Statistic specified in advance \(3 tests\) & [\d.]+ & [\d.]+ & \\textbf\{([\d.]+)\}'], 0.0006),
    # 2026-11-02 (T2-8). SI 문단을 다시 썼다. 종전 문장은 두 처리가 같은
    # 답을 낸다고 적었는데 앙상블 쪽이 3건과 4건으로 갈린다. 어법을 새
    # 문장에 맞추고, 갈리는 수 자체도 사실표에 걸어 둔다.
    ('conf_ens_lo',
     [r'runs from (\d+\.\d+) to \d+\.\d+ across the ensemble seeds'], 0.0006),
    ('conf_ens_hi',
     [r'runs from \d+\.\d+ to (\d+\.\d+) across the ensemble seeds'], 0.0006),
    ('conf_sgl_lo',
     [r'and from (\d+\.\d+) to\s*\n?\d+\.\d+ across the single-conformer'], 0.0006),
    ('conf_sgl_hi',
     [r'and from \d+\.\d+ to\s*\n?(\d+\.\d+) across the single-conformer'], 0.0006),
    ('conf_p_hi',
     [r'permutation value runs from 0\.0250 to (\d+\.\d+) under both'], 0.0006),
    # 'against four under the other three' 의 four 는 매칭 수이고 three 는
    # 시드 수다. 두 앵커를 섞지 않는다.
    # 2026-11-03 (T1-14). 본문 2.4 의 인코딩 단락과 안정성 단락.
    ('sp_n_singleton',
     [r'(\w+) of the thirty-eight retained coordinates'
      r'(?:[^.]{0,90}?)are binary indicators'], 0),
    ('sp_pedestal', [r'here (\d+), to the squared distance'], 0),
    ('sp_tabun_share',
     [r'[Ff]or tabun (?:it is|that addition is) ([\d.]+)\\%'], 0.05),
    ('sp_tabun_base', [r'moves tabun from\s*\n?([\d.]+) to [\d.]+ and the count'],
     0.0005),
    ('sp_tabun_drop1', [r'moves tabun from\s*\n?[\d.]+ to ([\d.]+) and the count'],
     0.0005),
    ('sp_match_drop1',
     [r'from three valid matches to (\w+)[.;] [Dd]ropping'], 0),
    ('sp_match_drop3', [r'three compounds or fewer gives (\w+)[,;]'], 0),
    ('sp_match_lowcard', [r'five or fewer distinct values gives (\w+),'], 0),
    ('sp_p_base', [r'one-sided \$p\$ falling from ([\d.]+) to [\d.]+'], 0.0005),
    ('sp_p_lowcard', [r'one-sided \$p\$ falling from [\d.]+ to ([\d.]+)'], 0.0005),
    ('sp_mantel_q_base',
     [r'value is ([\d.]+) under the baseline encoding',
      r'value between ([\d.]+) and [\d.]+'], 0.0005),
    ('sp_mantel_q_drop1',
     [r'baseline encoding and\s*\n?([\d.]+), [\d.]+, [\d.]+, [\d.]+, '
      r'[\d.]+ and [\d.]+ under the six'], 0.0005),
    ('sp_mantel_q_drop2',
     [r'baseline encoding and\s*\n?[\d.]+, ([\d.]+), [\d.]+, [\d.]+, '
      r'[\d.]+ and [\d.]+ under the six'], 0.0005),
    ('sp_mantel_q_drop3',
     [r'baseline encoding and\s*\n?[\d.]+, [\d.]+, ([\d.]+), [\d.]+, '
      r'[\d.]+ and [\d.]+ under the six'], 0.0005),
    ('sp_mantel_q_drop4',
     [r'baseline encoding and\s*\n?[\d.]+, [\d.]+, [\d.]+, ([\d.]+), '
      r'[\d.]+ and [\d.]+ under the six'], 0.0005),
    ('sp_mantel_q_drop5',
     [r'baseline encoding and\s*\n?[\d.]+, [\d.]+, [\d.]+, [\d.]+, '
      r'([\d.]+) and [\d.]+ under the six'], 0.0005),
    ('sp_mantel_q_lowcard',
     [r'baseline encoding and\s*\n?[\d.]+, [\d.]+, [\d.]+, [\d.]+, '
      r'[\d.]+ and ([\d.]+) under the six'], 0.0005),
    ('sp_mantel_q_max',
     [r'value between [\d.]+ and ([\d.]+)'], 0.0005),
    ('sp_mantel_nvar',
     [r'Mantel tests under all (\w+)\s*\n?encodings',
      r'All (\w+) are now recomputed'], 0),
    # 2026-11-03. 신설 문장들의 수치. 자리별 맹점을 만들지 않기 위해
    # 새로 쓴 모든 수치에 앵커를 붙인다.
    ('mantel_bh12_main',
     [r'over the twelve leaves a smallest adjusted\s*\n?value of ([\d.]+)\.'],
     0.0005),
    ('mantel_d2uv_struct_p',
     [r'structural space\s*\n?in the ultraviolet channel at ([\d.]+)'], 0.0005),
    ('mantel_d2uv_merged',
     [r'merged space in the same channel\s*\n?at ([\d.]+)'], 0.0005),
    ('mantel_d2nor_elec_p',
     [r'electronic space in the reflectance channel at ([\d.]+)\.'], 0.0005),
    ('d1_rho_max_abs',
     [r'in any space or channel, with \$\|\\rho\| \\le ([\d.]+)\$'], 0.0015),
    ('d1_p_min',
     [r'every \$p\$ at or above ([\d.]+)\.'], 0.0015),
    # 2026-11-03. 주입시험이 SI 의 타분 거리에 앵커가 없다는 것을 드러냈다.
    # 2.0038 은 세 군데(본문 2.4, SI 열순서 절, SI 응답공간 절)에 적혀 있다.
    ('t3_GA_full',
     [r'its normalised distance is\s*\n?([\d.]+) against a threshold of 2',
      r'crosses the threshold of 2:\s*\n?[\d.]+ under the old filter against '
      r'([\d.]+) under the new'], 0.00005),
    ('t3_GA_legacy',
     [r'crosses the threshold of 2:\s*\n?([\d.]+) under the old filter'],
     0.00005),
    # 2026-11-06 (9차). 커버레터를 본문과 같은 두 임계 서술로 고치면서
    # 'where the largest observed is' 앞 문구가 바뀌었다. 수치 주변만 남긴다.
    ('mantel_rho_max_abs',
     [r'(?:where the|the)\s+largest\s*\n?\s*observed is ([\d.]+)'], 0.005),
    ('conc_norm_rho',
     [r'at a Spearman correlation of exactly\s*\n?\$\+([\d.]+)\$ in every '
      r'one of the \w+ run-channel cells'], 0.0005),
    # 초록의 MDE 와 최대 관측 rho.
    ('mde_rho_mantel_bh',
     [r'resolves a rank\s*\n?correlation of ([\d.]+) at 80\\% power at the corrected',
      r'resolves a rank\s*\n?correlation of ([\d.]+) at 80\\% power at that corrected',
      r'that threshold the same simulation puts 80\\% power at\s*\n?\$\|\\rho\| = ([\d.]+)\$'], 0.006),
    ('mde_rho_mantel', [r'and ([\d.]+) at a nominal \$\\alpha\$'],
     0.005),
    ('conf_ens_match_hi',
     [r'against (\w+) under the\s*\n?other \w+\.'], 0),
    ('conf_ens_n4',
     [r'against \w+ under the\s*\n?other (\w+)\.'], 0),
    ('conf_sgl_match',
     [r'the count is (\w+) under all six single-conformer'], 0),
    ('conf_seeds', [r'under (\w+) seeds and rebuilding'], 0),
    ('mantel_bh12', [r'smallest adjusted value is\s*(\d+\.\d+) and no cell',
                     r'adjusted value from \d+\.\d+ to\s*(\d+\.\d+)\.'], 0.0006),
    ('mantel_bh8', [r'adjusted value from (\d+\.\d+) to'], 0.0006),
    # 2026-10-27. 종전 패턴은 'one-sided $p$ of (\d+\.\d+)' 였다. 병합공간
    # 대조 문장이 같은 어구를 쓰므로 0.0221 을 GFN2 값으로 읽었다. 앞 절을
    # 붙여 한 문장으로 묶는다.
    ('dft_es_gfn2_p',
     [r'one valid match of eleven and\s+a one-sided \$p\$ of (\d+\.\d+)',
      r'moves\s+it from (\d+\.\d+) to \d+\.\d+'], 0.006),
    ('dft_es_gfn2_match', [r'(\w+) valid match of eleven'], 0),
    ('dft_es_b3lyp_p',
     [r'dipole at GFN2 gives four valid matches and (\d+\.\d+)'], 0.006),
    ('dft_es_b3lyp_match',
     [r'dipole at GFN2 gives (\w+) valid matches'], 0),
    ('dft_es_full_p',
     [r'Replacing all three gives\s+five valid matches and (\d+\.\d+)',
      r'moves\s+it from \d+\.\d+ to (\d+\.\d+)',
      r'return (\d+\.\d+)\.\s+We do not promote'], 0.0006),
    ('dft_es_full_match',
     [r'Replacing all three gives\s+(\w+) valid matches',
      r'from one valid match to (\w+)'], 0),
    # 2026-11-06 (9차 B5). 종전 앵커는 'against four' 를 글자로 박고 있었다.
    # 예치물의 merged_gfn2.n_match 는 3 인데 원고가 four 라고 적고 있었고,
    # 앵커가 그 오기를 전제로 짜여 있어서 고치는 순간 두 사실이 미발견으로
    # 빠졌다. 개수를 글자로 박지 말고 사실로 묶는다.
    ('dft_mg_b3lyp_p',
     [r'a one-sided \$p\$ of (\d+\.\d+), against [a-z]+ and',
      r'holds\s+at (\d+\.\d+)'], 0.0006),
    ('dft_mg_b3lyp_match',
     [r'analysis and returns (\w+) valid matches and\s+a one-sided'], 0),
    ('dft_mg_gfn2_match',
     [r'of \d+\.\d+, against ([a-z]+) and \d+\.\d+ for the GFN2 control'], 0),
    ('dft_mg_gfn2_p',
     [r'against [a-z]+ and (\d+\.\d+) for the GFN2 control'], 0.0006),
    ('nopatho_struct_p', [r'to three decimal places \((\d+\.\d+),'], 0.006),
    ('nopatho_elec_p', [r'three decimal places \(\d+\.\d+, (\d+\.\d+) and'], 0.006),
    ('nopatho_merged_p', [r'three decimal places \(\d+\.\d+, \d+\.\d+ and (\d+\.\d+)\)'], 0.006),
    ('gfn2_omega_min', [r'they run from ([\d.]+) to\s*\n?[\d.]+~eV across this panel'], 0.006),
    ('gfn2_omega_max', [r'they run from [\d.]+ to\s*\n?([\d.]+)~eV across this panel'], 0.006),
    # SI 쌍극자 일치도 (2026-11-02 신설).
    ('dip_r', [r'agree at a Pearson \$r\$ of (\d+\.\d+)'], 0.0002),
    ('dip_mae', [r'mean absolute difference of (\d+\.\d+)~D'], 0.0006),
    ('dip_dmmp_gfn2', [r'DMMP is (\d+\.\d+)~D under GFN2'], 0.0006),
    ('dip_dmmp_b3lyp',
     [r'under GFN2 and\s*\n?(\d+\.\d+)~D under the density functional'], 0.0006),
    ('dip_maxdev', [r'a difference of (\d+\.\d+)~D, which is'], 0.0006),
    # ACS 본문 §2.4 전용 (2026-11-02). SI 와 어법이 달라 자리를 따로 건다.
    ('mt_dft_gfn2_p',
     [r'eleven and a one-sided \$p\$ of (\d+\.\d+) to \w+ matches'], 0.0005),
    ('mt_dft_full_p',
     [r'of \d+\.\d+ to \w+ matches and (\d+\.\d+)[.,]'], 0.0005),
    ('mt_dft_full_match',
     [r'of \d+\.\d+ to (\w+) matches and \d+\.\d+[.,]'], 0),
    ('mt_dft_mg_match',
     [r'coordinates, (\w+) valid matches and a one-sided',
      r'keeps (\w+) valid matches and a one-sided'], 0),
    ('mt_dft_mg_p',
     [r'valid matches and a one-sided \$p\$ of\s*\n?(\d+\.\d+) under both'],
     0.0005),
    ('mt_dft_mg_gfn2_p',
     [r'valid matches and a one-sided \$p\$ of\s*\n?(\d+\.\d+) under both'],
     0.0005),
    ('dft_omega_min', [r'range becomes ([\d.]+) to [\d.]+~eV'], 0.006),
    ('dft_omega_max', [r'range becomes [\d.]+ to ([\d.]+)~eV'], 0.006),
    ('dft_ps', [r'chloropicrin falling from [\d.]+ to ([\d.]+)'], 0.006),
    ('dft_cg', [r'phosgene from [\d.]+ to ([\d.]+)'], 0.006),
    ('dft_dmmp', [r'DMMP at ([\d.]+) against sarin'], 0.006),
    ('dft_gb', [r'DMMP at [\d.]+ against sarin at ([\d.]+)'], 0.006),
    ('dft_itf_above_gb', [r'(\w+) of the fourteen interferents still score above'], 0),
    ('dft_lo1', [r'three lowest agents become VX\s*\n?\(([\d.]+)\)'], 0.006),
    ('dft_lo2', [r'sulfur mustard \(([\d.]+)\) and cyclosarin'], 0.006),
    ('dft_lo3', [r'and cyclosarin\s*\n?\(([\d.]+)\) rather than'], 0.006),
    ('dft_rho_omega', [r'two sets of \$\\omega\$ values is ([\d.]+), against'], 0.006),
    ('dft_rho_lumo', [r'is [\d.]+, against\s*\n?([\d.]+) for the underlying LUMO'], 0.006),
    ('dft_n_dft', [r'recomputed the whole\s*\n?(\w+[-\w]*)-compound panel at'], 0),
    ('ascl3_a28_match', [r'Twenty-eight, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a28_dist', [r'Twenty-eight, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a35_match', [r'Thirty-five, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a35_dist', [r'Thirty-five, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a35l_dist', [r'Thirty-five, incumbent 3 \+ literature 5 & \d+/11 & CEES & ([\d.]+)'], 0.006),
    ('ascl3_a36_match', [r'Thirty-six, incumbent three & (\d+)/11 & dichloromethane'], 0),
    ('ascl3_a36_dist', [r'Thirty-six, incumbent three & \d+/11 & dichloromethane & ([\d.]+)'], 0.006),
    ('ascl3_a36l_match', [r'Thirty-six, incumbent 3 \+ literature 6 & (\d+)/11 & AsCl'], 0),
    ('null3_lower', [r'the three at \$p=([\d.]+)\$ over all',
                     r'zero sits at \$p=([\d.]+)\$ against an expected'], 0.006),
    ('null8_nested_mean', [r'gives a mean of ([\d.]+) rather',
                           r'against\s*\n?an expected ([\d.]+) once the null'], 0.006),
    ('null8_nested_p', [r'observed five sits at \$p=([\d.]+)\$ against it',
                        r'the five at \$p=([\d.]+)\$ against'], 0.006),
    ('n_eff_d2uv', [r'leaves (\d+) effective channels',
                    r'of which (\d+) are independent'], 0),
    ('n_eff_d3uv', [r'effective channels in the second\s*\n?set'
                    r"'s ultraviolet data and (\d+) in the third"], 0),
    ('mde_rho_compound', [r'as the unit, which is the unit the Mantel test uses,\s*\n?'
                          r'returns \$([\d.]+)\$'], 0.0015),
    ('mde_rho_mantel', [r'80\\% power is reached at \$\|\\rho\| = ([\d.]+)\$',
                        r'design detects \$\|\\rho\| \\ge ([\d.]+)\$ at 80\\% power',
                        # 2026-11-06 (9차). 주입시험이 커버레터의 MDE 를
                        # 무검사로 잡아냈다. 본문은 0.85/0.67 로 고쳤는데
                        # 커버레터만 0.67 단독으로 남아 있었다.
                        r'and ([\d.]+) at a nominal \$\\alpha\$'], 0.006),
    ('mantel_power_perm', [r'all \$6! = (\d+)\$ relabellings',
                           r'exact Mantel over all (\d+) relabellings'], 0),
    ('mantel_type_i', [r'type~I rate of the simulated test at ([\d.]+)'], 0.0),
    ('oi_legacy28_sets', [r'the old filter returned\s+([\w-]+) distinct retained sets'], 0),
    ('oi_legacy28_m3', [r'valid-match count of three in ([\w-]+)\s*\n?of the forty orderings'], 0),
    ('oi_legacy28_m4', [r'four in the remaining ([\w-]+)\.'], 0),
    ('oi_legacy35_sets', [r'differently: ([\w-]+) distinct retained sets'], 0),
    ('vol_withdrawn_p', [r'returned\s*\n?\$?p ?= ?(\d+\.\d+)\$? and was withdrawn',
                         r'It returns \$p = (\d+\.\d+)\$, and',
                         r'which returned\s*\n?\$p = (\d+\.\d+)\$'], 0.0005),
    ('vol_four_cell_bh', [r'one family of four\s*\n?cells it is (\d+\.\d+)',
                          r'read as one family give (\d+\.\d+)',
                          r'as one family of four cells instead, the smallest adjusted value is\s*\n?(\d+\.\d+)'], 0.0005),
    ('mantel_type_i_nrep', [r'over ([\d,\\\s]+) null replicates'], 0),
    ('mantel_d2nor_elec', [r'reflectance channel \(\$\\rho = \+\d+\.\d+\$,\s*\n?Mantel \$p = (\d+\.\d+)\$'], 0.0005),
    ('mantel_d2uv_merged_rho',
     [r"merged space in the\s*\n?second compound set's ultraviolet channel "
      r"\(\$\\rho = ([-+][\d.]+)\$"], 0.0015),
    # 2026-10-27. 종전 패턴은 rho 값(+0.518)을 박아 두었다. 재실행으로 rho 가
    # +0.514 가 되자 검사가 조용히 꺼졌다. 값이 아니라 문맥으로 잡는다.
    ('mantel_d2uv_merged', [r'ultraviolet channel \(\$\\rho = \+\d+\.\d+\$, Mantel \$p = (\d+\.\d+)\$'], 0.0005),
    ('mantel_d2uv_struct_rho', [r'that same set and channel \(\$\\rho = ([-+][\d.]+)\$'], 0.0015),
    # 2026-11-03 (T1-13). SI 문장을 여섯 셀로 늘리면서 어법이 바뀌었다.
    # 2026-11-03 (작업 C) 폐기. 중심화 '바닥' 논증을 두 본문 모두에서
    # 철회했으므로 바닥 여유 앵커 5건은 가리킬 문장이 없다. 사실표의
    # floor_* 값은 비교용으로 남겨 두되 본문 대조는 하지 않는다.
    ('prem_pairwise_viol',
     [r'in (\d+) of the 34\s*\n?compound-and-cell instances',
      r'in (\d+) of the 34 compound-and-cell instances'], 0),
    ('dmmp_ga_uv', [r'tabun, at \$\+?([\d.]+)\$ against sarin',
                    r'correlates with tabun at \$\+([\d.]+)\$'], 0.0015),
    ('dmmp_gb_uv', [r'against sarin\'s \$\+?([\d.]+)\$',
                    r'and with\s+sarin at \$\+([\d.]+)\$'], 0.0015),
    # 도판 (d) 에 박힌 p 도 같은 사실표와 대조한다.
    ('p_struct_fig', [r'A_rdkit[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_elec_fig', [r'B_electronic[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p_merged_fig', [r'C_merged[^\n]*exact p = ([\d.]+)'], 0.0005),
    ('p28_struct', [r'-dimensional value of (\d+\.\d+) sits at the',
                    r'from \$p=(\d+\.\d+)\$ to\s*\n?\$p=[\d.]+\$ on the twenty-eight',
                    r'structural \$p\$ of (\d+\.\d+),'], 0.0005),
    # 2026-11-04 (J-C). ACS본 어법으로 교정.
    ('rm_p28', [r'`p` of ([\d.]+) on the twenty-eight-compound panel'], 0.00005),
    ('rm_p35', [r'([\d.]+) on the primary\s*\n?thirty-five-compound panel'],
     0.00005),
    ('rm_p35lit', [r'([\d.]+) on the panel with eight simulants'], 0.00005),
    ('rm_q28', [r'Hochberg minima of ([\d.]+), [\d.]+ and [\d.]+'], 0.00005),
    ('rm_q35', [r'Hochberg minima of [\d.]+, ([\d.]+) and [\d.]+'], 0.00005),
    ('rm_q35lit', [r'Hochberg minima of [\d.]+, [\d.]+ and ([\d.]+)'], 0.00005),
    ('len_pages', [r'runs to ([\d.]+) journal pages against'], 0.05),
    ('len_body', [r'specify: ([\d,]+) words\s*\n?of text'], 0),
    ('len_refs', [r'and ([\d,]+) words of\s*\n?references'], 0),
    # 2026-11-07 (12차 7단계). 커버레터가 산출물 수치를 그대로 찍으므로
    # 숫자도 받는다.
    ('len_figs', [r'([a-z]+|\d+) double-column figures at 500 words'], 0),
    ('coll_merged_off',
     [r'the merged space from\s*\n?\$\\rho = \+([\d.]+)\$ to'], 0.00005),
    ('coll_merged_on',
     [r'\$\\rho = \+[\d.]+\$ to \$\+([\d.]+)\$ and the structural'], 0.00005),
    ('coll_struct_off',
     [r'the structural space from \$\+([\d.]+)\$ to'], 0.00005),
    ('coll_struct_on',
     [r'the structural space from \$\+[\d.]+\$ to\s*\n?\$\+([\d.]+)\$'],
     0.00005),
    ('coll_q_off',
     [r'over the twelve cells from\s*\n?([\d.]+) to [\d.]+\.'], 0.00005),
    ('coll_q_on',
     [r'over the twelve cells from\s*\n?[\d.]+ to ([\d.]+)\.'], 0.00005),
    ('thr_free_28',
     [r'twenty-eight-compound count from three valid matches to (\w+)',
      r'takes the twenty-eight-compound count from three to (\w+), while'], 0),
    ('thr_null_28', [r'any incumbent is ([\d.]+) of\s*\n?eleven on the '
                     r'twenty-eight-compound panel'], 0.005),
    ('thr_null_35', [r'and ([\d.]+) of eleven on the primary\s*\n?one'], 0.006),
    ('gate_checked', [r'(\d+) columns in each, out of \d+ to \d+ numeric',
                      r'reproduces from the code, (\d+) of the \d+ numeric'], 0),
    ('gate_numeric', [r'out of (\d+) to \d+ numeric columns',
                      r'reproduces from the code, \d+ of the (\d+) numeric'], 0),
    ('gate_numeric_max', [r'out of \d+ to (\d+) numeric columns'], 0),
    ('gate_caches', [r'the (\w+) caches the analysis reads'], 0),
    # 2026-11-07 (12차 A-9) 염료 패널.
    # 2026-11-07 (12차 A-9). 본문과 SI 가 어법이 다르다. 둘 다 건다.
    ('dye_n_nondye', [r'(\w+) of the\s*\n?twenty-nine are not chromogenic',
                      r'(\w+) are not chromogenic dyes: quinine'], 0),
    ('dye_n_fluor', [r'twenty-nine are not chromogenic dyes: (\w+) are'], 0),
    ('dye_top_uv', [r'maximum \$\\Delta E_\{00\}\$ of ([\d.]+) in the'], 0.05),
    ('dye_top_normal', [r'against ([\d.]+) in reflectance, a ratio'], 0.05),
    ('dye_top_ratio', [r'a ratio of ([\d.]+), the largest'], 0.05),
    ('dye_anthracene_ratio', [r'pyrene follow at ([\d.]+) and [\d.]+'], 0.05),
    ('dye_pyrene_ratio', [r'pyrene follow at [\d.]+ and ([\d.]+)'], 0.05),
    # 2026-11-07 (12차 A-5) 양수 가상궤도.
    ('lumo_b3lyp_pos',
     [r'In (\w+) of the twenty-eight compounds the B3LYP LUMO'], 0),
    ('lumo_gfn2_pos',
     [r'(\w+)\s*\n?of the twenty-eight GFN2-xTB LUMO energies are positive'], 0),
    ('lumo_b3lyp_max', [r'\$n\$-hexane at \$\+([\d.]+)\$~eV'], 0.005),
    ('gate_exempt', [r'The (\w+)\s*\n?columns it does not recompute'], 0),
    ('pct_struct_caption', [r'is the (\d+)th percentile of those draws',
                            r'observed value at the (\d+\.\d+)th percentile',
                            r'place it at the\s*\n?(\d+\.\d+)th percentile'], 0.06),
    ('pct_struct', [r'observed value at the ([\d.]+)(?:st|nd|rd|th) percentile of that',
     r'-dimensional value of [\d.]+ sits at the ([\d.]+)(?:nd|th|st|rd)\s*\n?percentile'], 0.06),
    ('pct_elec', [r'electronic value of [\d.]+ (?:sits )?at the ([\d.]+)(?:nd|th|st|rd)'], 0.06),
    ('pct_pc3', [r'three-principal-component value of [\d.]+ sits at the ([\d.]+)(?:st|nd|rd|th)\s*\n?percentile'], 0.06),
    ('mantel_minq', [r'smallest adjusted value is\s*\n?([\d.]+) and no cell',
                     r'over the twelve leaves a smallest adjusted\s+value of ([\d.]+)',
                     r'smallest \$q\$ = ([\d.]+)\)'], 0.0005),
    ('white_frac28', [r'the rate falls from [\d.]+\\% to ([\d.]+)\\%'], 0.06),
    # 2026-11-03 (2차) 폐기. 이 앵커는 폐기 수치 0.094 를 닻으로 썼다.
    # 0.094 를 0.149 로 교정하자 꺼졌다. white_med28_2 가 새 어법 자리에서
    # 같은 값을 검사한다.
    ('q_by_space28', [r'over 6, per space & [\d.]+ & [\d.]+ & ([\d.]+)'], 0.0005),
    ('q_prespec28', [r'returns a value below 0\.05, at ([\d.]+)',
     r'and ([\d.]+) for the only scheme that fixes',
     r'would survive at\s*\n?\$q=([\d.]+)\$'], 0.0005),
    ('meff28', [r'estimates ([\d.]+) effective tests'], 0.005),
    ('meff_p28', [r'and gives\s*\n?([\d.]+), conservative here by a factor'], 0.0005),
    ('match35_inc', [r'twenty-eight-compound panel and (\w+) on the primary one'], 0),
    ('match35_lit', [r'Admitting the\s*\n?five as simulants gives (\w+) matched agents'], 0),
    ('match28_inc', [r'agents of eleven, against\s*\n?(\w+) on the twenty-eight-compound panel'], 0),
    ('dfp_gb', [r'fluorophosphate is present, at ([\d.]+),'], 0.005),
    ('dfp_gd', [r'is present, at [\d.]+, ([\d.]+) and'], 0.005),
    ('dfp_gf', [r'is present, at [\d.]+, [\d.]+ and ([\d.]+) from'], 0.005),
    ('dmmp_gb_35', [r'on that panel sarin stands\s*\n?([\d.]+) from DMMP'], 0.005),
    ('dmmp_gd_35', [r'from DMMP, soman ([\d.]+) and'], 0.005),
    ('dmmp_gf_35', [r'soman [\d.]+ and cyclosarin ([\d.]+)\.'], 0.005),
    ('cees_hd', [r'sulfide at ([\d.]+) from sulfur mustard'], 0.005),
    # 최소검출효과 (2026-10-17). 본문 Limitations + SI + 도판 주석 3경로 대조.
    ('mde_delta35', [r'moved ([\d.]+) typical nearest-neighbour',
                     r'returns \$p<0\.05\$ at\s*\n?\$\\delta = ([\d.]+)\$',
                     r'"mde_delta35": "([\d.]+)"'], 0.005),
    ('mde_gap_median', [r'carry has a median of\s*\n?([\d.]+)',
                        r'other compound has a median of (0\.\d+)',
                        r'median of ([\d.]+), more than twice',
                        r'"gap_median35": "([\d.]+)"'], 0.0005),
    ('mde_gap_min', [r'only at sarin, where it is\s*\n?([\d.]+)',
                     r'it is ([\d.]+), does it fall below',
                     r'"gap_min35": "([\d.]+)"'], 0.0005),
    ('mde_p0_35', [r'it returns ([\d.]+) on\s*\n?the primary panel',
                   r'"p_at_zero35": "([\d.]+)"'], 0.0005),
    ('mde_p0_28', [r'and ([\d.]+) on the twenty-eight-compound panel, reproducing',
                   r'"p_at_zero28": "([\d.]+)"'], 0.0005),
    # 9.625 의 소수 둘째자리 반올림 경계이므로 허용오차를 0.006 으로 둔다
    ('null_exp35', [r'expectation of ([\d.]+) set by the 21:3'], 0.006),
    ('null_exp28', [r'smaller panel against ([\d.]+)'], 0.005),
    ('wy35', [r'Westfall-Young min-\$P\$ over 18\s+& (\d+\.\d+) &'], 0.0005),
    ('wy28', [r'handles the dependence between tests exactly,\s+gives ([\d.]+)',
     r'Westfall-Young min-\$P\$ over 18\s+& [\d.]+ & [\d.]+ & ([\d.]+)'], 0.0005),
    ('wy35lit', [r'Westfall-Young min-\$P\$ over 18\s+& [\d.]+ & (\d+\.\d+) &'], 0.0005),
    # 2026-11-06 (9차). Figure 1(b) 의 재표집 수치를 본문으로 올리면서
    # 면제를 해제하고 정규식을 새 문장에 맞춘다.
    ('boot_inc_zero', [r'no match at all in ([\d.]+)\\% of replicates',
                       r'exactly zero\s*\n?in ([\d.]+)\\% of replicates'], 0.05),
    ('boot_inc_le1', [r'at most one in ([\d.]+)\\%'], 0.05),
    ('boot_diff_pos', [r'positive in (\d+) of 1000'], 0),
    ('ret_gb3', [r'as its nearest compound in ([\d.]+)\\% of replicates'], 0.05),
    ('ret_gd3', [r'replicates, soman in ([\d.]+)\\%'], 0.05),
    ('ret_gf3', [r'and\s*\n?cyclosarin in ([\d.]+)\\%'], 0.05),
    ('ret_ga8', [r'cyanophosphonate in ([\d.]+)\\%'], 0.05),
    ('phorate_fold', [r'is predicted (\d[\d,]*) times\s*\n?too safe'], 1),
    ('malathion_fold', [r'is predicted (\d[\d,]*) times too toxic'], 1),
    # 표 3 (tab:match) 전수 대조. 2026-10-17 신설. 이 표는 그때까지 검증기에
    # 등록되지 않아 11행 중 9행이 구판 값인 채로 144건 검사를 통과했다.
    ('t3_GB', [r'GB & \S+ & (\d+\.\d+) &', r'Sarin \(([\d.]+)\), soman',
              r'where sarin sits at (\d+\.\d+)'], 0.0005),
    ('t3_GD', [r'GD & \S+ & (\d+\.\d+) &',
               r'Sarin \([\d.]+\), soman \(([\d.]+)\)'], 0.0005),
    ('t3_CK', [r'CK & \S+ & (\d+\.\d+) &'], 0.0005),
    ('t3_GF', [r'GF & \S+ & (\d+\.\d+) &',
               r'cyclosarin \(([\d.]+)\) are matched'], 0.0005),
    ('t3_AC', [r'AC & \S+ & (\d+\.\d+) &'], 0.0005),
    ('t3_CG', [r'CG & \S+ & (\d+\.\d+) &',
               r'merged space \([\d.]+ against ([\d.]+)\)'], 0.0005),
    ('t3_HD', [r'HD & \S+ & (\d+\.\d+) &'], 0.0005),
    ('t3_GA', [r'GA & \S+ & (\d+\.\d+) &', r'tabun at ([\d.]+) and VX',
               r'passes its value of (\d+\.\d+)'], 0.0005),
    ('t3_L', [r'L & \S+ & (\d+\.\d+) &'], 0.0005),
    ('t3_PS', [r'PS & \S+ & (\d+\.\d+) &'], 0.0005),
    ('t3_VX', [r'VX & \S+ & (\d+\.\d+) &', r'VX at\s*\n?([\d.]+) are not'], 0.0005),
    # res:match 산문의 공간별 수치. 종전 산문은 세 공간의 원거리와 정규화
    # 거리를 섞어 적었고 어느 값도 산출물과 대조되지 않았다.
    ('struct_GD', [r'soman at ([\d.]+)\s*\n?and cyclosarin at'], 0.0005),
    ('struct_GF', [r'cyclosarin at ([\d.]+) (?:falling )?outside'], 0.0005),
    ('struct_VX', [r'and ([\d.]+) in the structural\s*\n?space that is not a match'],
     0.0005),
    ('merged_VX', [r'at ([\d.]+) typical separations in the merged space'], 0.0005),
    ('elec_AC', [r'dichloromethane at ([\d.]+),\s*\n?because acetonitrile'], 0.0005),
    ('struct_CG', [r'structural space \([\d.]+ typical separations\s*\n?against ([\d.]+)\)'],
     0.005),
    ('minp35', [r'Smallest nominal \$p\$ \(18 tests\)\s+& (\d+\.\d+) &',
                r'eighteen-test family there is ([\d.]+)\s*\n?with equal weight per'],
     0.0005),
    ('minp35lit', [r'Smallest nominal \$p\$ \(18 tests\)\s*&\s*'
                   r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0005),
    ('minp28', [r'Smallest nominal \$p\$ \(18 tests\)\s*&\s*'
                r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?(?:[\d.]+|\$<\$0\.001)\}?\s*&\s*' r'(?:\\textbf\{)?([\d.]+)\}?'], 0.0005),
    ('elec_CG', [r'electronic space \([\d.]+ against ([\d.]+)\)'], 0.005),
    ('struct_CG_rival', [r'structural space \(([\d.]+) typical separations'], 0.005),
    ('merged_CG_rival', [r'and in the merged space \(([\d.]+) against'], 0.005),
    ('elec_CG_rival', [r'chloroform in the electronic space \(([\d.]+) against'], 0.005),
    # 블록 가중 민감도 (2026-10-17). SI 표 + SI 산문 + 본문 한 문장 3경로.
    ('bb_nstruct', [r'Weighting the (\d+) structural', r'\$n_\{s\}=(\d+)\$',
                    r'puts (\d+) structural coordinates'], 0),
    ('bb_nelec', [r'and (\d+) electronic descriptors as two equal blocks',
                  r'\$n_\{e\}=(\d+)\$'], 0),
    # 2026-11-07 (12차 A-3). 종전 두 패턴이 0\.151 을 리터럴로 들고 있었고
    # 그래서 그 값은 검사 대상 밖이었다. 리터럴을 지우고 두 값을 각각
    # 자기 사실에 건다.
    ('bb_minp35', [r'lowers that value to ([\d.]+)',
                   r'nominal value from [\d.]+ to ([\d.]+)',
                   r'Thirty-five, incumbent three & \d+/11 & \d+/11 & [\d.]+ & (\d+\.\d+) &',
                   r'with equal weight per descriptor and ([\d.]+) with'],
     0.0005),
    ('bb_minpdesc35',
     [r'nominal value from ([\d.]+) to [\d.]+',
      r'Thirty-five, incumbent three & \d+/11 & \d+/11 & (\d+\.\d+) &',
      r'([\d.]+) with equal weight per descriptor and'],
     0.0005),
    ('bb_minp35lit',
     [r'Thirty-five, incumbent 3 \+ literature 5 & \d+/11 & \d+/11 & [\d.]+ & (\d+\.\d+) &'], 0.0005),
    ('bb_minp28',
     [r'Twenty-eight, incumbent three & \d+/11 & \d+/11 & [\d.]+ & (\d+\.\d+) &'], 0.0005),
    ('bbq_35', [r'Thirty-five, incumbent three & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_35lit',
     [r'Thirty-five, incumbent 3 \+ literature 5 & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    ('bbq_28',
     [r'Twenty-eight, incumbent three & \d+/11 & \d+/11 & [\d.]+ & [\d.]+ & [\d.]+ & ([\d.]+)'],
     0.0005),
    # 2026-10-24. 유효매칭 열을 per descriptor / per block 두 칸으로 쪼갰다.
    # 둘 다 걸지 않으면 한쪽을 망가뜨리는 주입이 통과한다 (실제로 통과했다).
    ('bb_cur_match35', [r'Thirty-five, incumbent three & (\d+)/11'], 0),
    ('bb_cur_match35lit', [r'Thirty-five, incumbent 3 \+ literature 5 & (\d+)/11'], 0),
    ('bb_cur_match28', [r'Twenty-eight, incumbent three & (\d+)/11'], 0),
    ('bb_match35', [r'Thirty-five, incumbent three & \d+/11 & (\d+)/11'], 0),
    ('bb_match35lit', [r'Thirty-five, incumbent 3 \+ literature 5 & \d+/11 & (\d+)/11'], 0),
    ('bb_match28', [r'Twenty-eight, incumbent three & \d+/11 & (\d+)/11'], 0),
    ('bb_meandiff28', [r'there falls from 0\.025 to ([\d.]+)'], 0.0005),
    ('bbpre_28', [r'which stays at ([\d.]+) on the twenty-eight-compound panel'],
     0.0005),
    # 네 번째 구성 (2026-10-20). 표 1 네 번째 블록 + §3.1 산문 + Limitations.
    # 표 3 사건의 교훈대로 표를 먼저 등록한다.
    ('pp_st_md', [row_pat('35 cpd., ten designated', 'structural', 0)], 0.0005),
    ('pp_st_lr', [row_pat('35 cpd., ten designated', 'structural', 2)], 0.0005),
    ('pp_st_ds', [row_pat('35 cpd., ten designated', 'structural', 3)], 0.0005),
    # 2026-11-04 (작업 D). 라벨을 '지정 후보 총수' 로 통일했다. 종전
    # '+ literature 5' 는 2행과 4행이 같은 문자열이 되어 앵커가 엉뚱한
    # 행에 붙었다. 농약 구성은 4행이고 이제 라벨이 유일하다.
    ('pp_el_md', [row_pat('\\quad (3 incumbent + 5 literature', 'electronic', 0)], 0.0005),
    ('pp_el_mdn', [row_pat('\\quad (3 incumbent + 5 literature', 'electronic', 1)], 0.0005),
    ('pp_el_nv', [row_pat('\\quad (3 incumbent + 5 literature', 'electronic', 5)], 0.0005),
    ('pp_mg_md', [row_pat('\\quad + 2 pesticides)', 'merged', 0)], 0.0005),
    ('pp_mg_lr', [row_pat('\\quad + 2 pesticides)', 'merged', 2)], 0.0005),
    ('pp_mg_ds', [row_pat('\\quad + 2 pesticides)', 'merged', 3)], 0.0005),
    # pp_n_perm 은 UNQUOTED 면제 항목이다. 앵커를 중복으로 두면
    # 영구 '문맥 미발견' 이 남는다. 2026-11-03 제거.
    ('pp_min_p', [r'smallest nominal value of ([\d.]+) and a'], 0.00001),
    # 2026-11-03. tab:family 를 생성으로 돌리면서 칸 사이 공백이 바뀌어
    # 공백 개수를 고정한 이 패턴이 꺼졌다. \s* 로 느슨하게 한다.
    ('pf_meff', [r'\$m_\{\\mathrm\{eff\}\}\$ of 18 \(Li-Ji\)\s*&'
                 r'[\s\d.]+&[\s\d.]+&[\s\d.]+&\s*([\d.]+)'], 0.005),
    ('pf_prespec', [r'Statistic specified in advance \(3 tests\) & [\d.]+ & [\d.]+ & '
                    r'\\textbf\{[\d.]+\} & \\textbf\{([\d.]+)\}'], 0.0005),
    ('struct_removed', [r'remove (\d+) of the 51 on the primary panel'], 0),
    ('struct_removed_28', [r'of the 51 on the primary panel and\s*(\d+) on the'], 0),
    ('struct_kept_35', [r'leaving (\d+) and \d+ structural coordinates'], 0),
    ('struct_kept_28', [r'leaving \d+ and (\d+) structural coordinates'], 0),
    ('struct_candidates', [r'removes \d+ of the (\d+)',
                          r'physicochemical descriptors,\s*\n?(\d+) after excluding',
                           r'descriptors \((\d+) after excluding run metadata\)'], 0),
    ('pp_bh_q', [r'Benjamini-Hochberg value of ([\d.]+) over the same',
                 r'produces an adjusted\s*\n?(?:value of )?(\d+\.\d+)',
                 r'into one with an adjusted value of\s*\n?(\d+\.\d+)'], 0.00005),
    # 인용 문헌의 수치 (2026-10-20). outputs/cite_facts.json 이 단일 출처다.
    ('j18_identified', [r'identified\s*\n?(\d+) simulants'], 0),
    ('j18_unsuitable', [r'judged (\d+) unsuitable'], 0),
    ('j21_identified', [r'a later one screened (\d+) chemicals'], 0),
    ('j21_discounted', [r'discounted (\d+) on the risk'], 0),
    ('j21_suitable', [r'repository of (\d+) suitable'], 0),
    ('j21_possible', [r'suitable and (\d+) provisional'], 0),
    # ------------------------------------------- 휘발성 축 (2026-11-03 신설)
    # 본문 res:panels 말미 단락과 SI res:volatility, 그리고 생성표 두 개.
    ('vol_p_md', [r'one-sided \$p\$ of (0\.\d+) over all\s+680',
                  r'statistic gives a one-sided \$p\$ of (0\.\d+)'], 0.0005),
    # 'ranking 316th of 680' 과 구분해야 한다. 뒤따르는 문장을 닻으로 쓴다.
    ('vol_rank_md', [r'ranking the incumbent trio (\w+) of 680',
                     r'680 designations,\s*\n?ranking (\w+)\.'], 0),
    ('vol_obs_md', [r'margin is (0\.\d+) in standardised',
                    r'margin itself is (0\.\d+) in standardised'], 0.00005),
    ('vol_p_nvalid', [r'valid matches returns one match with a \$p\$ of '
                      r'(0\.\d+)',
                      r'one match with a \$p\$ of (0\.\d+)',
                      # 2026-11-03 (작업 B). 초록 편입분.
                      r'one valid match at \$p=(0\.\d+)\$'], 0.0005),
    ('vol_rank_nvalid', [r'ranking (\d+)\w*\s+of 680',
                         r'(\d+)\w*\s+of 680 and therefore'], 0),
    # 2026-11-03. 'all 680 designations' 는 본문 2.1 의 'all 2024
    # designations' 와 같은 어법이라 패널을 특정해야 한다. 휘발성 절의
    # 전수열거 수는 28종 패널의 것이므로 nperm28 이 이미 검사한다. 중복
    # 앵커를 두는 대신 휘발성 문장 쪽만 따로 가리킨다.
    ('vol_nperm', [r'exact enumeration of all (\d+) designations'], 0),
    ('vol_range', [r'spans (\d+\.\d+) decades'], 0.0005),
    # \w 는 하이픈을 먹지 않아 'twenty-five' 에서 'twenty' 만 잡혔다.
    ('vol_t1', [r'([\w-]+) of them at tier one'], 0),
    ('vol_gap', [r'gap of (0\.\d+) decades between',
                 r'on a gap of (0\.\d+) decades'], 0.00005),
    ('vol_bp_ratio_gap', [r'ratio condition once the gap exceeds (0\.\d+)'],
     0.00005),
    ('vol_bp_ratio_vp', [r'vapour pressure below\s+(0\.\d+)~mmHg'], 0.00005),
    ('vol_bp_closer_gap', [r'condition beyond\s+(0\.\d+) decades'], 0.00005),
    ('vol_bp_closer_vp', [r'decades, or (0\.\d+)~mmHg'], 0.00005),
    # 2026-11-03 (작업 B). 형태좌표 축소 민감도 문단.
    ('sr_nshape', [r'(\w+) of the thirty-eight coordinates'], 0),
    ('sr_dims28', [r'Nine of the (thirty-\w+) coordinates'], 0),
    ('sr_dims35', [r'nine of the ([\w-]+) on the primary'], 0),
    ('sr_maxr28', [r'block is (0\.\d+) on the twenty-eight'], 0.00005),
    ('sr_maxr35', [r'and (0\.\d+) on the primary\s*\n?panel'], 0.00005),
    ('sr_nderived', [r'with the ([\w-]+)\s*\n?derived coordinates dropped'], 0),
    ('sr_id_lo', [r'discrepancies of (\d+\.\d+) to [\d.]+\\%'], 0.05),
    ('sr_id_hi', [r'discrepancies of [\d.]+ to (\d+\.\d+)\\%'], 0.05),
    ('sr_share_full', [r'distance falls from (\d+\.\d+) to [\d.]+\\%'], 0.05),
    ('sr_share_red', [r'falls from [\d.]+ to (\d+\.\d+)\\%'], 0.05),
    ('sr_acn_hi_full', [r'it falls from (\d+\.\d+) and [\d.]+\\% to'], 0.05),
    ('sr_acn_lo_full', [r'falls from [\d.]+ and (\d+\.\d+)\\% to'], 0.05),
    ('sr_acn_hi_red', [r'to\s*\n?(0\.\d+) and 0\.\d+\\%\.'], 0.005),
    ('sr_acn_lo_red', [r'to\s*\n?0\.\d+ and (0\.\d+)\\%\.'], 0.005),
    ('sr_match_full', [r'matches there falls from ([\w-]+) to'], 0),
    ('sr_match_red', [r'falls from \w+ to ([\w-]+)\. On'], 0),
    ('sr_mantel_q_red', [r'shape coordinates raises it to (0\.\d+)'], 0.00005),
    # 전제검정 라벨순열 귀무 (초록).
    # 2026-11-07 (12차 Step 5). 초록에서 그 문장을 뺐다. 같은 사실을
    # 적는 자리는 SI 이고 거기는 'All six cells clear $q = 0.05$' 로
    # 띄어쓰기가 다르다.
    ('pn_cells', [r'[Aa]ll ([\w-]+) cells? clear \$q\s*=\s*0\.05\$'], 0),
    # 2026-11-04 (E2). 루이사이트 기하이성질체 민감도.
    ('lew_dnorm_E', [r'primary panel from (\d\.\d+) to'], 0.00005),
    ('lew_dnorm_Z', [r'from \d\.\d+ to (\d\.\d+), leaves CEES'], 0.00005),
    # 2026-11-03 (작업 C). 바닥 논증 대체 문단.
    ('pn_null_max', [r'difference of \$\+[\d.]+\$ to\s*\n?\$\+(0\.\d+)\$',
                     r'between \$\+[\d.]+\$ and\s*\n?\$\+(0\.\d+)\$'],
     0.00005),
    ('pn_null_min', [r'difference of \$\+(0\.\d+)\$ to',
                     r'between \$\+(0\.\d+)\$ and'], 0.00005),
    ('pn_obs_min', [r'observed differences of \$\+(0\.\d+)\$',
                    r'observed values of \$\+(0\.\d+)\$'], 0.00005),
    ('pn_obs_max', [r'observed differences of \$\+[\d.]+\$ to\s*\n?'
                    r'\$\+(0\.\d+)\$',
                    r'observed values of \$\+[\d.]+\$ to \$\+(0\.\d+)\$'],
     0.00005),
    ('conc_norm_cells',
     [r'\$\+1\.000\$ in every one of the ([a-z]+) run-channel cells'], 0),
    # WORDS 가 six -> 6 으로 바꿔 준다.
    ('pn_exact_n', [r'by enumeration of all ([\d,]+)\s*\n?relabellings'], 0),
    ('pn_within_min', [r'separates them, to \$\+(0\.\d+)\$'], 0.00005),
    ('pn_within_max', [r'to \$\+[\d.]+\$ to\s*\n?\$\+(0\.\d+)\$ within'],
     0.00005),
    ('pn_between_hi', [r'within compounds against \$(-0\.\d+)\$'], 0.00005),
    ('pn_between_lo', [r'against \$-[\d.]+\$ to\s*\n?\$(-0\.\d+)\$\s*\n?'
                       r'between them'], 0.00005),
    ('pn_zs_min', [r'between (0\.\d+) and 0\.\d+ after the'], 0.00005),
    ('pn_zs_max', [r'between 0\.\d+ and (0\.\d+) after the'], 0.00005),
    # MDE 의 BH18 문턱.
    ('mde_bh18_35', [r'moves the threshold to\s*\n?\$\\delta = (0\.\d+)\$'],
     0.005),
    ('mde_bh18_28', [r'corrected\s*\n?threshold is \$\\delta = (0\.\d+)\$'],
     0.005),
    ('mde_ratio_bh35', [r'then (\d\.\d+) times what the design'], 0.005),
    ('mde_ratio_bh28', [r'a ratio of (0\.\d+)'], 0.005),
    ('mde_gap_median_28', [r'a median gap of (0\.\d+)'], 0.0005),
    # 재표집 회복률. 41.0 을 DMMP 자리에 쓴 오류의 재발을 막는다.
    ('nb_nboot', [r'agent in [\d.]+\\% of (\d+) bootstrap replicates'], 0),
    ('nb_dmmp', [r'recovers sarin as its nearest\s*\n?agent in (\d+\.\d+)\\%'],
     0.05),
    ('nb_dpgme', [r'DPGME recovers sulfur mustard in (\d+\.\d+)\\%'], 0.05),
    ('nb_dcm', [r'recovers cyanogen chloride in (\d+\.\d+)\\%'], 0.05),
    # 2026-11-03 (작업 A-3) 폐기. chi2 와 4셀 p 는 철회된 합성의 값이다.
    # 교체된 자리는 본문 재집필(작업 C) 때 새 앵커를 붙인다.
    # 2026-11-03 (작업 A-3). 네 셀을 채널별로 묶어 적도록 고쳤으므로
    # 나열 순서가 d1_uv, d2_uv, d1_normal, d2_normal 로 바뀌었다.
    ('vol_rho_d1_uv', [r'four: \$\+(0\.\d+)\$ and'], 0.0005),
    ('vol_rho_d2_uv',
     [r'\$ and \$\+(0\.\d+)\$ in the ultraviolet'], 0.0005),
    ('vol_rho_d1_normal',
     [r'ultraviolet channel, \$\+(0\.\d+)\$ and'], 0.0005),
    # 2026-11-06. 'in the reflectance channel' 만으로는 DMMP 대비 문장도
    # 걸린다. 휘발성 절의 짝 문장으로 한정한다.
    ('vol_rho_d2_normal',
     [r'\$\+0\.\d+\$ and \$\+(0\.\d+)\$\s*\n?in the reflectance channel'], 0.0005),
    # 2026-11-03 (작업 A-3) 신설. 채널내 런간 Fisher + 채널 BH.
    ('vol_fisher_nor', [r'reflectance channel gives\s+\$p = (0\.\d+)\$'],
     0.00005),
    ('vol_fisher_uv', [r'ultraviolet channel does not\s*\n?\(\$q = (0\.\d+)\$\)',
                       r'ultraviolet channel gives \$p = (0\.\d+)\$'],
     0.00005),
    ('vol_fisher_minq', [r'reflectance\s+channel reaches \$q = (0\.\d+)\$',
                         r'\$p = 0\.0215\$ and \$q = (0\.\d+)\$',
                         # 2026-11-03 (작업 B). 초록 편입분.
                         r'reflectance\s*\n?channel at \$q=(0\.\d+)\$'],
     0.00005),
    # 2026-11-07. 본문 자리 전용 앵커. SI 자리는 vol_fisher_minq /
    # vol_p_nvalid 가 따로 지킨다.
    ('vol_fisher_minq_main',
     [r'reaches\s*\n?\$q = (0\.\d+)\$, in the\s*\n?reflectance channel'],
     0.00005),
    ('vol_p_nvalid_main',
     [r'valid-match count returns one match at \$p = (0\.\d+)\$'], 0.00005),
    ('n_freedoms',
     [r'(\w+) analyst choices therefore move',
      r'the (\w+) analyst freedoms in full'], 0),
    ('match_range_hi',
     [r'arbitrary choices move between zero and\s*\n?(\w+)'], 0),
    ('vol_gap_cg', [r'sits (0\.\d+) decades from phosgene',
                    r'(0\.\d+) decades from phosgene'], 0.0005),
    ('vol_gap_ac', [r'(0\.\d+) decades\s+from hydrogen cyanide',
                    r'but\s+(0\.\d+) from hydrogen cyanide',
                    r'sits (0\.\d+) decades\s*from hydrogen cyanide'],
     0.0005),
    ('vol_gap_gb', [r'sits (0\.\d+) decades from sarin',
                    r'sits\s+(0\.\d+) from sarin'], 0.0005),
    ('vol_gap_gb_pool', [r'an interferent, sits (0\.\d+)',
                         r'while an interferent sits (0\.\d+)'], 0.0005),
    ('vol_gap_hd', [r'(0\.\d+) decades from sulfur\s+mustard'], 0.0005),
    ('vol_gap_vx', [r'within (\d+\.\d+) decades of VX',
                    r'nothing within (\d+\.\d+) decades'], 0.0005),
    ('vol_n_int_nearest', [r'For (\w+) of the\s+eleven agents the nearest'],
     0),
    # ------ 2026-11-03 (2차). 리뷰어 1 이 적발한 분모·범위 오류 계열.
    ('dmmp_g4mean_uv', [r'over all four G agents is \$\+([\d.]+)\$ in the'],
     0.0005),
    ('dmmp_g4mean_nor',
     [r'ultraviolet channel and \$(-[\d.]+)\$\s*\n?in reflectance'], 0.0005),
    # 2026-11-07 (12차 7단계). 분량 조정에서 'it is' 를 뺐다.
    ('dmmp_gmean_uv_3',
     [r'other than tabun(?: it is)? \$(-[\d.]+)\$'], 0.0005),
    ('dmmp_gmean_nor_3',
     [r'other than tabun(?: it is)? \$-[\d.]+\$ and\s*\n?\$(-[\d.]+)\$'],
     0.0005),
    ('vol_lew_cis_p',
     [r'from \$p = [\d.]+\$ at rank\s*\n?fourteen to \$p = ([\d.]+)\$'],
     0.00005),
    ('vol_lew_cis_rank', [r'at rank\s*\n?fourteen to \$p = [\d.]+\$ at rank '
                          r'(\w+)\.'], 0),
    # ------ 2026-11-03 (2차). 교정한 자리에 앵커가 없어 주입이 통과했다.
    #        주입 시험이 그것을 잡아냈고, 여기서 메운다.
    ('rand3_frac28_2', [r'rate falls from ([\d.]+)\\% to'], 0.05),
    ('rand3_med28_2', [r'median \$p\$ rises from ([\d.]+) to'], 0.0005),
    ('white_med28_2', [r'median \$p\$ rises from [\d.]+ to ([\d.]+)'],
     0.0005),
    ('pct_struct_2', [r'the ([\d.]+)th percentile of the subspaces one'], 0.06),
    ('fam_bh18_p28_2',
     [r'eighteen tests whose smallest adjusted value is ([\d.]+)'], 0.0006),
    ('pairs_d2uv_struct', [r'structural cell at \$p = ([\d.]+)\$ rather'],
     0.0006),
    ('mantel_d2uv_struct',
     [r'structural cell at \$p = [\d.]+\$ rather than ([\d.]+)'], 0.0006),
    ('pairs_d2uv_merged', [r'the merged cell at ([\d.]+) rather than'], 0.0006),
    ('mantel_d2uv_merged',
     [r'the merged cell at [\d.]+ rather than ([\d.]+)'], 0.0006),
    ('pairs_d2nor_elec',
     [r'electronic cell it runs the other way, to ([\d.]+) from'], 0.0006),
    ('mantel_d2nor_elec_2',
     [r'runs the other way, to [\d.]+ from ([\d.]+)'], 0.0006),
    ('pooled_struct_uv', [r'gives \$\\rho = \+([\d.]+)\$ \(95'], 0.0006),
    ('pooled_struct_uv_lo', [r'\(95\\% CI \$(-[\d.]+)\$ to \$\+[\d.]+\$\) in '
                             r'the ultraviolet'], 0.0006),
    ('pooled_struct_uv_hi', [r'CI \$-[\d.]+\$ to \$\+([\d.]+)\$\) in the '
                             r'ultraviolet'], 0.0006),
    ('pooled_struct_nor', [r'and \$\+([\d.]+)\$ \(\$-[\d.]+\$ to \$\+[\d.]+\$\)'
                           r' in reflectance'], 0.0006),
    ('struct_kept_28_2', [r'The structural space has (\d+) coordinates'], 0),
]

WORDS = {'zero': 0, 'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6,
         'seven': 7, 'eight': 8, 'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12,
         'thirteen': 13, 'fourteen': 14, 'fifteen': 15, 'sixteen': 16,
         'seventeen': 17, 'eighteen': 18, 'nineteen': 19, 'twenty': 20,
         'twenty-one': 21, 'twenty-two': 22,
         'twenty-six': 26, 'twenty-seven': 27, 'twenty-eight': 28,
         'twenty-nine': 29, 'thirty': 30,
         'twenty-three': 23, 'twenty-four': 24, 'twenty-five': 25,
         'thirty-three': 33, 'thirty-four': 34,
         # 2026-11-03 (작업 B). 병합공간 차원을 본문이 영문으로 적는다.
         'thirty-one': 31, 'thirty-two': 32, 'thirty-five': 35,
         'thirty-six': 36, 'thirty-seven': 37, 'thirty-eight': 38,
         'thirty-nine': 39, 'forty': 40, 'forty-one': 41}

# 서수. 2026-11-03. 휘발성 절이 '순위 열네 번째'를 영문 서수로 적는다.
# 서수를 숫자로 바꿔 쓰면 영어가 어색해지므로 맵을 넓힌다.
WORDS.update({
    'first': 1, 'second': 2, 'third': 3, 'fourth': 4, 'fifth': 5,
    'sixth': 6, 'seventh': 7, 'eighth': 8, 'ninth': 9, 'tenth': 10,
    'eleventh': 11, 'twelfth': 12, 'thirteenth': 13, 'fourteenth': 14,
    'fifteenth': 15, 'sixteenth': 16, 'seventeenth': 17,
    'eighteenth': 18, 'nineteenth': 19, 'twentieth': 20,
    'twenty-first': 21, 'twenty-second': 22, 'twenty-third': 23,
    'twenty-fourth': 24, 'twenty-fifth': 25, 'twenty-sixth': 26,
    'twenty-seventh': 27, 'twenty-eighth': 28})


# 폐기값 블랙리스트. 재분석으로 값이 바뀔 때마다 구판을 여기에 넣는다.
# 본문이 구판임을 밝히고 인용하는 문장은 면제 패턴으로 통과시킨다.
SUPERSEDED = [
    # 2026-11-04 (작업 D). 패널 구성 명명. 8 은 문헌이 제안한 화합물 수가
    # 아니라 지정 후보의 총수다 — 기존 3종(DMMP, DPGME, 디클로로메탄) +
    # 문헌 5종(CEES, DCNP, DCP, DFP, 살리실산메틸). 'the eight literature
    # simulants' 는 문헌이 8종을 제안했다고 읽힌다. 생성표 라벨도
    # '여덟 지정' 으로 통일했다.
    # 2026-11-04 (E2). 이성질체 교체는 판정을 바꾸지 않는다(둘 다 음성).
    # 바뀐다고 적으면 민감도의 결론이 뒤집힌다.
    ('이성질체가 매칭 판정을 바꾼다는 주장 (실측: E·Z 모두 음성)',
     [r'changes the match verdict, which is positive',
      r'the \$Z\$ form makes it a match'],
     []),
    # 2026-11-04 (E1). 런 1·3 은 어떤 순도 수치도 보고하지 않는다. 원문을
    # 읽어 확인했다. 공급자 등급·역가·확인분석 어느 것도 없다.
    ('런 1·3 에 순도 수치가 있다는 주장 (원문 확인: 전무)',
     [r'Both report supplier grades',
      r'first and third runs report purit',
      r'purity (?:was )?confirmed for all three runs'],
     []),
    # 2026-11-04 (E1/E6). 원논문이 말하지 않는 출처 주장. 런 1 은 작용제를
    # 연구소에서 '받았다'(obtained from)고만 적고, 런 3 은 출처를 아예 적지
    # 않는다. 어느 쪽도 '생산했다' 고 말하지 않는다.
    ('작용제 생산 주장 (원논문 미기재 — 런1 obtained from, 런3 무기재)',
     [r'All agents were produced\s*\n?and handled at',
      r'agents were produced and handled'],
     []),
    # 2026-11-04 (E3). 원논문은 A계열 코드가 어느 화합물인지 동정하지 않고
    # 구조 그림을 'representative' 로 단다. 우리가 구조를 배정하지 않는다.
    ('A계열 코드에 구조·화합물 배정 (원논문이 동정하지 않는다)',
     [r'A-230 is\b', r'A-232 is\b', r'A-234 is\b', r'A-242 is\b',
      r'A-series compounds are identified as'],
     []),
    # 2026-11-04 (G3). S-disc 의 두 수는 서로 다른 것을 센다. 207종 중
    # 일곱만 DCM 보다 안전하고(경계값 DCM 자신은 분면 8에 포함), 한 자릿수
    # 아래 완화 기준에서 후보가 여덟이다. 'the seven are a lower bound' 는
    # 바로 앞 문장의 여덟과 뭉개져 읽힌다.
    ('S-disc 의 두 수 지시대상 뭉갬 (일곱=더 안전, 여덟=완화기준 후보)',
     [r'the seven are a lower bound',
      r'safer than all but seven of the'],
     []),
    ('패널 구성 명명 구판 (8 = 기존 3 + 문헌 5, 문헌이 8종이 아니다)',
     [r'eight literature simulants',
      r'literature eight',
      r'\+ literature 8'],
     []),
    # 2026-11-03 (작업 C). 중심화 '바닥' 논증 전량 철회.
    #
    # 두 가지가 겹쳐 있었다. (1) centred() 가 unit(v-mu) 로 재정규화하므로
    # Sum u = 0 이 깨진다. 실측: 평균을 뺀 직후 7.2e-16, 정규화 후 0.1292
    # ~0.7430 (premise_null.py zero_sum_residual). 항등식의 전제가 없다.
    # (2) 영합 벡터에서도 -1/(n-1) 은 노름가중 평균 코사인의 상한이고
    # 등노름에서만 등호다. '바닥' 이라는 이름부터 방향이 반대였다.
    # 대체: 농도 내 라벨 순열 귀무. 여섯 셀 전부 유의, BH 최소 q 0.0001.
    ('중심화 바닥 논증 (항등식 전제 붕괴 + 부등호 방향 반대)',
     [r'algebraic floor',
      r'clear that floor by',
      r'clears? that floor',
      r'mean cosine of at least',
      r'fixes the mean cosine between',
      r'floor its compound count\s+allows',
      r'above the floor, which is close enough'],
     []),
    # 부호뒤집기 전수열거는 여섯 셀 전부에서 자기 설계 바닥을 돌려준다.
    # p 가 데이터가 아니라 설계의 함수이므로 유효한 검정으로 쓰지 않는다.
    ('부호뒤집기 검정을 증거로 사용 (설계 바닥을 돌려준다)',
     [r'cannot reach 0\.05 under any test',
      r'and not\s+a loss of evidence',
      r'every cell returns the smallest value its compound count allows:'],
     []),
    # 2026-11-03 (작업 B). 분석선택 수와 매칭 범위가 구판으로 되돌아가는 것을
    # 막는다. 형태 블록이 다섯째 선택이고 그 선택에서 28종 패널 매칭이 0 이
    # 되므로 범위 하단이 셋에서 영으로 내려갔다. 검증기 사실이 아니라
    # 열거의 결과이므로 금지항목으로 지킨다.
    ('분석선택 수 구판 (정본 다섯. 형태 블록이 다섯째다)',
     [r'Four analyst choices therefore move',
      r'moved by all four and is not a measurement',
      r'four arbitrary choices can move'],
     []),
    ('유효매칭 범위 구판 (정본 영에서 여섯)',
     [r'can move between\s+three and six'],
     []),
    # 2026-11-03 (작업 B). 제목의 Sixteen. A계열 4종은 어느 분석에도
    # 기여하지 않고 1차 패널은 11종, 중심 검정은 셀당 4~7종이다. 제목에
    # 들어갈 수 있는 수가 없으므로 수를 적지 않는다.
    ('제목의 작용제 수 (B3. 측정 16 대 1차패널 11 대 셀당 4~7)',
     [r'Reanalysis of Sixteen Authentic'],
     []),
    # 2026-11-03 (2차). 철회한 두 주장. 되살아나면 멈춘다.
    #
    # (a) DPGME 의 지정 계열. 등재부의 'blister' 는 PMC3231691 본문에서
    #     온 라벨이고 지정 문서가 아니다(2026-09-03 기록, 승인 미회신).
    #     SI 가 특허 문언이라고 적었던 'a preferred G agent simulant' 는
    #     특허 원문을 확인할 수 없었다. 두 주장 모두 쓰지 않는다. 본문이
    #     사실로 적는 것은 '특허가 요구한 성질은 IMS 검출성' 하나다.
    ('DPGME 지정 계열 주장 (특허 원문 미확인, 등재부 라벨은 지정문서 아님)',
     [r'a preferred G agent simulant',
      r'blister-agent simulant with a nerve agent',
      r'designated for the blister agent',
      r'designated for the choking agent'],
     [r'An\s+earlier version of this section wrote',
      r'We withdraw that\s+statement']),
    # 2026-11-02 신설. 세 원논문의 Methods 를 직접 읽기 전에, 어레이 농도를
    # ppm 으로, 같은 화합물의 서로 다른 농도를 '반복측정(replicate)' 으로
    # 적었다. 둘 다 수행하지 않은 측정을 서술한 것이므로 [C-3] 위반이다.
    # 정본: 농도는 아세토니트릴 용액 µM, 반복측정은 없다.
    ('어레이 농도 단위, ppm 오기 (정본 µM)',
     [r'(?:\b10|\b50|\b100|\b500)~?\\?,?\s*ppm\b'],
     [r'AEGL', r'mg\\,m\^\{-3\}', r'converted using', r'thousands separator',
      r'bracketed by the agency']),
    ('어레이 반복측정 주장 (정본: 조합당 1행, 반복측정 없음)',
     [r'[Rr]eplicates? of (?:one|the same) compound',
      r'within replicates of',
      r'replicate agreement exceeds',
      # 긍정 주장 꼴. 정본은 'not agreement between repeated measurements'
      # 라는 부정문이므로 아래 면제 패턴으로 통과시킨다. 2026-11-02 에
      # 주입시험이 이 꼴을 통과시켜 가드의 구멍이 드러났다.
      # 2026-11-07 (11차 주입시험). 종전 꼴은 'between repeated measurements'
      # 전부를 금지하고, 정본의 부정문을 면제 패턴 두 개로 통과시켰다.
      # 면제 창이 ±240자라서 같은 단락 안의 주입('that is between repeated
      # measurements')도 함께 면제됐다. 부정어를 패턴 자체에 넣어 정본이
      # 애초에 걸리지 않게 하고, 단락 면제를 없앤다.
      r'(?<!not )between repeated\s+measurements'],
     [r'descriptor[- ]resampling', r'resampl']),
    # 2026-11-02 (T2-1). 이 항목은 구판 정본이 0.0265 (3자리로 0.027) 일 때
    # 0.028 을 금지했다. 결정적 상관필터로 교체한 뒤 정본이 0.0279 (3자리로
    # 0.028) 가 되었으므로, 이 항목을 그대로 두면 블랙리스트가 정본 값을
    # 금지한다. 아래 _blacklist_staleness() 가 이 상황을 자동으로 잡는다.
    # 금지 대상은 그 이전 판본인 0.0265 쪽으로 옮긴다.
    ('28종 구조공간 명목 p, 결정적 필터 이전 (정본 0.028)',
     [r'structural[\s\S]{0,60}0\.0265\b', r'0\.0265[\s\S]{0,60}percentile'],
     [r'earlier version', r'superseded', r'single-conformer',
      r'column order', r'input order']),
    ('부분공간 백분위 (정본 22.4)',
     [r'25\.8th percentile', r'25\.4th percentile', r'\b28th percentile',
      r'22\.7th percentile', r'\b23rd percentile', r'\b25th percentile',
      r'24\.3rd percentile', r'25\.9th percentile'], []),
    # 2026-11-02. 종전 패턴은 벗은 0.031 이었다. 병리지수 제거 변종의 정본
    # 구조공간 p 가 0.0309 이고 3자리로 0.031 이므로, 그 값을 적을 자리마다
    # 실패가 났다. 28종 본 분석을 가리키는 어법으로 좁힌다.
    ('28종 구조공간 명목 p, 앙상블 이전 (정본 0.028)',
     [r'structural[\s\S]{0,60}0\.031\b', r'0\.031[\s\S]{0,40}one-sided'],
     [r'pathological', r'without the two']),
    ('28종 구조좌표 수, 메타열 제외 이전 (정본 36)',
     [r'37-dimensional', r'Weighting the 37 structural',
      r'leaving \d+ and 37 structural'], []),
    ('전자공간 백분위 (정본 76.9)', [r'74\.2nd'], []),
    ('Li-Ji 유효검정수 (정본 5.00)', [r'6\.00 effective tests'], []),
    ('도판2 캡션의 통계량 서술 (실제로는 평균차 하나)',
     [r'Smallest \$p\$ over\s+the six statistics for 2000'], []),
    ('부분공간 0.05 미만 비율 (정본 36.0 / 백색화 16.5)',
     [r'36\.3\\?%', r'32\.6\\?%', r'33\\% of them fall below',
      r'19\\% after whitening', r'18\.9\\?%'], []),
    # 2026-11-02. 정본이 0.0419 (3자리로 0.042) 가 되었으므로 0.042 를
    # 금지할 수 없다. 금지 대상을 그 앞 두 판본으로 옮긴다.
    ('사전지정 BH q, 결정적 필터 이전 (정본 0.042)',
     [r'0\.046\b', r'0\.0397\b'], []),
    ('28종 BH over 18, 결정적 필터 이전 (정본 0.163)',
     [r'adjusted value is 0\.159', r'eighteen tests is 0\.159',
      r'over 18\s+& [\d.]+ & [\d.]+ & 0\.159',
      r'adjusted value is 0\.1809', r'gives 0\.1809'], []),
    ('28종 Westfall-Young (정본 0.079)',
     [r'approximating it, gives 0\.075',
      r'Westfall[\s\S]{0,160}0\.075', r'0\.075[\s\S]{0,160}Westfall'], []),
    # 2026-11-02 (T2-1). 정본 매치 수가 4 에서 3 으로 내려갔다. 열 순서를
    # 40회 치환하면 구판 필터가 3건 25회 / 4건 15회를 내므로 '4' 가 소수
    # 결과였다. 이제 금지해야 하는 쪽은 '네 건' 서술이다.
    ('28종 매치 수, 결정적 필터 이전 (정본 3)',
     [r'restores four matches', r'counted four matches',
      r'four agents are matched', r'supplies all four',
      r'[Ff]our agents acquire a valid match'],
     [r'earlier version', r'superseded', r'column order', r'input order']),
    ('Mantel 다중성 (정본 12셀 0.117)',
     [r'0\.177', r'adjusted value from [\s\S]{0,30}to\s*0\.124',
      r'adjusted value is 0\.100\b'], []),
    ('B3LYP 매칭 4건을 전부 문헌 주장 짝이라고 서술',
     [r'the four are the pairings the sensor',
      r'to DPGME, which are the pairings'], []),
    # 2026-10-27. 쌍극자를 GFN2 로 남겼을 때의 비유의 서술. 세 좌표 전부
    # B3LYP 이면 p = 0.0368 이므로 이 문장은 데이터와 반대다.
    ('B3LYP 전자공간 비유의 주장 (정본 p=0.0368)',
     [r'permutation test is still not\s+significant at any conventional',
      r'not significant at any conventional level',
      r'not evidence that the\s+incumbents are close to the agents'], []),
    ('B3LYP 3좌표 매칭 5건을 전부 기록 지지 짝이라고 서술 (정본 4/5)',
     [r'[Tt]he five are pairings the designation',
      r'all five are pairings', r'five of them pairings the designation'], []),
    ('28종 구조공간 유의 통계량 수 (정본 2)', [r'one of six statistics in the structural'], []),
    ('28종 구조좌표 수 (정본 36)',
     [r'38-dimensional value', r'thirty-eight structural coordinates are not'], []),
    # 2026-11-02 (T2-1). 종전 항목은 정본이 1.989 일 때 2.003 을 금지했다.
    # 결정적 필터에서 정본이 2.0038 이 되었고, 금지 패턴 '2\.003' 이 그
    # 부분문자열에 걸려 옳은 값이 실패를 냈다. B2a 의 반올림 대조로는
    # 이 꼴(부분문자열 충돌)을 잡을 수 없으므로 여기에 기록해 둔다.
    # 금지 대상은 그 이전 두 판본이다. 구판을 명시해 인용하는 자리는
    # 면제한다. 열 순서 절이 1.9948 을 대조값으로 인용한다.
    ('tabun 정규화거리, 결정적 필터 이전 (정본 2.0038)',
     [r'\b1\.989\b', r'tabun[\s\S]{0,60}\b1\.995\b',
      r'\b1\.995\b[\s\S]{0,60}threshold'],
     [r'earlier', r'superseded', r'put tabun at', r'old filter',
      r'previous filter', r'column order']),
    ('두 패널 감소가 같다는 서술 (정본: 13 대 15)',
     [r'attrition is the same on\s+both panels'], []),
]


def tonum(s):
    s = s.strip()
    if s.lower() in WORDS:
        return float(WORDS[s.lower()])
    # LaTeX 천단위 구분자(\,)와 쉼표를 모두 제거한다. 순서를 뒤집으면
    # '41\\,685' 가 '41\\685' 로 남아 파싱이 깨진다.
    # 2026-11-01. 캡처가 문장 끝 마침표를 함께 먹으면 ValueError 로 검증기
    # 전체가 죽었다. 숫자 뒤 마침표 하나는 문장부호로 보고 떼어낸다.
    s = re.sub(r'(?<=\d)\.$', '', s)
    return float(s.replace('\\,', '').replace(',', ''))


def main():
    if not preflight():
        return 1
    f = load_facts()
    # 문헌 인용 수치를 같은 사실표에 올린다. 키 충돌은 즉시 멈춘다 (산출물
    # 값이 문헌값으로 조용히 덮이면 검사 전체가 무의미해진다).
    for _k, (_v, _src, _quote) in LITERATURE.items():
        if _k in f:
            raise SystemExit(f'[중단] 문헌 키 {_k} 가 산출물 키와 충돌한다')
        f[_k] = _v
    # 2026-11-02. 같은 양을 두 군데서 센다. (B) 검사는 어법당 한 번만
    # 맞추므로, 2.4 의 어법이 맞으면 Conclusions 의 다른 어법은 보이지
    # 않는다. 주입시험에서 실제로 Conclusions 의 BH 값만 바꿔도 통과했다.
    # 자리마다 키를 따로 두고 같은 산출물 값에 건다.
    # 2026-11-02. 커버레터도 본문 수치를 다시 센다. 주입시험에서 0.964,
    # +0.76, +0.77 세 자리가 전부 통과했다 (어법이 본문과 달라 (B) 검사가
    # 본문 쪽에서 이미 맞고 끝났다). 자리마다 키를 따로 둔다.
    for _dst, _src_key in (('bh_q35_concl', 'bh_q35'),
                           ('cl_dmmp_ga', 'dmmp_ga_uv'),
                           ('cl_dmmp_gb', 'dmmp_gb_uv'),
                           ('cl_dmmp_gd', 'dmmp_gd_uv'),
                           ('cl_dmmp_gf', 'dmmp_gf_uv'),
                           ('cl_dmmp_vx', 'dmmp_vx_uv'),
                           # 2026-11-02 (GATE 8 7차). ACS 본문 §2.4 가 SI 의
                           # DFT 수치를 다시 센다. SI 어법에 앵커가 이미
                           # 걸려 있어 본문 사본은 무검사로 남았다.
                           ('mt_dft_gfn2_p', 'dft_es_gfn2_p'),
                           ('mt_dft_full_p', 'dft_es_full_p'),
                           ('mt_dft_full_match', 'dft_es_full_match'),
                           ('mt_dft_mg_p', 'dft_mg_b3lyp_p'),
                           ('mt_dft_mg_match', 'dft_mg_b3lyp_match'),
                           ('mt_dft_mg_gfn2_p', 'dft_mg_gfn2_p')):
        if _src_key in f:
            f[_dst] = f[_src_key]
    def strip_comments(s):
        # LaTeX 주석은 PDF 에 들어가지 않는다. 주석을 읽으면 '보관용' 문구가
        # 검증 통과의 그림자 출처가 되므로 대조 전에 제거한다. \% 는 보존.
        out = []
        for line in s.split('\n'):
            if line.lstrip().startswith('%'):
                continue
            out.append(re.sub(r'(?<!\\)%.*$', '', line))
        return '\n'.join(out)

    tex = strip_comments(open(TEX).read())
    body = tex[tex.index('\\begin{abstract}'):]
    # 2026-10-13. 종전에는 main.tex 만 읽었다. 그래서 커버레터와 하이라이트가
    # 본문과 다른 Spearman·후보수를 싣고도 PASS 가 나왔다. 투고 패키지에 들어가는
    # 글은 전부 같은 산출물과 대조한다.
    # 2026-11-01. 목록이 아래 정규화본 쪽에도 똑같이 적혀 있어 한쪽만 고치면
    # 조용히 어긋났다. 상수 하나로 묶는다. paper_acs/main.tex 를 넣는 이유는
    # ACS Sensors 재집필본이 지금부터 커지는데, 검증기가 elsarticle 판만 보면
    # 새 원고가 무검사 상태로 자라기 때문이다. 이 세션에서 반복해 겪은
    # '한쪽 문서만 구판' 사고를 미리 막는다.
    for extra in EXTRA_TEX:
        if os.path.exists(extra):
            body += '\n' + strip_comments(open(extra).read())
    # 도판 주석도 같은 규칙으로 검사한다. 본문만 보던 설계 때문에 Figure 2 의
    # 구판 문자열('the least acutely toxic of the thirteen')이 통과했다.
    import fig_annot
    annot = fig_annot.all_text()
    if not annot:
        print('[verify] 경고: outputs/figure_annotations.json 이 비어 있다. '
              '도판 스크립트를 먼저 실행해야 주석이 검사된다.')
    body = body + '\n' + annot
    # 공백 정규화본. 문서와 도판 주석 레코드 사이에는 \x00 을 넣는다. 패턴의
    # \s 는 이 문자를 먹지 않으므로 경계를 넘는 오탐이 생기지 않는다. 넣지
    # 않았더니 주석 두 레코드에 걸친 문자열이 한 칸으로 읽혀 세 건이 잘못
    # 걸렸다 (2026-10-24).
    _segs = [re.sub(r'\s+', ' ', seg) for seg in
             ([tex] + [strip_comments(open(x).read()) for x in EXTRA_TEX
                       if os.path.exists(x)]
              )]
    # 도판 주석은 레코드가 이미 한 줄이라 줄바꿈 문제가 없다. 정규화본에
    # 넣으면 레코드끼리 붙어 오탐만 생기므로 원문 쪽에만 남긴다.
    body_flat = '\x00'.join(_segs)
    fails, checked, unseen = [], 0, []
    for _g, _why in ((_ordinal_typos(), '서수 접미사 오류'),
                     (_handmade_tables(), '생성기 없는 수치 표'),
                     (_generated_table_sync(), '생성 표 사본 시차'),
                     (_title_sync(), '제목 비동기'),
                     (_figure_panel_refs(), '없는 도판 패널 인용'),
                     (_si_section_labels(), 'SI 항해성'),
                     (_filter_dependent_staleness(), '구판 예치물'),
                     (_label_scope_guards(), '수치는 맞고 설명이 틀린 자리')):
        if _g:
            fails.append(f'{_why} (2026-11-03 신설 가드):')
            fails.extend(_g)
    _dp_hits = _dropped_periods()
    if _dp_hits:
        fails.append('문장 끝 마침표가 빠진 자리 (되쓰기 사고 가드):')
        fails.extend(_dp_hits)

    # (D) 도판 신선도. 2026-10-13 신설. 도판은 outputs/ 에 생성되고 paper/ 로
    # 복사되는데 그 복사가 수작업이었다. 분석을 다시 돌린 뒤 복사를 잊으면
    # 원고에 옛 그림이 박히고, 수치 대조는 전부 통과한다. GATE 8 3차에서
    # 도판과 본문이 어긋난 경로가 이것이므로 시각을 직접 비교한다.
    FIGS = ['Fig_main.png', 'Fig_plocal.png', 'Fig_h1.png', 'Fig_mde.png',
            'Fig_talk_response.png', 'Fig_talk_tox.png']
    for fig in FIGS:
        # 2026-11-04 (J-C). 사본 대상을 paper_acs/ 로 옮긴다.
        src, dst = f'outputs/{fig}', f'paper_acs/{fig}'
        if not os.path.exists(src):
            unseen.append(f'{fig}: outputs/ 에 없음 (도판 스크립트 미실행)')
            continue
        if not os.path.exists(dst):
            fails.append(f'  도판 누락 {fig}: paper_acs/ 에 없다')
            continue
        if os.path.getmtime(dst) < os.path.getmtime(src) - MTIME_TOL:
            fails.append(f'  도판 구판 {fig}: paper_acs/ 사본이 outputs/ 보다 '
                         f'오래됐다 (sh build_paper_acs.sh 로 복사할 것)')
        checked += 1
    # 분석 산출물보다 도판이 오래되면 그림이 낡은 수치를 그리고 있다.
    newest = max(os.path.getmtime(x) for x in
                 ['outputs/experiment_summary.json', 'outputs/rollback_c.json']
                 if os.path.exists(x))
    for fig in FIGS:
        src = f'outputs/{fig}'
        if os.path.exists(src) and os.path.getmtime(src) < newest - MTIME_TOL:
            fails.append(f'  도판 구판 {fig}: 분석 산출물보다 오래됐다 '
                         f'(sh build_paper_acs.sh --figs 로 재생성할 것)')

    for key, pats, tol in CHECKS:
        if key not in f:
            unseen.append(f'{key}: 사실표에 없음')
            continue
        want = float(f[key])
        found = []
        for pat in pats:
            for m in re.finditer(pat, body):
                found.append((tonum(m.group(1)), m.group(0)[:62].replace('\n', ' '), m.group(1)))
        # 2026-10-24. 패턴이 리터럴 공백 하나를 요구하는데 그 자리에서 TeX 줄이
        # 바뀌면 문맥을 놓친다. 처음에는 '못 찾았을 때만' 공백을 눌러 다시
        # 보게 했는데, 그러면 같은 수치가 두 군데 적혀 있고 한쪽만 줄바꿈에
        # 걸린 경우 그 한쪽이 검사되지 않는다. 실제로 본문의 백분위 두 개가
        # 그 경로로 구판인 채 통과했다. 항상 둘 다 훑고 중복만 걷어낸다.
        for pat in pats:
            for m in re.finditer(pat, body_flat):
                found.append((tonum(m.group(1)), m.group(0)[:62], m.group(1)))
        seen, uniq = set(), []
        for val, ctx, raw in found:
            k = (round(val, 12), ' '.join(ctx.split()))
            if k in seen:
                continue
            seen.add(k)
            uniq.append((val, ctx, raw))
        found = uniq
        if not found:
            unseen.append(f'{key}: 본문에서 문맥을 찾지 못함 (값 {want:g})')
            continue
        for got, ctx, raw in found:
            checked += 1
            # 2026-10-27. 절대 허용오차만 쓰면 반올림 오류가 빠져나간다.
            # 0.0264706 을 본문이 0.027 로 적었는데 |차| 0.00053 이 허용오차
            # 0.0006 보다 작아 통과했다. 엄격한 검사(tol <= 0.0006)에서는
            # 본문이 적은 자릿수를 세어, 그 자릿수에서 정본의 유효한
            # 반올림인지를 묻는다. 반올림 방향(올림/내림)은 따지지 않는다.
            if tol <= 0.0006:
                _m = re.search(r'\.(\d+)', raw)
                _dp = len(_m.group(1)) if _m else 0
                _half = 0.5 * (10 ** -_dp) + 1e-12
                if abs(want - got) > _half:
                    fails.append(
                        f'  반올림 불일치 {key}: 본문 {raw} vs 산출물 '
                        f'{want:.6g} (소수 {_dp}자리면 {want:.{_dp}f})'
                        f'   [{ctx}]')
                continue
            if abs(got - want) > tol:
                fails.append(f'  불일치 {key}: 본문 {got:g} vs 산출물 {want:g}'
                             f'   [{ctx}]')

    # (B2) 폐기값 블랙리스트 (2026-10-25 신설).
    #
    # 왜 필요한가. (B) 의 검사는 어법 단위다. 'from $p=0.031$ to' 라는 패턴
    # 하나가 어딘가에서 맞으면 통과하므로, 같은 양을 다른 문장으로 적은
    # 자리에 구판 값이 남아 있어도 보이지 않는다. 2026-10-25 에 그 경로로
    # 본문 수치 16건이 PASS 를 받은 채 남아 있었다. 여기서는 방향을 뒤집어,
    # 폐기된 값이 원고 어디에든 나타나면 실패로 본다.
    #
    # 각 항목은 (사유, 금지 패턴들, 면제 패턴들)이다. 면제는 구판임을
    # 명시하고 인용하는 문장을 위한 것이다. 금지 패턴이 맞은 자리 앞뒤
    # 240자 안에서 면제 패턴이 맞으면 의도적 인용으로 본다.
    # (B2a) 블랙리스트 자체의 신선도 (2026-11-02 신설, T2-1 에서 드러남).
    #
    # 왜 필요한가. SUPERSEDED 는 '이 값은 폐기됐다' 는 금지 목록이다. 정본이
    # 바뀌면 어제 폐기값이던 수가 오늘 정본이 될 수 있다. 실제로 상관필터를
    # 결정적 규칙으로 교체하자 28종 구조공간 명목 p 의 정본이 0.027 에서
    # 0.028 로 올라갔는데, 블랙리스트는 0.028 을 금지하고 있었다. 그러면
    # 원고가 옳은 값을 적을 때마다 검증기가 실패를 낸다. 더 나쁜 것은
    # 분석자가 그 실패를 피하려고 틀린 값을 되돌려 적는 경로다.
    #
    # 그래서 금지 패턴 중 '벗은 소수 리터럴' 꼴인 것을 골라, 그 값이 현재
    # 사실표의 어떤 값과 그 자릿수에서 일치하는지 본다. 일치하면 그 항목은
    # 낡았다고 보고한다. 어법이 섞인 패턴은 판정하지 않는다.
    _bare = re.compile(r'^(?:\\b)?(\d+\.\d+)(?:\\b)?$')
    _stale = []
    for why, pats, exempt in SUPERSEDED:
        for pat in pats:
            _m = _bare.match(pat.replace('\\.', '.'))
            if not _m:
                continue
            _v = float(_m.group(1))
            _dp = len(_m.group(1).split('.')[1])
            for _k, _fv in f.items():
                try:
                    _fvf = float(_fv)
                except (TypeError, ValueError):
                    continue
                if abs(round(_fvf, _dp) - _v) < 1e-12:
                    _stale.append(f'  낡은 금지 항목 [{why}]: {pat} 가 현재 '
                                  f'정본 {_k} = {_fvf:.6g} 의 {_dp}자리 '
                                  f'반올림과 같다. 항목을 갱신할 것')
    if _stale:
        fails.append('블랙리스트가 정본 값을 금지한다 (검사 B2a):')
        fails.extend(sorted(set(_stale)))

    superseded = []
    for why, pats, exempt in SUPERSEDED:
        hits = []
        for src in (body, body_flat):
            for pat in pats:
                for m in re.finditer(pat, src):
                    a, b = max(0, m.start() - 240), m.end() + 240
                    win = src[a:b]
                    if any(re.search(e, win, re.I) for e in exempt):
                        continue
                    hits.append(' '.join(m.group(0).split())[:64])
        if hits:
            uniq = sorted(set(hits))
            superseded.append(f'  폐기값 잔존 {why}: ' + ' / '.join(uniq[:3]))
    if superseded:
        fails.append('폐기값이 원고에 남아 있다 (검사 B2):')
        fails.extend(superseded)

    # (A-pdf) 올라가는 PDF 가 원고보다 오래됐는지 (2026-11-07 신설).
    #
    # 왜 필요한가. build_paper_acs.sh 가 본문과 SI 만 빌드하고 커버레터는
    # 부르지 않았다. 그래서 paper_acs/cover_letter.pdf 가 나흘 전(11-03)
    # 판본으로 멈춰 있었고, 11차 교정 전체 — 제목 교체와 분량 사유 단락
    # 신설 — 가 빠진 PDF 가 투고될 상태였다. 검증기는 .tex 를 읽으므로
    # 글은 맞았고 올라가는 파일만 틀렸다. 도판에는 시차 검사가 있었는데
    # PDF 에는 없었다.
    # mtime 이 아니라 내용 해시로 본다. 주입시험은 원고를 스냅샷으로
    # 되돌리므로 mtime 이 새로워지지만 내용은 같다 (2026-11-07 실측).
    import hashlib as _hl
    import os as _o2
    _stamp = {}
    if _o2.path.exists('paper_acs/.buildstamp'):
        for _ln in open('paper_acs/.buildstamp'):
            _h, _, _n = _ln.strip().partition(' ')
            _stamp[_o2.path.basename(_n.strip())] = _h
    for _b in ('main', 'supporting', 'cover_letter'):
        _t, _p = f'paper_acs/{_b}.tex', f'paper_acs/{_b}.pdf'
        if not _o2.path.exists(_p):
            fails.append(f'{_p} 가 없다 (검사 A-pdf): '
                         'sh build_paper_acs.sh 를 돌릴 것')
            continue
        if not _stamp:
            fails.append('paper_acs/.buildstamp 가 없다 (검사 A-pdf): '
                         'sh build_paper_acs.sh 를 돌릴 것')
            break
        _cur = _hl.sha256(open(_t, 'rb').read()).hexdigest()
        if _stamp.get(f'{_b}.tex') != _cur:
            fails.append(f'{_p} 가 {_b}.tex 와 다른 판본이다 (검사 A-pdf): '
                         'sh build_paper_acs.sh 를 돌릴 것 — 올라가는 '
                         'PDF 가 원고와 다르다')

    # (A-tag) 원고가 가리키는 예치 태그가 실제로 공개돼 있는지
    # (2026-11-07 신설, 12차 A-1).
    #
    # 왜 필요한가. 11차에 Data Availability 를 GitHub 태그로 확정했는데,
    # 그 태그와 그것이 가리키는 커밋이 로컬에만 있었다. 즉 원고가 존재하지
    # 않는 공개 상태와 존재하지 않는 라이선스를 주장하는 채로 verify 500건을
    # 전부 통과했다. 어느 검사도 '이 태그가 원격에 있는가' 를 묻지 않았기
    # 때문이다. 산문을 고치는 것으로는 다음 번에 또 난다.
    #
    # 망 없이 돌리는 재현자를 막지 않으려고, 원격을 물을 수 없을 때는
    # 로컬 태그만 확인하고 그 사실을 적는다.
    import re as _r3
    import subprocess as _sp
    _tagm = _r3.search(r'\\texttt\{(acs-submission-[\w-]+)\}',
                       open(TEX).read())
    if not _tagm:
        print('[verify] 검사 A-tag: 원고가 예치 태그를 적지 않는다. 투고 '
              '직전에 태그를 만들어 푸시하고 Data Availability 에 적을 것 '
              '(main.tex 의 [확인 필요] 주석).')
    if _tagm:
        _tag = _tagm.group(1)
        try:
            _loc = _sp.run(['git', 'tag', '-l', _tag], capture_output=True,
                           text=True, timeout=20).stdout.split()
        except Exception:
            _loc = []
        if _tag not in _loc:
            fails.append(
                f'원고가 가리키는 예치 태그 {_tag} 가 로컬에 없다 '
                '(검사 A-tag): Data Availability 가 존재하지 않는 상태를 '
                '주장한다. 태그를 만들거나 원고에서 태그명을 뺄 것')
        else:
            try:
                _rem = _sp.run(['git', 'ls-remote', '--tags', 'origin'],
                               capture_output=True, text=True,
                               timeout=60)
                _ok = _rem.returncode == 0
            except Exception:
                _ok = False
            if not _ok:
                print('[verify] 검사 A-tag: 원격을 물을 수 없어 로컬 태그만 '
                      '확인했다. 투고 전에 git ls-remote --tags origin 으로 '
                      '직접 확인할 것')
            elif f'refs/tags/{_tag}' not in _rem.stdout:
                fails.append(
                    f'예치 태그 {_tag} 가 원격에 없다 (검사 A-tag): '
                    '원고 Data Availability 가 공개돼 있지 않은 상태를 '
                    '가리킨다. 태그를 푸시하거나 원고에서 뺄 것')

    # (B2d) 투고 패키지 사본이 작업본과 어긋났는지 (2026-11-07 신설).
    #
    # 왜 필요한가. submission/ 은 sync_submission.sh 가 채우는 사본이고,
    # 그 스크립트를 돌리는 것은 사람이 기억해야 하는 일이었다. 11차 푸시
    # 뒤에 확인해 보니 submission/ 은 9차 스냅샷이었고, main.tex 가
    # 11차에 철회한 제목을 들고 있었다. 스크립트 82개 중 17개가 구판이고
    # 루트에서 archive 로 옮긴 make_ga.py 의 사본이 남아 있었다.
    # 저장소를 여는 사람이 둘 중 어느 쪽을 열지는 우리가 정하지 못한다.
    #
    # 이 프로젝트에서 '두 번째 사본이 조용히 늙는' 실패가 네 번째다
    # (paper/ 대 paper_acs/, 생성표 목록, 도판 생성기 목록, README).
    # 한 번 맞추는 것으로는 닫히지 않으므로 어긋나면 멈추게 한다.
    #
    # submission/ 이 없는 상태(공개 zip 을 푼 경우)는 정상이므로 건너뛴다.
    if _o2.path.isdir('submission'):
        import filecmp as _fc
        _sd = []
        for _n in ('main.tex', 'supporting.tex', 'cover_letter.tex',
                   'refs.bib'):
            _a, _b = f'paper_acs/{_n}', f'submission/{_n}'
            if not _o2.path.exists(_b):
                _sd.append(f'  submission/{_n} 없음')
            elif not _fc.cmp(_a, _b, shallow=False):
                _sd.append(f'  submission/{_n} 가 paper_acs/{_n} 와 다르다')
        if _o2.path.isdir('submission/code'):
            _root = {_o2.path.basename(x) for x in
                     _glob.glob('*.py') + _glob.glob('*.sh')
                     if _o2.path.basename(x) != 'sync_submission.sh'}
            _root |= {'README.md', 'requirements.txt'}
            _pkg = {_o2.path.basename(x) for x in
                    _glob.glob('submission/code/*.py')
                    + _glob.glob('submission/code/*.sh')}
            _pkg |= {_n for _n in ('README.md', 'requirements.txt')
                     if _o2.path.exists(f'submission/code/{_n}')}
            for _n in sorted(_root - _pkg):
                _sd.append(f'  submission/code/{_n} 없음 (루트에 있다)')
            for _n in sorted(_pkg - _root):
                _sd.append(f'  submission/code/{_n} 가 루트에 없다 '
                           f'(옮겼거나 지운 파일의 사본)')
            _nd = [_n for _n in sorted(_root & _pkg)
                   if not _fc.cmp(_n, f'submission/code/{_n}', shallow=False)]
            if _nd:
                _sd.append(f'  submission/code 의 스크립트 {len(_nd)}개가 '
                           f'루트와 다르다: ' + ', '.join(_nd[:6])
                           + (' 등' if len(_nd) > 6 else ''))
        if _sd:
            fails.append('투고 패키지 사본이 작업본과 어긋났다 (검사 B2d): '
                         'sh sync_submission.sh 를 돌릴 것')
            fails.extend(_sd[:10])

    # (B2c) TOC 머리글이 제목과 같은 주장을 하는지 (2026-11-07 신설).
    #
    # 왜 필요한가. 9차에서 제목을 'does not predict' 에서 'no resolvable
    # agreement' 로 내렸는데 TOC 머리글만 구판으로 남았다. 독자가 가장
    # 먼저 보는 자리에 철회한 주장이 있었고, 눈으로 보고서야 찾았다.
    # 두 문장이 거명하는 작용제 집합이 같아야 한다.
    # 2026-11-07. 처음 쓸 때 _fa 를 지역 import 없이 썼고 try/except 가
    # 그 NameError 를 먹어 가드가 조용히 꺼졌다. 손으로 머리글을 틀리게
    # 바꿔 보고서야 알았다. 예외를 먹지 않고, 등재가 없으면 실패로 본다.
    import fig_annot as _fa2
    import re as _r2
    _ann = _fa2.load()
    _hl = (_ann.get('toc_acs') or {}).get('headline', '')
    if not _hl:
        fails.append('TOC 머리글이 등재부에 없다 (검사 B2c): '
                     'python3 make_toc_acs.py 를 돌릴 것')
    else:
        _ttl = _r2.search(r'\\title\{(.+?)\}', open(TEX).read(), _r2.S)
        _ttl = ' '.join(_ttl.group(1).split()) if _ttl else ''
        # 2026-11-07. 'cyclosarin' 이 'sarin' 을 품으므로 부분문자열로
        # 세면 둘을 구별하지 못한다. 단어 경계를 요구한다.
        _AG = ('cyclosarin', 'sarin', 'soman', 'tabun', 'VX')

        def _agents(_t):
            return {_a for _a in _AG
                    if _r2.search(rf'\b{_a}\b', _t, _r2.I)}

        _in_hl, _in_ttl = _agents(_hl), _agents(_ttl)
        if _in_hl != _in_ttl:
            fails.append(
                'TOC 머리글과 제목이 다른 작용제를 거명한다 (검사 B2c): '
                f'머리글 {sorted(_in_hl)} vs 제목 {sorted(_in_ttl)}')

    # (B2b) 철회한 검정이 유효한 것으로 되살아나는지 (2026-11-07 신설).
    #
    # 왜 필요한가. 11차 주입시험에서 '부호뒤집기 검정을 유효한 것으로
    # 되돌림' 이 통과했다. 부호뒤집기 검정은 설계의 성질을 보고하는 것이어서
    # 철회했는데, 철회 문장이 어떤 앵커에도 묶여 있지 않아 '우리는 이것을
    # 사전지정 검정으로 유지한다' 를 주입해도 아무도 보지 않았다. 블랙리스트는
    # 어법을 열거해야 하므로 되살아나는 꼴을 다 적을 수 없다. 대신 철회
    # 선언의 존재를 요구한다: 'sign-flip' 이 나오는 자리마다 그 근처에
    # 철회를 말하는 말이 있어야 한다.
    # 2026-11-07. 'property of the design' 을 철회어로 세면 안 된다.
    # 주입문('reports a property of the design and of the data, and we keep
    # it as the prespecified test')이 그 구절을 그대로 들고 있어서 가드가
    # 자기 창 안에서 면제돼 버렸다. 철회를 말하는 말만 센다.
    _WD = re.compile(r'withdraw|withdrew|withdrawn|'
                     r'recorded at the time|no longer|not of the data', re.I)
    _sf = []
    for m in re.finditer(r'sign-flip', body_flat):
        win = body_flat[max(0, m.start() - 200):m.end() + 420]
        if not _WD.search(win):
            _sf.append(' '.join(win[180:300].split())[:70])
    if _sf:
        fails.append('철회한 부호뒤집기 검정이 유효한 것으로 서술된다 '
                     '(검사 B2b):')
        for _t in sorted(set(_sf))[:4]:
            fails.append('  ' + _t)

    # (B3) 작용제 거리 전수 훑기 (2026-10-27 신설).
    #
    # 왜. (B) 는 등록된 어법만 본다. 2026-10-27 에 res:match 산문이
    # "cyclosarin (1.492) and tabun (1.989) ... VX at 3.076" 으로 구판을 싣고도
    # 332건 PASS 를 받았다. 같은 문장의 sarin·soman 은 전파됐는데 세 개가
    # 남았고, 그 어법을 덮는 패턴이 없었다. 같은 계열 실패가 이 세션에서
    # 네 번째다.
    #
    # 여기서는 어법을 세지 않는다. 작용제 이름 가까이에 있는 소수를 전부 모아,
    # 그 작용제가 보고된 어느 패널에서도 갖지 않는 값이면 실패로 본다.
    AGENT_ALIASES = {
        'sarin': 'sarin (GB)', 'soman': 'soman (GD)',
        'cyclosarin': 'cyclosarin (GF)', 'tabun': 'tabun (GA)',
        'sulfur mustard': 'sulfur mustard (HD)', 'lewisite': 'lewisite 1 (L)',
        'hydrogen cyanide': 'hydrogen cyanide (AC)',
        'cyanogen chloride': 'cyanogen chloride (CK)',
        'phosgene': 'phosgene (CG)', 'chloropicrin': 'chloropicrin (PS)',
        'VX': 'VX',
    }
    allowed = {}
    for _src in ('outputs/match_by_panel.json', 'outputs/match_by_panel_ascl3.json'):
        if not os.path.exists(_src):
            continue
        for _pan in json.load(open(_src)).values():
            if not isinstance(_pan, dict) or 'agents' not in _pan:
                continue
            for _a in _pan['agents']:
                allowed.setdefault(_a['agent'], set()).add(round(_a['d_norm'], 3))
    # 병합공간 말고 구조·전자 공간의 정규화거리도 본문에 나온다.
    # 2026-11-04 (작업 D). 종전에는 glob 이 접미사 붙은 사본까지 끌어왔다.
    # 그 사본 넷은 11-02 08:17 의 상관필터 교체 전 실행이고 파이프라인이
    # 다시 만들지 않는다. C_merged 정규화자가 3.5013/3.5188/3.5331 로
    # 정본 3.5590 과 다르고 1차 패널 매칭을 3~8개로 적는데 정본은 0 이다.
    # 구판 거리가 화이트리스트를 넓혀 본문의 구판 인용을 통과시킨다.
    # 파이프라인이 생성하는 한 벌만 읽는다.
    for _dt in ['outputs/distance_table_normalized.csv']:
        _d = pd.read_csv(_dt)
        for _r in _d.itertuples():
            allowed.setdefault(_r.real_agent, set()).add(round(_r.d_sim_norm, 3))
    # 거리로 읽힐 수 있는 꼴만 본다. 괄호 안의 수와 'at N' 꼴.
    NEAR = re.compile(
        r'\b(' + '|'.join(sorted(AGENT_ALIASES, key=len, reverse=True)) +
        r')\b[^.\n]{0,24}?(?:\((\d\.\d{2,3})\)|at (\d\.\d{2,3})\b)')
    # 'N from sarin' 처럼 수가 이름 앞에 오는 꼴도 같은 자격으로 본다.
    # 2026-10-27 에 커버레터의 'at 1.07 from sarin, 0.99 from soman' 이 이
    # 경로로 빠져나갔다.
    BEFORE = re.compile(
        r'(\d\.\d{2,3}) from\s+(' + '|'.join(
            sorted(AGENT_ALIASES, key=len, reverse=True)) + r')\b')
    # 같은 꼴로 omega(eV)와 Fukui 지수도 적히므로, 거리 어휘가 가까이 있는
    # 자리만 본다. 없으면 거리 진술이 아니라고 보고 넘긴다.
    DIST_CUE = re.compile(
        r'distance|normalis|separation|nearest|match|threshold|stands at|sits at',
        re.I)
    stray = []
    _hits = [(m, AGENT_ALIASES[m.group(1)], float(m.group(2) or m.group(3)),
              m.end(1)) for m in NEAR.finditer(body_flat)]
    _hits += [(m, AGENT_ALIASES[m.group(2)], float(m.group(1)), m.end(2))
              for m in BEFORE.finditer(body_flat)]
    for m, name, val, _name_end in _hits:
        ok = allowed.get(name, set())
        if not ok:
            continue
        win = body_flat[max(0, m.start() - 160):m.end() + 160]
        if not DIST_CUE.search(win):
            continue
        # '2-chloroethyl ethyl sulfide, a sulfur mustard simulant, at 2.35'
        # 처럼 작용제 이름이 모사체의 역할을 설명하는 자리는 그 작용제의
        # 거리가 아니다.
        if re.match(r'\s*simulant', body_flat[_name_end:_name_end + 12]):
            continue
        # 구판임을 밝히고 인용하는 문장은 면제한다.
        if re.search(r'earlier version|earlier analysis|superseded|put tabun at',
                     win, re.I):
            continue
        # 2026-11-03. 휘발성 축의 간격은 log10 증기압 자수이고 서술자
        # 공간의 거리가 아니다. 'decades' 가 창에 있으면 그 자리는
        # vol_gap_* 앵커가 검사한다.
        if re.search(r'decades', win):
            continue
        if not any(abs(val - v) <= 0.006 for v in ok):
            stray.append(f'{name} = {val:g}  [{" ".join(m.group(0).split())[:58]}]')
    if stray:
        fails.append(f'작용제 거리 전수 훑기 (검사 B3): 어느 패널 산출물과도 '
                     f'맞지 않는 값 {len(stray)}건')
        for _t in sorted(set(stray))[:8]:
            fails.append('  ' + _t)

    # (B4) 표 2 의 최근접 모사체 이름 (2026-10-27 신설).
    #
    # 왜. 표 행 패턴이 모사체 이름을 박아 두고 있었다. 재실행으로 HD 의
    # 최근접이 dichloromethane 에서 DPGME 로 바뀌었는데, 패턴이 'HD & DCM &'
    # 를 요구했으므로 검사가 조용히 꺼졌고 표에는 틀린 이름이 남았다.
    # 이름을 수치와 같은 자격으로 검사한다.
    SIM_ABBR = {'DMMP': 'DMMP', 'DCM': 'dichloromethane', 'DPGME': 'DPGME',
                'CEES': '2-chloroethyl ethyl sulfide', 'DFP': 'DFP',
                'DCNP': 'diethyl cyanophosphonate'}
    AGENT_ABBR = {'GB': 'sarin (GB)', 'GD': 'soman (GD)', 'GF': 'cyclosarin (GF)',
                  'GA': 'tabun (GA)', 'HD': 'sulfur mustard (HD)',
                  'L': 'lewisite 1 (L)', 'AC': 'hydrogen cyanide (AC)',
                  'CK': 'cyanogen chloride (CK)', 'CG': 'phosgene (CG)',
                  'PS': 'chloropicrin (PS)', 'VX': 'VX'}
    if os.path.exists('outputs/match_by_panel.json'):
        _mb = json.load(open('outputs/match_by_panel.json'))['panel28_incumbent3']
        _near = {a['agent']: a['nearest_simulant'] for a in _mb['agents']}
        _rowpat = re.compile(
            r'\b(' + '|'.join(AGENT_ABBR) + r') & (\S+) & \d+\.\d+ &')
        _bad = []
        for m in _rowpat.finditer(body_flat):
            _ag = AGENT_ABBR[m.group(1)]
            _got = m.group(2)
            _want = _near.get(_ag)
            if _want is None:
                continue
            if SIM_ABBR.get(_got, _got) != _want:
                _bad.append(f'{m.group(1)}: 표 {_got} vs 산출물 {_want}')
        if _bad:
            fails.append('표 2 의 최근접 모사체 이름 (검사 B4):')
            for _t in sorted(set(_bad)):
                fails.append('  ' + _t)

    # (C) 데이터가 반증하는 문구. 수치가 맞아도 주장이 틀릴 수 있다.
    # Figure 2 의 구판 주석은 숫자를 하나도 틀리지 않으면서
    # 'the least acutely toxic of the thirteen' 이라고 주장했다.
    contradictions = []
    if not f.get('_dcm_is_safest', True):
        # 부정·정정 문맥은 제외한다. '더 이상 최소독성이 아니다' 같은 서술은
        # 반증된 주장이 아니라 그 주장을 거둔 문장이다.
        NEG = re.compile(r'no longer|is not the|are not the|rather than|'
                         r'but is not|corrected|withdraw', re.I)
        for m in re.finditer(r'[^.\n]*least\s+acutely\s+toxic[^.]*\.', body):
            seg = m.group(0)
            if NEG.search(seg):
                continue
            if 'dichloromethane' in seg or 'simulant' in seg:
                contradictions.append(
                    '  반증된 주장: 본문/도판이 dichloromethane 을 '
                    f'최소독성이라 하지만 실제 순위는 {f["dcm_panel_rank"]}위이고 '
                    f'{", ".join(f["_panel_safer_than_dcm"])} 가 더 안전하다.'
                    f'  [{" ".join(seg.split())[:90]}]')

    # (B) 근사 불일치: 사실표 값에 0.5% 안쪽으로 가깝지만 같지 않은 수치
    # 근사 검사는 1000 이상의 값에만 적용한다. 반올림 분기가 실제로 문제가 된
    # 것은 천단위 구분자가 붙는 환산값이었고(41,683 대 41,685), 그 아래에서는
    # 서로 무관한 양끼리 우연히 0.5% 안쪽으로 가까워져 오경보만 낸다
    # ('270 chemical entries' 와 malathion 의 268.8배).
    targets = {k: float(v) for k, v in f.items()
               if isinstance(v, (int, float)) and abs(float(v)) >= 1000}
    near = []
    for m in re.finditer(r'(?<![\d.])(\d[\d,]*(?:\\,\d{3})*(?:\.\d+)?)', body):
        try:
            v = tonum(m.group(1))
        except ValueError:
            continue
        if v < 1000:
            continue
        raw = m.group(1)
        dec = len(raw.split('.')[1]) if '.' in raw else 0
        sig = len(re.sub(r'[^0-9]', '', raw).lstrip('0')) or 1
        for k, t in targets.items():
            if t == 0:
                continue
            # 본문이 쓴 정밀도로 반올림해 일치하면 정상 표기다.
            # 41,684.91 을 41,685 로 적는 것은 결함이 아니고,
            # 41,683 으로 적는 것은 결함이다.
            if round(t, dec) == round(v, dec):
                continue
            try:
                from decimal import Decimal
                if float(f'%.{sig}g' % t) == v:
                    continue
            except Exception:
                pass
            d = abs(v - t) / t
            if 0 < d <= 0.005:
                near.append(f'  근사 불일치 {k}: 본문 {v:,.4g} vs 산출물 {t:,.4g} '
                            f'(차이 {100*d:.3f}%)  [...{body[max(0,m.start()-40):m.end()+12]}]'
                            .replace('\n', ' '))
    near = sorted(set(near))

    print(f'[verify] 문맥 대조 {checked}건')
    # 2026-10-24. 종전에는 '문맥 미발견' 이 경고였고 판정은 PASS 였다. 그래서
    # 본문을 다시 쓰거나 표를 다시 생성해 정규식이 안 맞게 되면 검사가 조용히
    # 꺼지고 화면에는 통과가 찍혔다. 이 세션에서만 그 경로로 세 번 빠뜨렸다.
    # 꺼진 검사는 없는 검사보다 나쁘다. 실패로 올린다. 본문에 싣지 않기로
    # 한 사실은 UNQUOTED 에 이유와 함께 적어 명시적으로 면제한다.
    unexpected = [u for u in unseen
                  if u.split(':')[0] not in UNQUOTED and '도판' not in u]
    if unseen:
        print(f'[verify] 문맥 미발견 {len(unseen)}건 (검사 누락 — 정규식 보강 필요)')
        for u in unseen:
            print('  ', u)
    if unexpected:
        fails.append(f'  검사 누락 {len(unexpected)}건: '
                     + ', '.join(u.split(':')[0] for u in unexpected)
                     + ' — 정규식을 고치거나 UNQUOTED 에 면제 사유를 적을 것')
    if fails:
        print(f'[verify] 불일치 {len(fails)}건')
        print('\n'.join(fails))
    if near:
        print(f'[verify] 근사 불일치 {len(near)}건 (반올림 분기 의심)')
        print('\n'.join(near))
    if contradictions:
        print(f'[verify] 데이터가 반증하는 주장 {len(contradictions)}건')
        print('\n'.join(sorted(set(contradictions))))
    ok = not fails and not near and not contradictions
    print(f'[verify] 판정: {"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


def apply_number_fixes(dry=True):
    """--fix. 불일치한 수치를 정본 값으로 원고에 되쓴다 (2026-11-02 신설).

    왜 필요한가. 이 프로젝트의 반복 결함은 계산 오류가 아니라 전파 실패다.
    T2-1 에서 상관필터 하나를 바꾸자 검증기가 원고 수치 138건의 불일치를
    보고했다. 그만한 수를 손으로 옮기면 몇 건은 반드시 틀리고, 틀린 자리가
    어디인지도 알 수 없다. 검증기는 이미 (어느 어법, 본문 값, 정본 값,
    본문이 쓴 자릿수) 를 알고 있으므로, 같은 정보로 되쓸 수 있다.

    보수적으로 짠다.
      - 원문 파일에서 그대로 맞는 어법만 고친다. 주석 제거본이나 공백
        정규화본에서만 맞은 자리는 고치지 않고 목록으로 보고한다. 그 자리는
        줄바꿈이 끼어 있어 손으로 봐야 한다.
      - 캡처가 순수한 숫자가 아니면(단어 수사 'seventeen', 천단위 '41\\,685')
        건드리지 않는다.
      - 본문이 쓴 소수 자릿수를 유지한다. 0.0265 -> 0.0279, 0.169 -> 0.192.
      - 자릿수가 0 이면 정수로 쓴다.
    고친 뒤에는 반드시 검증기를 다시 돌려 PASS 를 확인해야 한다. 이 함수는
    판정을 내리지 않는다.
    """
    f = load_facts()
    files = [TEX] + [x for x in EXTRA_TEX if os.path.exists(x)]
    num_re = re.compile(r'^\d[\d.]*$')
    total, skipped = 0, []
    for path in files:
        s = open(path).read()
        edits = []                       # (a, b, old, new, key)
        for key, pats, tol in CHECKS:
            if key not in f:
                continue
            try:
                want = float(f[key])
            except (TypeError, ValueError):
                continue
            for pat in pats:
                for m in re.finditer(pat, s):
                    raw = m.group(1)
                    if not num_re.match(raw.strip()):
                        continue
                    try:
                        got = tonum(raw)
                    except Exception:
                        continue
                    _core = raw[:-1] if raw.endswith('.') else raw
                    _dm = re.search(r'\.(\d+)', _core)
                    _dp = len(_dm.group(1)) if _dm else 0
                    if tol <= 0.0006:
                        bad = abs(want - got) > 0.5 * (10 ** -_dp) + 1e-12
                    else:
                        bad = abs(got - want) > tol
                    if not bad:
                        continue
                    # 2026-11-03. 종전에는 캡처 전체를 새 수치로 덮었다.
                    # 패턴 다수가 ([\d.]+) 이므로 캡처가 문장 끝 마침표를
                    # 함께 먹는데(tonum 이 그 마침표를 떼어내므로 검사는
                    # 통과한다), 되쓸 때 마침표가 사라졌다. ACS 본문 네
                    # 자리에서 문장이 붙어 버렸다. 캡처의 꼬리 마침표를
                    # 보존한다.
                    _tail = '.' if raw.endswith('.') else ''
                    new = f'{want:.{_dp}f}' + _tail
                    if new == raw:
                        continue
                    a, b = m.span(1)
                    edits.append((a, b, raw, new, key))
        # 겹치는 편집은 뒤에서 앞으로 적용한다.
        edits.sort(key=lambda t: -t[0])
        done, seen = [], set()
        for a, b, old, new, key in edits:
            if any(a < pb and b > pa for pa, pb in seen):
                skipped.append(f'{path}: {key} 편집 범위 중복 — 손으로 확인')
                continue
            seen.add((a, b))
            s = s[:a] + new + s[b:]
            done.append((key, old, new))
        if done:
            total += len(done)
            print(f'  [{path}] {len(done)}건')
            for key, old, new in sorted(done):
                print(f'      {key:<22} {old} -> {new}')
            if not dry:
                open(path, 'w').write(s)
    if skipped:
        print('\n  손으로 봐야 하는 자리:')
        for x in sorted(set(skipped)):
            print(f'      {x}')
    print(f'\n[fix] {"적용 예정" if dry else "적용"} {total}건'
          + ('  (--dry-run: 파일을 쓰지 않았다)' if dry else ''))
    print('[fix] 되쓰기는 어법이 걸린 자리만 고친다. 매칭 수·차원 수처럼 '
          '문장을 고쳐야 하는 변경은 잡지 못한다. 반드시 검증기를 다시 '
          '돌릴 것.')
    return 0


if __name__ == '__main__':
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument('--fix', action='store_true',
                    help='불일치 수치를 정본 값으로 원고에 되쓴다')
    _p.add_argument('--dry-run', action='store_true',
                    help='--fix 와 함께. 무엇을 고칠지만 보고한다')
    _a, _rest = _p.parse_known_args()
    if _a.fix:
        sys.exit(apply_number_fixes(dry=_a.dry_run))
    sys.exit(main())
