#!/usr/bin/env python3
"""verify_consistency.py 가 실제로 결함을 잡는지 시험한다 (2026-10-11 신설).

통과하는 검증기는 그 자체로는 아무것도 보증하지 않는다. 종전 검증기는
핵심 수치가 본문에 '존재하는지'만 보았으므로 늘 통과했고, 같은 양이 208 과
218 로 동시에 적혀 있어도 통과했다. 그래서 과거에 실제로 발생한 결함을
본문에 한 건씩 주입하고, 검증기가 그때마다 FAIL 하는지 확인한다.

여기 실린 여섯 건은 모두 GATE 8 심사에서 지적된 실제 사례다.
"""
import atexit
import re
import shutil
import signal
import subprocess
import sys
import tempfile

TEX = 'paper/main.tex'

_NUM = r'\d+(?:\.\d+)?'


def resolve(a, b, src):
    """주입 앵커 (a, b) 를 src 의 현재 글에 맞춘다 (2026-11-02 신설).

    왜 필요한가. 상관필터를 결정적 규칙으로 교체하자 원고 수치 157건이
    갱신되었고, 이 시험의 앵커 42개가 한꺼번에 글자 그대로 사라져 전부
    SKIP 되었다. 가드 42개가 조용히 무력화된 것이고 그것은 가드가 없는
    것과 같다. 같은 종류의 무력화가 이 프로젝트에서 반복됐다.

    어법은 그대로 두고 숫자만 바뀌는 것이 정상 경로이므로, 앵커를 숫자에
    묶지 않는다. a 를 글자 그대로 찾지 못하면 a 의 숫자를 느슨하게 풀어
    같은 어법을 찾고, a 와 b 가 다른 숫자 자리를 그 찾은 글에 그대로
    옮긴다. 어법이 유일하게 걸리지 않으면 손으로 보라고 넘긴다.

    단어로 적힌 수(four, three)는 숫자가 아니므로 이 완화가 덮지 못한다.
    그런 항목은 여전히 SKIP 으로 드러나며, 그래야 한다.

    2026-11-03 추가. 표를 생성으로 돌리면 칸 사이 공백 폭이 바뀐다. 손표
    시절의 앵커는 공백 개수까지 글자 그대로 들고 있었으므로 숫자 완화만으로
    는 걸리지 않았고 앵커 5개가 SKIP 으로 남았다. 공백 연속은 길이를 묻지
    않는다.
    """
    if a in src:
        return a, b
    pa = re.split(f'({_NUM})', a)
    if len(pa) < 2:
        return None, None

    def _lit(t):
        """공백 연속은 \\s+ 로, 나머지는 글자 그대로."""
        return ''.join(r'\s+' if p.isspace() else re.escape(p)
                       for p in re.split(r'(\s+)', t) if p)

    # 엄격(공백까지 글자 그대로) 먼저, 유일하지 않으면 공백 완화. 순서를
    # 바꾸면 안 된다. 완화를 먼저 쓰면 종전에 유일하게 걸리던 앵커 두 개가
    # 복수 일치로 바뀌어 거꾸로 SKIP 된다 (2026-11-03 실측).
    ms = []
    for _lax in (False, True):
        pat = ''.join((_lit(t) if _lax else re.escape(t))
                      if i % 2 == 0 else _NUM
                      for i, t in enumerate(pa))
        ms = list(re.finditer(pat, src))
        if len(ms) == 1:
            break
    if len(ms) != 1:
        return None, None
    found = ms[0].group(0)
    na, nb = re.findall(_NUM, a), re.findall(_NUM, b)
    if len(na) != len(nb) or na == nb:
        return None, None            # 숫자 아닌 곳이 달라졌다. 손으로.
    pf = re.split(f'({_NUM})', found)
    for i, (x, y) in enumerate(zip(na, nb)):
        if x != y:
            pf[2 * i + 1] = y
    return found, ''.join(pf)

# 2026-11-02. 이 시험이 타임아웃으로 중단되었을 때 주입값(0.0221 자리에
# 0.0321)이 paper/main.tex 에 그대로 남았다. finally 는 SIGTERM 을 막지
# 못한다. 원본을 프로세스 수명에 걸어 두고, 어떤 경로로 죽든 되돌린다.
_PRISTINE = {}


def _guard(path):
    """path 의 현재 내용을 원본으로 등록한다. 중복 등록은 무시한다."""
    if path not in _PRISTINE:
        _PRISTINE[path] = open(path).read()


def _restore_all():
    for path, text in _PRISTINE.items():
        try:
            if open(path).read() != text:
                open(path, 'w').write(text)
                print(f'  [복원] {path} 를 원본으로 되돌렸다', file=sys.stderr)
        except OSError as e:
            print(f'  [복원 실패] {path}: {e}', file=sys.stderr)


atexit.register(_restore_all)


def _on_signal(signum, frame):
    _restore_all()
    sys.exit(128 + signum)


for _sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
    try:
        signal.signal(_sig, _on_signal)
    except (ValueError, AttributeError):
        pass

# (설명, 원문, 주입할 오류)
CASES = [
    ('한 양을 두 값으로 (208/218 계열)',
     'association across the 207',
     'association across the 218'),
    ('초록이 Results 와 다른 Spearman 을 싣는다',
     '(Spearman $r=+0.062$',
     'independent (Spearman $r=-0.022$'),
    # --- 롤백 C 헤드라인 (2026-10-13 추가). 이 수치들이 주장을 떠받친다.
    ('BH 보정값 오기 (0.159 대 0.059)',
     'adjusted value of 0.181',
     'adjusted value of 0.059'),
    ('무작위 3종 부분공간 비율 오기 (36.3 대 6.3)',
     'and 33.7\\% of the subspaces reach',
     'and 3.7\\% of the subspaces reach'),
    ('주성분 3개 절단 p 오기 (0.0044 대 0.0440)',
     'from $p=0.026$ to\n$p=0.0029$ on the twenty-eight',
     'from $p=0.026$ to\n$p=0.0290$ on the twenty-eight'),
    ('Mantel p 를 쌍독립 가정 p 로 바꿔 씀 (0.011 대 0.048)',
     '$\\rho = +0.514$, Mantel $p = 0.015$',
     '$\\rho = +0.514$, Mantel $p = 0.048$'),
    ('omega 상관계수 오기 (-0.899 대 -0.959)',
     'r(E_{\\mathrm{LUMO}}, \\omega) =\n-0.899',
     'r(E_{\\mathrm{LUMO}}, \\omega) =\n-0.959'),
    ('GA 코사인 오기 (0.764 대 0.564)',
     'tabun, at $+0.764$ against sarin',
     'tabun, at $+0.564$ against sarin'),
    ('백분위 오기 — 구조 관측값을 더 특별해 보이게 (25.8 대 4.2)',
     'value of 0.028 sits at the\n22.4th percentile',
     'value of 0.026 sits at the 4.2nd percentile'),
    ('주성분3 백분위 오기 (2.1 대 0.4)',
     'value of 0.028 sits at the\n22.4th',
     'value of 0.026 sits at the 0.4th'),
    ('주 패널 유효매칭 수 오기 (zero 대 three)',
     'twenty-eight-compound panel and zero on the primary one',
     'twenty-eight-compound panel and three on the primary one'),
    ('28종 유효매칭 수 오기 (three 대 five)',
     'against\nthree on the twenty-eight-compound panel',
     'against\nseven on the twenty-eight-compound panel'),
    # 2026-11-03. 종전 앵커 '& 0.192 & 0.039 & 0.025' 는 숫자만으로 되어
    # 있어서, tab:family 를 생성으로 돌린 뒤 생성 표 안에서 다섯 자리에
    # 걸렸다. 복수 일치는 SKIP 이므로 가드가 꺼진다. 행 이름을 닻으로
    # 넣어 유일하게 만든다.
    ('WY 조정 p 오기 (0.079 대 0.045) — 유의/비유의가 뒤집힌다',
     'Westfall-Young min-$P$ over 18           & 0.457 & 0.112 & 0.079',
     'Westfall-Young min-$P$ over 18           & 0.457 & 0.112 & 0.045'),
    ('부트스트랩 0건 비율 오기 (58.8 대 98.8)',
     'exactly zero\nin 58.6\\% of replicates',
     'exactly zero\nin 98.6\\% of replicates'),
    ('부트스트랩 차 양수 복제수 오기 (1000 대 600)',
     'positive in 1000 of 1000 replicates',
     'positive in 600 of 1000 replicates'),
    ('GB 회복률 오기 (37.8 대 7.8) — 가장 약한 지점을 숨김',
     'as its nearest compound in 38.0\\% of replicates',
     'as its nearest compound in 18.0\\% of replicates'),
    ('35종 귀무기대 오기 (9.63 대 6.63)',
     'expectation of 9.63 set by the 21:3',
     'expectation of 6.63 set by the 21:3'),
    ('최소검출효과 오기 (0.22 대 0.42) — 검정력 주장이 뒤집힌다',
     'moved 0.22 typical nearest-neighbour',
     'moved 0.42 typical nearest-neighbour'),
    ('실제 결손 중위 오기 (0.463 대 0.163) — 검출 가능 폭보다 작아진다',
     'carry has a median of\n0.445',
     'carry has a median of\n0.145'),
    ('사린 결손 오기 (0.117 대 0.317)',
     'only at sarin, where it is 0.110',
     'only at sarin, where it is 0.317'),
    ('CEES-HD 거리 오기 (0.57 대 1.57)',
     'sulfide at 0.54 from sulfur mustard',
     'sulfide at 1.57 from sulfur mustard'),
    ('DFP 변위 거리 오기 (1.07 대 2.07) — 변위 논증이 무너진다',
     'fluorophosphate is present, at 1.07, 1.00 and 1.12',
     'fluorophosphate is present, at 2.07, 1.00 and 1.12'),
    ('사전지정 q 오기 (0.046 대 0.146) — 유의/비유의가 뒤집힌다',
     'would survive at $q=0.042$',
     'would survive at $q=0.142$'),
    ('백색화 통제 수치 오기 (16.9 대 36.9)',
     'the rate falls from 33.7\\% to 18.4\\%',
     'the rate falls from 33.7\\% to 38.4\\%'),
    ('meff 오기 (5.00 대 15.00)',
     'estimates\n5.00 effective tests',
     'estimates\n15.00 effective tests'),
    # 2026-11-03. 앵커에 행 이름이 없어 생성 표의 전자공간 행 네 개에
    # 모두 걸렸다. 복수 일치는 SKIP 이고 SKIP 은 꺼진 가드다.
    ('표 1 의 28종 전자공간 행이 구판 (0.275 대 0.297)',
     '\\quad (omits the seven) & electronic  &  3 & 0.275 & 0.175',
     '\\quad (omits the seven) & electronic  &  3 & 0.297 & 0.175'),
    # --- 네 번째 구성 (2026-10-20). 이 논문에서 다중검정 보정을 통과하는
    # 유일한 구성이므로 수치 하나가 틀리면 주장이 통째로 바뀐다.
    # 2026-11-03. 숫자 완화가 '28' 과 '35' 를 함께 풀어 두 블록의 구조공간
    # 행에 모두 걸렸다. 다섯째 칸까지 넣어 모양으로 구분한다 (35종 행은
    # 그 자리가 \textbf{$<$0.001} 이다).
    ('표 1 셋째 블록 구조공간 행이 틀림 (0.028 대 0.128)',
     '28 cpd., incumbent 3   & structural  & 35 & \\textbf{0.02794} & 0.1103'
     ' & 0.1132',
     '28 cpd., incumbent 3   & structural  & 35 & \\textbf{0.1279} & 0.1103'
     ' & 0.1132'),
    ('표 1 네 번째 블록 병합공간 행이 틀림 (0.005 대 0.055)',
     '+ pesticides 2   & merged      & 42 & \\textbf{0.006}',
     '+ pesticides 2   & merged      & 42 & \\textbf{0.056}'),
    ('네 번째 구성 BH q 를 비유의로 오기 (0.0006 대 0.0600) — 결론이 뒤집힌다',
     'Benjamini-Hochberg value of 0.0011 over the same',
     'Benjamini-Hochberg value of 0.0600 over the same'),
    ('표 2 네 번째 열의 사전지정 값 오기 (0.015 대 0.150) — 유의성이 뒤집힌다',
     'Statistic specified in advance (3 tests) & 0.386 & 0.096 & \\textbf{0.040} & \\textbf{0.017}',
     'Statistic specified in advance (3 tests) & 0.386 & 0.096 & \\textbf{0.040} & \\textbf{0.170}'),
    ('표 2 네 번째 열의 meff 오기 (5.00 대 15.00)',
     '$m_{\\mathrm{eff}}$ of 18 (Li-Ji)     & 6.00  & 5.00  & 5.00  & 5.00',
     '$m_{\\mathrm{eff}}$ of 18 (Li-Ji)     & 6.00  & 5.00  & 5.00  & 15.00'),
    ('서술자 제거 수 오기 (13 대 23) — 51 과 38 이 안 맞게 된다',
     'remove 13 of the 51 on the primary panel',
     'remove 23 of the 51 on the primary panel'),
    ('초록의 네 번째 구성 보정값 오기 (0.0006 대 0.0060)',
     'produces an adjusted 0.0011.',
     'produces an adjusted 0.0110.'),
    ('결론의 네 번째 구성 보정값 오기 (0.0006 대 0.0006 유지 실패)',
     'into one with an adjusted value of\n0.0011.',
     'into one with an adjusted value of\n0.0260.'),
    ('네 번째 구성 최소 nominal p 오기 (0.00006 대 0.00600)',
     'smallest nominal value of 0.00006 and a',
     'smallest nominal value of 0.00600 and a'),
    # 2026-11-03 폐기. C(24,10) 전수열거 횟수는 농약 구성 절을 줄이면서
    # 본문에서 빠졌고, verify_consistency 의 UNQUOTED 에 사유와 함께
    # 면제로 올렸다. 본문에 없는 문장에 대한 주입 시험은 영구 SKIP 이
    # 되어 '꺼진 가드'로 남으므로 여기서도 함께 내린다. 그 수치를 본문에
    # 다시 싣는다면 UNQUOTED 항목과 이 항목을 같이 되살릴 것.
    # ('네 번째 구성 순열 수 오기 (1,961,256 대 1,916,256)',
    #  '\\binom{24}{10}=1\\,961\\,256',
    #  '\\binom{24}{10}=1\\,916\\,256'),
    # --- 인용 문헌에서 옮겨 적은 수치 (2026-10-20).
    ('James 2018 부적합 판정 수 오기 (17 대 12)',
     'judged 17 unsuitable',
     'judged 12 unsuitable'),
    ('James 2021 배제 수 오기 (78 대 87)',
     'discounted 78 on the risk',
     'discounted 87 on the risk'),
    ('James 2021 적합 수 오기 (39 대 93)',
     'repository of 39 suitable',
     'repository of 93 suitable'),
    ('블록균등 최소 p 오기 (0.076 대 0.176) — 가중 의존성 진술이 무력화된다',
     'nominal value from 0.151 to 0.061',
     'nominal value from 0.151 to 0.176'),
    # --- 표 3 (2026-10-17 추가). 이 표는 그때까지 검증기에 등록되지 않아
    # 11행 중 9행이 구판 값인 채로 144건 검사를 통과했다. 동료 검토에서
    # 지적됐다. 반올림으로 설명되지 않던 두 행과 산문·표 불일치를 시험한다.
    ('표 3 PS 행이 구판 (2.530 대 2.70) — 반올림으로 설명 안 됨',
     'PS & DMMP & 2.433 & no',
     'PS & DMMP & 2.70 & no'),
    ('표 3 CG 행이 구판 (1.600 대 1.64)',
     'CG & DCM & 1.600 & no',
     'CG & DCM & 1.64 & no'),
    ('표 3 GA 행이 임계값 2 를 넘지 않게 보이도록 오기',
     'GA & DMMP & 1.995 & yes',
     'GA & DMMP & 1.903 & yes'),
    ('산문이 표 3 과 다른 소만 거리를 싣는다 (1.395 대 1.495)',
     'soman (1.445),',
     'soman (1.495),'),
    ('전자공간 유효매칭 개수 주장이 산문에만 있고 거리는 구판 (1.407 대 1.335)',
     'dichloromethane at 1.407,',
     'dichloromethane at 1.335,'),
    ('포스겐 경쟁 거리를 다른 공간의 값으로 섞어 씀 (1.12 대 1.14)',
     'and in the merged space (1.120 against',
     'and in the merged space (1.240 against'),
    # 2026-10-23. 역방향 최근접 작용제 부트스트랩.
    ('DPGME 의 G계열 비율이 구판이다 (46.7 대 61.2)',
     'the\nG series in 56.9\\%', 'the\nG series in 26.9\\%'),
    ('DMMP 의 사린 최근접 비율 오기 (93.3 대 88.1)',
     'DMMP is nearest sarin in\n92.8\\% of replicates',
     'DMMP is nearest sarin in\n82.8\\% of replicates'),
    ('DCM 의 G계열 비율 오기 (2.3 대 12.3)',
     'Dichloromethane is nearest a G-series agent in 4.2\\%',
     'Dichloromethane is nearest a G-series agent in 14.2\\%'),
    # 2026-10-20. 독성 도판이 본문 Figure 3 으로 올라와 주입 대상도 본문이다.
    ('Fig3 라벨 작용제 수 오기 (ten 대 nine)',
     "name the ten agents of this paper's panel",
     "name the nine agents of this paper's panel"),
    ('Fig3 질소머스터드 수 오기 (three 대 two)',
     'three nitrogen mustards', 'two nitrogen mustards'),
    ('Fig3 무라벨 작용제 수 오기 (seven 대 six)',
     'the remaining seven held-out agents are plotted',
     'the remaining six held-out agents are plotted'),
    ('도판 캡션의 각주 유래 개수 (5 대 9)',
     'Five of the nineteen values come from the',
     'Nine of the nineteen values come from the'),
    ('적용영역 내 최대오차를 전체 최대로 오기 (2.68 대 3.42)',
     'wrong by 2.68 log units. The dotted line',
     'wrong by 3.42 log units. The dotted line'),
]

# 커버레터·하이라이트도 검사 대상이다 (2026-10-13 추가).
# 2026-11-03 신설. SI 의 네 표(ablation / family / blockbalance / sparse)와
# plocal 을 make_si_tables.py 생성으로 돌렸다. 그 결과 이 시험의 앵커 8개가
# 원고 .tex 에서 사라져 SKIP 되었다 — 가드 8개가 조용히 꺼진 것이다.
#
# 생성된 표에 주입할 때 주의할 점이 하나 있다. 사본 한쪽만 고치면 검증기의
# '생성 표 사본 시차' 가드가 먼저 걸리고, 시험은 FAIL 을 받아 통과한 것처럼
# 보인다. 값 검사가 아니라 다른 가드가 잡은 것이므로 거짓 보증이다. 그래서
# outputs/ 원본과 두 사본 세 곳을 함께 고친다. 그러면 시차 가드는 조용하고
# 값 검사만 발동할 수 있다.
GEN_TABLES = ('tab_ablation.tex', 'tab_family.tex', 'tab_blockbalance.tex',
              'tab_sparse.tex', 'tab_plocal.tex')


def gen_paths(name):
    return (f'outputs/{name}', f'paper_acs/{name}', f'paper/{name}')


def find_gen(a, b):
    """생성 표에서 앵커를 찾는다. (이름, 해결된 a, 해결된 b) 또는 None."""
    for name in GEN_TABLES:
        try:
            src = open(f'outputs/{name}').read()
        except OSError:
            continue
        ra, rb = resolve(a, b, src)
        if ra:
            return name, ra, rb
    return None


SIDECAR = 'paper/cover_letter.tex'
ACS = 'paper_acs/main.tex'

# 2026-11-02. ACS 판본 전용 주입 사례. 둘 다 실제로 발생한 결함이다.
# 세 원논문(hazref, multiplexref, group3ref)의 Methods 를 직접 읽기 전에,
# 어레이 농도를 ppm 으로 적었고 (실제 µM), 같은 화합물의 서로 다른 농도를
# '반복측정' 으로 적었다 (원자료는 조합당 한 행). 수행하지 않은 측정을
# 서술한 것이므로 [C-3] 위반이고, 되돌아오면 잡혀야 한다.
ACS_CASES = [
    # 2026-11-03. 종전 앵커는 표지서신의 농도 추세 문장이었고 그 주장을
    # 뺐으므로 자리가 없어졌다. 단위가 틀리면 가장 치명적인 자리, 즉
    # 실험절의 노출 농도 목록으로 옮긴다.
    ('어레이 농도를 ppm 으로 되돌림 (정본 µM)',
     'final concentrations of 0, 10, 50, 100 and\n500~$\\mu$M',
     'final concentrations of 0, 10, 50, 100 and\n500~ppm'),
    ('농도축 추세를 ppm 으로 되돌림 (정본 µM)',
     'from $+0.363$ at 10~$\\mu$M to\n$+0.771$ at 500~$\\mu$M',
     'from $+0.363$ at 10~ppm to\n$+0.771$ at 500~ppm'),
    ('수행하지 않은 반복측정 주장 (본문)',
     "One compound's patterns at two\nconcentrations correlate at",
     'Replicates of one compound correlate at'),
    # 2026-11-02 재집필로 방법 문장이 2.1 에서 Experimental 로 옮겨갔다.
    # 앵커를 현행 본문에 다시 건다.
    ('수행하지 않은 반복측정 주장 (전제 문단)',
     'itself at a different dose and not agreement between repeated\nmeasurements',
     'itself at a different dose, that is between repeated\nmeasurements'),
    ('수행하지 않은 반복측정 주장 (도판 캡션)',
     'within one compound across concentrations\n  and between different compounds',
     'within replicates of the same compound and\n  between different compounds'),
    # 2026-11-02. ACS 2.1 전제가 세 원논문의 분류 성능 위에 서 있다. 그
    # 수치가 조용히 바뀌면 전제가 근거 없이 바뀌는 것이므로 잡아야 한다.
    ('문헌값 d1 정확도 오기 (정본 100%)',
     'dimethyl methylphosphonate at 100\\%\naccuracy',
     'dimethyl methylphosphonate at 98\\%\naccuracy'),
    ('문헌값 d1 가시광 단독 오기 (정본 87.5%)',
     'against 87.5\\% under visible illumination alone',
     'against 82.5\\% under visible illumination alone'),
    ('문헌값 d2 교차검증 정확도 오기 (정본 0.964)',
     'all seven agents at an accuracy of\n0.964',
     'all seven agents at an accuracy of\n0.984'),
    ('문헌값 d3 LOO 정확도 오기 (정본 1.000)',
     'four Novichok agents at an accuracy of 1.000',
     'four Novichok agents at an accuracy of 0.980'),
    ('문헌값 d3 검출한계 오기 (정본 10 µM)',
     'with a detection limit near 10~$\\mu$M',
     'with a detection limit near 1~$\\mu$M'),
    # 2026-11-02. Conclusions 가 2.4·2.5 의 수치를 다시 센다. 두 자리가
    # 어긋나면(한쪽만 갱신되면) 잡혀야 한다. 과거 실패 계열이 바로 이것이다.
    ('결론부 BH 보정값이 2.4 와 어긋남 (정본 0.6299)',
     'an incumbent simulant, at a\nBenjamini--Hochberg adjusted value of 0.6299',
     'an incumbent simulant, at a\nBenjamini--Hochberg adjusted value of 0.1695'),
    ('결론부 AEGL 훈련 화합물 수 오기 (정본 230)',
     'random forest trained on 230 AEGL chemicals is wrong in one direction',
     'random forest trained on 320 AEGL chemicals is wrong in one direction'),
]

# 2026-11-02. ACS SI (paper_acs/supporting.tex) 전용.
# 이 파일은 paper/supplementary.tex 의 평행 사본이다. 한쪽만 고치면 조용히
# 갈라지므로 ACS 쪽 단독 오기도 잡히는지 확인한다.
ACS_CL = 'paper_acs/cover_letter.tex'

# 2026-11-02. ACS 커버레터 전용. 커버레터가 본문 수치를 다시 세므로 한쪽만
# 바뀌면 투고 패키지 안에서 모순이 난다. 처음 앵커를 달기 전에는 네 자리
# 가운데 셋이 통과했다.
ACS_CL_CASES = [
    # 2026-11-03 (T1-10). 표지서신에서 분류 정확도 주장을 뺐으므로 그
    # 케이스는 주입할 자리가 없다. 가드 수가 조용히 줄지 않도록, 표지서신이
    # 지금 싣는 주장으로 교체한다.
    ('커버레터 MDE 오기 (정본 0.67)',
     'resolves a rank correlation of 0.67 at 80',
     'resolves a rank correlation of 0.37 at 80'),
    ('커버레터 DMMP-GA 코사인 오기 (정본 +0.76)',
     'correlates at $+0.76$ with tabun', 'correlates at $+0.86$ with tabun'),
    ('커버레터 DMMP-GF 코사인 부호자리 오기 (정본 -0.72)',
     '$-0.72$ with cyclosarin', '$-0.62$ with cyclosarin'),
    # 2026-11-03 (T1-9). 농도 추세 주장도 표지서신에서 뺐다. 같은 이유로
    # 최대 관측 rho 주장으로 교체한다.
    ('커버레터 최대 관측 rho 오기 (정본 0.57)',
     'where the largest\nobserved is 0.57', 'where the largest\nobserved is 0.87'),
]

ACS_SUPP = 'paper_acs/supporting.tex'
ACS_SUPP_CASES = [
    ('ACS SI 단독으로 모델 성능 오기 (정본 R^2 0.458)',
     r'R^{2}=0.458 \pm', r'R^{2}=0.389 \pm'),
    ('ACS SI 단독으로 적용영역 거리 오기 (정본 3.28)',
     "DMMP's domain distance, 3.28,", "DMMP's domain distance, 3.09,"),
]
# 2026-10-24. 공저자 심사가 SI 의 검증 문장 두 개를 반박했고 둘 다 사실이었다.
# 그 문장들이 다시 틀어지면 잡히도록 건다.
SUPP_CASES_EXTRA = [
    # 2026-10-24 신설. 중복 채널 축약은 저자 확인(형광 미발생, ROI 가 배경을
    # 읽음)을 거쳐 들어간 새 수치다. 본문이 유효 채널 수를 구판 29 로 되돌리거나
    # 화합물 단위 검정력을 쌍 단위 값으로 되돌리면 잡혀야 한다.
    ('유효 채널 수를 선언 채널 수로 되돌림 (23 대 29)',
     'leaves 23 effective channels',
     'leaves 29 effective channels'),
    ('세 번째 세트 유효 채널 오기 (25 대 29)',
     "set's ultraviolet data and 25 in the third",
     "set's ultraviolet data and 29 in the third"),
    ('검정력을 쌍 단위 값으로 되돌림 (0.924 대 0.669)',
     'compounds as the unit, which is the unit the Mantel test uses, returns $0.924$',
     'compounds as the unit, which is the unit the Mantel test uses, returns $0.669$'),
    ('SI 의 양방향 대조 건수 오기',
     'agree on all 250 entries', 'agree on all 240 entries'),
    ('SI 의 1000 이상 값 개수 오기',
     'Thirty-four values reach 1000', 'Thirty-three values reach 1000'),
]

SIDECAR_CASES = [
    ('커버레터가 본문과 다른 Spearman 을 싣는다',
     '$r = +0.062$ (95\\% CI $-0.077$ to\n$+0.199$)',
     '$r = -0.022$ (95\\% CI $-0.156$ to\n$+0.111$)'),
]

# 도판 주석 결함은 outputs/figure_annotations.json 에 주입한다.
ANNOT = 'outputs/figure_annotations.json'
ANNOT_CASES = [
    ('도판 주석이 데이터가 반증하는 주장을 한다',
     'Fig_toxicity', 'panel_a',
     'dichloromethane, an incumbent simulant with no designating document, is the least acutely toxic '
     'of the thirteen.'),
    ('도판 주석의 개수가 구판이다',
     'Fig_toxicity', 'panel_b_header',
     'no regulator-grade acute inhalation value exists for any of these 15'),
    # 2026-10-22. 그래픽 초록이 구판 수치를 들고 투고 직전까지 남아 있었다.
    ('그래픽 초록의 보정값이 구판이다',
     'graphical_abstract', 'panel_a',
     'corrected over the 18 tests: $q$ = 0.212'),
    ('그래픽 초록의 Spearman 이 구판이다',
     'graphical_abstract', 'panel_c',
     'Spearman $r$ = -0.022   95% CI [-0.156, +0.111]   $n$ = 214\n'
     "faithful + safe: 0 at the incumbent simulant's own bar,\n"
     '8 one order of magnitude below it'),
]


# 아래 8건은 2026-10-14 분량 재배치로 supplementary.tex 로 옮겨간 문장이다.
# main.tex 에만 주입하던 종전 설계로는 전부 SKIP 되어, 검증기가 SI 를 덮는지
# 확인할 수 없었다. 주입 대상을 SI 로 돌린다.
SUPP = 'paper/supplementary.tex'
SUPP_CASES = [
    ('천단위 환산값 반올림 분기 (41,685 대 41,683)',
     r'dichloromethane (41\,685)', r'dichloromethane (41\,683)'),
    ('적용영역 거리 오기 (3.28 대 3.09)',
     "DMMP's domain distance, 3.28,", "DMMP's domain distance, 3.09,"),
    ('구판 모델 성능 잔존 (R^2 0.458 대 0.389)',
     r'R^{2}=0.458 \pm', r'R^{2}=0.389 \pm'),
    ('후보 수 불일치 (8 대 7)',
     'Eight compounds then combine', 'Seven compounds then combine'),
    ('held-out 작용제 수 불일치 (20 대 17)',
     'twenty agents held out', 'seventeen agents held out'),
    ('충실도 비율 오기 (3.1 대 3.6)',
     'trichloroethylene, is 3.1 times closer',
     'trichloroethylene, is 3.6 times closer'),
    # --- 블록 가중 민감도 (2026-10-17 추가).
    ('SI 블록균등 표의 28종 q 가 0.05 밑으로 오기',
     'Twenty-eight, incumbent three & 4/11 & 3/11 & 0.025 & 0.009 & 0.185 & 0.146',
     'Twenty-eight, incumbent three & 4/11 & 3/11 & 0.025 & 0.009 & 0.185 & 0.046'),
    ('블록균등에서 유효매칭이 변한 것처럼 오기 (0/11 대 2/11)',
     'Thirty-five, incumbent three & 0/11',
     'Thirty-five, incumbent three & 2/11'),
]
SUPP_CASES_DFT = [
    # 2026-10-25 신설. DFT 대조 수치는 산출물(outputs/dft_vs_gfn2.csv)에 묶여
    # 있다. 원고가 그 값에서 벗어나면 잡혀야 한다. 특히 omega 축의 '결론은
    # 유지, 순서는 철회' 서술은 세 숫자(1.90, 2.08, 7)에 전부 걸려 있다.
    ('DFT omega 상한 오기 (7.92 대 9.92)',
     'range becomes 1.43 to 7.92~eV',
     'range becomes 1.43 to 9.92~eV'),
    ('GFN2 omega 하한 오기 (1.71 대 2.71)',
     'they run from 1.71 to',
     'they run from 2.71 to'),
    ('DMMP 가 사린보다 높아지도록 뒤집음 (1.90 대 2.90)',
     'DMMP at 1.90 against sarin at 2.08',
     'DMMP at 2.90 against sarin at 2.08'),
    ('사린 omega 오기 (2.08 대 2.80)',
     'DMMP at 1.90 against sarin at 2.08',
     'DMMP at 1.90 against sarin at 2.80'),
    ('사린 상위 간섭물 수 오기 (seven 대 nine)',
     'and seven of the fourteen interferents still score above sarin',
     'and nine of the fourteen interferents still score above sarin'),
    ('최저 작용제 VX 값 오기 (1.54 대 1.24)',
     '(1.54), sulfur mustard (1.73)',
     '(1.24), sulfur mustard (1.73)'),
    ('omega 순위상관 오기 (0.628 대 0.828)',
     'values is 0.628, against',
     'values is 0.828, against'),
    ('LUMO 순위상관 오기 (0.898 대 0.698)',
     '0.898 for the underlying LUMO energies',
     '0.698 for the underlying LUMO energies'),
]

SUPP_CASES_STALE = [
    # 2026-10-25 신설. 재분석 뒤 전파되지 않아 16건이 PASS 를 받은 채 남아
    # 있던 수치들. (B) 는 어법 단위라 다른 문장에 남은 구판을 못 봤다.
    # 여기서는 각 자리를 실제로 구판으로 되돌려 검사 B2 또는 (B) 가 잡는지
    # 확인한다. 구판 값은 커밋 aa1d01a 의 원고에서 가져왔다.
    ('28종 구조공간 p 를 구판으로 (0.031 대 0.028)',
     '0.028 sits at the\n22.4th percentile',
     '0.028 sits at the 24.3rd percentile'),
    ('부분공간 백분위를 구판으로 (25.4 대 22.7)',
     '0.028 sits at the\n22.4th percentile of the resulting distribution',
     '0.026 sits at the 21.3rd percentile of that\ndistribution'),
    ('도판2 캡션 백분위 오기 (26th 대 21st)',
     '0.026 is the 24th percentile of those draws',
     '0.026 is the 21th percentile of those draws'),
    ('28종 매치 수를 구판으로 (four 대 three)',
     'Omitting those five restores three matches',
     'Omitting those five restores four matches'),
    ('구조공간 유의 통계량 수를 구판으로 (two 대 one)',
     'compounds under two of six statistics in the structural space',
     'compounds under one of six statistics in the structural space'),
    ('Westfall-Young 을 구판으로 (0.079 대 0.075)',
     'rather than approximating it, gives 0.078.',
     'rather than approximating it, gives 0.075.'),
    ('사전지정 q 를 구판으로 (0.046 대 0.042)',
     'merged results would survive at $q=0.042$',
     'merged results would survive at $q=0.046$'),
    ('도판1 캡션 BH 를 구판으로 (0.185 대 0.159)',
     'value over its eighteen tests is 0.165.}',
     'value over its eighteen tests is 0.159.}'),
    ('28종 구조좌표 수를 구판으로 (36 대 38)',
     'the observed 36-dimensional value of',
     'the observed 38-dimensional value of'),
    ('두 패널 보존 좌표를 같게 (38/36 대 38/38)',
     'leaving 38 and 36 structural coordinates',
     'leaving 38 and 38 structural coordinates'),
    ('SI tabun 거리를 구판으로 (2.0038 대 1.9948)',
     'its normalised distance is\n2.0038 against a threshold of 2',
     'its normalised distance is\n1.9948 against a threshold of 2'),
    ('SI 서술 BH 를 구판으로 (0.185 대 0.159)',
     'eighteen tests whose smallest adjusted value is 0.185',
     'eighteen tests whose smallest adjusted value is 0.159'),
    ('도판2 캡션 통계량 서술을 구판으로 (평균차 대 여섯 중 최소)',
     '(c) The mean-difference\nstatistic for 2000 random three-descriptor subspaces',
     '(c) Smallest $p$ over\nthe six statistics for 2000 random three-descriptor subspaces'),
]

SUPP_CASES_CONF = [
    # 2026-10-26 신설. 형태이성질체 민감도 수치는 산출물에 묶여 있다.
    # 원고가 그 범위를 좁히거나 매칭 수를 바꾸면 잡혀야 한다.
    ('앙상블 tabun 하한 오기 (1.896 대 1.796)',
     'runs from 1.9943 to 2.0106 across the ensemble seeds',
     'runs from 1.7943 to 2.0106 across the ensemble seeds'),
    ('앙상블 tabun 상한 오기 (1.905 대 1.805)',
     'runs from 1.9943 to 2.0106 across the ensemble seeds',
     'runs from 1.9943 to 1.8106 across the ensemble seeds'),
    ('단일선택 tabun 범위 오기 (1.876 대 1.976)',
     'and from 2.0172 to\n2.0273 across the single-conformer',
     'and from 1.7172 to\n2.0273 across the single-conformer'),
    ('단일선택 순열 p 상한 오기 (0.028 대 0.038)',
     'permutation value runs from 0.0250 to 0.0279 under both',
     'permutation value runs from 0.0250 to 0.0379 under both'),
    ('시드 수 오기 (six 대 four)',
     'under six seeds and rebuilding',
     'under four seeds and rebuilding'),
    ('구조좌표 수를 메타열 포함 시절로 (36 대 37)',
     'value of 0.028 sits at the\n22.4th',
     'value of 0.029 sits at the 24.3rd'),
    ('구조공간 명목 p 를 앙상블 이전으로 (0.026 대 0.031)',
     '0.028 sits at the\n22.4th percentile',
     '0.031 sits at the 24.3rd percentile'),
    ('유효검정수를 구판으로 (5.00 대 7.00)',
     '5.00 effective tests of eighteen',
     '7.00 effective tests of eighteen'),
]

SUPP_CASES_DIPOLE = [
    # 2026-10-27 신설. 쌍극자까지 B3LYP 로 바꾼 3좌표 배치와 병합공간 대조.
    # 종전 원고는 쌍극자를 GFN2 로 남긴 채 '어떤 통상 수준에서도 유의하지
    # 않다'고 적었고, 그 문장은 세 좌표를 다 바꾸면 거짓이 된다. 이 숫자
    # 다섯 개가 그 교정을 떠받친다.
    ('3좌표 B3LYP p 오기 (0.0368 대 0.1368)',
     'Replacing all three gives\nfive valid matches and 0.0368',
     'Replacing all three gives\nfive valid matches and 0.1368'),
    ('3좌표 B3LYP 매칭 수 오기 (five 대 four)',
     'Replacing all three gives\nfive valid matches',
     'Replacing all three gives\nfour valid matches'),
    ('비유의 주장 복원 (정본 p=0.0368)',
     'return 0.0368. We do not promote',
     'return 0.0368. The permutation test is still not\nsignificant at any '
     'conventional level. We do not promote'),
    ('병합공간 B3LYP p 오기 (0.0221 대 0.1221)',
     'a one-sided $p$ of 0.0221, against four and',
     'a one-sided $p$ of 0.1221, against four and'),
    ('병합공간 B3LYP 매칭 수 오기 (three 대 five)',
     'analysis and returns three valid matches and\na one-sided',
     'analysis and returns five valid matches and\na one-sided'),
    ('병합공간 GFN2 대조값 오기 (0.0250 대 0.0350)',
     'against four and 0.0250 for the GFN2 control',
     'against four and 0.0350 for the GFN2 control'),
    ('5건 전부 기록 지지 짝이라고 서술 (정본 4/5)',
     'Four of those five are sarin',
     'The five are pairings the designation record supports. Those five '
     'are sarin'),
]

CASES += SUPP_CASES_CONF
CASES += SUPP_CASES_STALE
CASES += [
    # 2026-10-27. 본문도 같은 세 숫자를 싣는다. SI 만 검사하면 본문이 조용히
    # 구판으로 남는 경로가 열린다. 이번 세션에 네 번 겪은 실패 양상이다.
    ('본문 3좌표 B3LYP p 오기 (0.0368 대 0.128)',
     'it from 0.275 to 0.0368 and from one valid match to five',
     'it from 0.275 to 0.128 and from one valid match to five'),
    ('본문 3좌표 B3LYP 매칭 수 오기 (five 대 four)',
     'it from 0.275 to 0.0368 and from one valid match to five',
     'it from 0.275 to 0.0368 and from one valid match to four'),
    ('본문 병합공간 값 오기 (0.0221 대 0.0321)',
     'rebuilt the same way, holds\nat 0.0221',
     'rebuilt the same way, holds\nat 0.0321'),
]
SUPP_CASES += SUPP_CASES_DFT
SUPP_CASES += SUPP_CASES_EXTRA
SUPP_CASES += SUPP_CASES_DIPOLE


def run():
    r = subprocess.run([sys.executable, 'verify_consistency.py'],
                       capture_output=True, text=True)
    return r.returncode


def selfcheck():
    """주입 자체가 무효인 케이스를 찾는다 (2026-11-02 신설).

    왜 필요한가. 2026-11-02 에 원고 수치 157건을 갱신하면서 시험 앵커를
    일괄 치환했고, 그 과정에서 두 케이스의 원문과 주입값이 같아졌다. 같은
    글을 같은 글로 바꾸면 검증기는 당연히 통과하므로 시험은 'FAIL 주입했는데
    통과했다' 를 보고하는데, 원인은 검증기가 아니라 시험이다. 그 구분을
    사람이 매번 하게 두지 않는다.
    """
    bad = []
    for name, cases in (('CASES', CASES), ('SUPP_CASES', SUPP_CASES),
                        ('SIDECAR_CASES', SIDECAR_CASES),
                        ('ACS_CASES', ACS_CASES),
                        ('ACS_SUPP_CASES', ACS_SUPP_CASES),
                        ('ACS_CL_CASES', ACS_CL_CASES)):
        for desc, a, b in cases:
            if a == b:
                bad.append(f'  [시험 결함] {name}: {desc} — 원문과 주입값이 '
                           '같다. 주입이 무효다')
    return bad


def main():
    _self = selfcheck()
    if _self:
        print('\n'.join(_self))
        print(f'[test_verify] 시험 자체의 결함 {len(_self)}건. 먼저 고칠 것')
        return 1
    for _p in (TEX, SUPP, SIDECAR, ANNOT, ACS, ACS_SUPP, ACS_CL):
        _guard(_p)
    good = open(TEX).read()
    backup = tempfile.NamedTemporaryFile('w', suffix='.tex', delete=False)
    backup.write(good)
    backup.close()
    bad = 0
    try:
        if run() != 0:
            print('  FAIL 정상 본문이 통과하지 못한다 — 먼저 그것을 고쳐야 한다')
            return 1
        print('  ok   정상 본문 통과')
        # 2026-10-15. 절을 SI 로 옮기면 앵커가 main.tex 에서 사라져 시험이
        # 조용히 SKIP 되었다. 두 파일을 모두 뒤져 앵커가 있는 쪽에 주입한다.
        supp_good = open(SUPP).read()
        for desc, a, b in CASES:
            ra, rb = resolve(a, b, good)
            if ra:
                tgt, src = TEX, good
            else:
                ra, rb = resolve(a, b, supp_good)
                if ra:
                    tgt, src = SUPP, supp_good
                else:
                    _g = find_gen(a, b)
                    if _g is None:
                        print(f'  SKIP {desc}: 원문을 찾지 못함 '
                              f'({a[:40]!r})')
                        bad += 1
                        continue
                    _name, ra, rb = _g
                    _paths = gen_paths(_name)
                    _bak = {q: open(q).read() for q in _paths}
                    try:
                        for q in _paths:
                            open(q, 'w').write(_bak[q].replace(ra, rb, 1))
                        rc = run()
                    finally:
                        for q in _paths:
                            open(q, 'w').write(_bak[q])
                    if rc == 0:
                        print(f'  FAIL {desc}: 주입했는데 통과했다 (생성표)')
                        bad += 1
                    else:
                        print(f'  ok   {desc}: 잡았다 (생성표 {_name})')
                    continue
            open(tgt, 'w').write(src.replace(ra, rb, 1))
            rc = run()
            open(tgt, 'w').write(src)
            where = '' if tgt == TEX else ' (SI)'
            if rc == 0:
                print(f'  FAIL {desc}: 주입했는데 통과했다{where}')
                bad += 1
            else:
                print(f'  ok   {desc}: 잡았다{where}')
        pbak = open(SUPP).read()
        try:
            for desc, a, b in SUPP_CASES:
                ra, rb = resolve(a, b, pbak)
                if not ra:
                    _g = find_gen(a, b)
                    if _g is None:
                        print(f'  SKIP {desc}: SI 원문을 찾지 못함 '
                              f'({a[:40]!r})')
                        bad += 1
                        continue
                    _name, _ra, _rb = _g
                    _paths = gen_paths(_name)
                    _bak = {q: open(q).read() for q in _paths}
                    try:
                        for q in _paths:
                            open(q, 'w').write(_bak[q].replace(_ra, _rb, 1))
                        _rc = run()
                    finally:
                        for q in _paths:
                            open(q, 'w').write(_bak[q])
                    if _rc == 0:
                        print(f'  FAIL {desc}: 주입했는데 통과했다 (생성표)')
                        bad += 1
                    else:
                        print(f'  ok   {desc}: 잡았다 (생성표 {_name})')
                    continue
                open(SUPP, 'w').write(pbak.replace(ra, rb, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다 (SI)')
                open(SUPP, 'w').write(pbak)
        finally:
            open(SUPP, 'w').write(pbak)

        # 2026-11-02 신설. ACS 판본 전용. 세 원논문 Methods 를 읽기 전에
        # 적은 두 서술(ppm, 반복측정)이 모두 수행하지 않은 측정을 가리켰다.
        # ACS 커버레터 점검 (2026-11-02).
        cbak = open(ACS_CL).read()
        try:
            for desc, a, b in ACS_CL_CASES:
                ra, rb = resolve(a, b, cbak)
                if not ra:
                    print(f'  SKIP {desc}: 커버레터 원문을 찾지 못함 ({a[:40]!r})')
                    bad += 1
                    continue
                open(ACS_CL, 'w').write(cbak.replace(ra, rb, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다 (커버레터)')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다 (커버레터)')
                open(ACS_CL, 'w').write(cbak)
        finally:
            open(ACS_CL, 'w').write(cbak)

        # ACS SI 평행 사본 점검 (2026-11-02).
        qbak = open(ACS_SUPP).read()
        try:
            for desc, a, b in ACS_SUPP_CASES:
                ra, rb = resolve(a, b, qbak)
                if not ra:
                    print(f'  SKIP {desc}: ACS SI 원문을 찾지 못함 ({a[:40]!r})')
                    bad += 1
                    continue
                open(ACS_SUPP, 'w').write(qbak.replace(ra, rb, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다 (ACS SI)')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다 (ACS SI)')
                open(ACS_SUPP, 'w').write(qbak)
        finally:
            open(ACS_SUPP, 'w').write(qbak)

        abak = open(ACS).read()
        try:
            for desc, a, b in ACS_CASES:
                ra, rb = resolve(a, b, abak)
                if not ra:
                    print(f'  SKIP {desc}: ACS 원문을 찾지 못함 ({a[:40]!r})')
                    bad += 1
                    continue
                open(ACS, 'w').write(abak.replace(ra, rb, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다 (ACS)')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다 (ACS)')
                open(ACS, 'w').write(abak)
        finally:
            open(ACS, 'w').write(abak)

        sbak = open(SIDECAR).read()
        try:
            for desc, a, b in SIDECAR_CASES:
                ra, rb = resolve(a, b, sbak)
                if not ra:
                    print(f'  SKIP {desc}: 원문을 찾지 못함')
                    bad += 1
                    continue
                open(SIDECAR, 'w').write(sbak.replace(ra, rb, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다')
                open(SIDECAR, 'w').write(sbak)
        finally:
            open(SIDECAR, 'w').write(sbak)
        import json
        abak = open(ANNOT).read()
        try:
            for desc, fig, key, text in ANNOT_CASES:
                d = json.load(open(ANNOT))
                d[fig][key] = text
                json.dump(d, open(ANNOT, 'w'), ensure_ascii=False, indent=1)
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다')
                open(ANNOT, 'w').write(abak)
        finally:
            open(ANNOT, 'w').write(abak)
    finally:
        shutil.copyfile(backup.name, TEX)
    _n = (len(CASES) + len(SUPP_CASES) + len(ANNOT_CASES)
          + len(SIDECAR_CASES) + len(ACS_CASES) + len(ACS_SUPP_CASES)
          + len(ACS_CL_CASES))
    print(f'[test_verify] 시험 {_n}건, 실패 {bad}건')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
