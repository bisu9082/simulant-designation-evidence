#!/usr/bin/env python3
"""verify_consistency.py 가 실제로 결함을 잡는지 시험한다 (2026-10-11 신설).

통과하는 검증기는 그 자체로는 아무것도 보증하지 않는다. 종전 검증기는
핵심 수치가 본문에 '존재하는지'만 보았으므로 늘 통과했고, 같은 양이 208 과
218 로 동시에 적혀 있어도 통과했다. 그래서 과거에 실제로 발생한 결함을
본문에 한 건씩 주입하고, 검증기가 그때마다 FAIL 하는지 확인한다.

여기 실린 여섯 건은 모두 GATE 8 심사에서 지적된 실제 사례다.
"""
import shutil
import subprocess
import sys
import tempfile

TEX = 'paper/main.tex'

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
     'value of 0.026 sits at the 24.3rd percentile',
     'value of 0.026 sits at the 4.2nd percentile'),
    ('주성분3 백분위 오기 (2.1 대 0.4)',
     'value of 0.026 sits at the 24.3rd',
     'value of 0.026 sits at the 0.4th'),
    ('주 패널 유효매칭 수 오기 (zero 대 three)',
     'twenty-eight-compound panel and zero on the primary one',
     'twenty-eight-compound panel and three on the primary one'),
    ('28종 유효매칭 수 오기 (three 대 five)',
     'against\nfour on the twenty-eight-compound panel',
     'against\nseven on the twenty-eight-compound panel'),
    ('WY 조정 p 오기 (0.075 대 0.045) — 유의/비유의가 뒤집힌다',
     '& 0.169 & 0.038 & 0.025', '& 0.169 & 0.038 & 0.045'),
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
     'would survive at $q=0.040$',
     'would survive at $q=0.142$'),
    ('백색화 통제 수치 오기 (16.9 대 36.9)',
     'the rate falls from 33.7\\% to 18.4\\%',
     'the rate falls from 33.7\\% to 38.4\\%'),
    ('meff 오기 (5.00 대 15.00)',
     'estimates\n5.00 effective tests',
     'estimates\n15.00 effective tests'),
    ('표 1 의 28종 전자공간 행이 구판 (0.275 대 0.297)',
     '& electronic  &  3 & 0.275  &',
     '& electronic  &  3 & 0.297  &'),
    # --- 네 번째 구성 (2026-10-20). 이 논문에서 다중검정 보정을 통과하는
    # 유일한 구성이므로 수치 하나가 틀리면 주장이 통째로 바뀐다.
    ('표 1 네 번째 블록 구조공간 행이 틀림 (0.013 대 0.113)',
     '28 cpd., incumbent 3   & structural  & 37 & \\textbf{0.0265}',
     '28 cpd., incumbent 3   & structural  & 37 & \\textbf{0.1265}'),
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
    ('네 번째 구성 순열 수 오기 (1,961,256 대 1,916,256)',
     '\\binom{24}{10}=1\\,961\\,256',
     '\\binom{24}{10}=1\\,916\\,256'),
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
SIDECAR = 'paper/cover_letter.tex'
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
     '0.026 sits at the 24.3rd percentile',
     '0.028 sits at the 24.3rd percentile'),
    ('부분공간 백분위를 구판으로 (25.4 대 22.7)',
     '0.026 sits at the 24.3rd percentile of that\ndistribution',
     '0.026 sits at the 21.3rd percentile of that\ndistribution'),
    ('도판2 캡션 백분위 오기 (26th 대 21st)',
     '0.026 is the 24th percentile of those draws',
     '0.026 is the 21th percentile of those draws'),
    ('28종 매치 수를 구판으로 (four 대 three)',
     'Omitting those five restores four matches',
     'Omitting those five restores three matches'),
    ('구조공간 유의 통계량 수를 구판으로 (two 대 one)',
     'compounds under two of six statistics in the structural space',
     'compounds under one of six statistics in the structural space'),
    ('Westfall-Young 을 구판으로 (0.079 대 0.075)',
     'rather than approximating it, gives 0.078.',
     'rather than approximating it, gives 0.075.'),
    ('사전지정 q 를 구판으로 (0.046 대 0.042)',
     'merged results would survive at $q=0.040$',
     'merged results would survive at $q=0.042$'),
    ('도판1 캡션 BH 를 구판으로 (0.185 대 0.159)',
     'value over its eighteen tests is 0.165.}',
     'value over its eighteen tests is 0.159.}'),
    ('28종 구조좌표 수를 구판으로 (36 대 38)',
     'the observed 36-dimensional value of',
     'the observed 38-dimensional value of'),
    ('두 패널 보존 좌표를 같게 (38/36 대 38/38)',
     'leaving 38 and 36 structural coordinates',
     'leaving 38 and 38 structural coordinates'),
    ('SI tabun 거리를 구판으로 (1.989 대 2.003)',
     'its normalised distance is 1.995 against a threshold of 2',
     'its normalised distance is 2.003 against a threshold of 2'),
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
     'sits between 1.896 and\n1.905 across the ensemble',
     'sits between 1.796 and\n1.905 across the ensemble'),
    ('앙상블 tabun 상한 오기 (1.905 대 1.805)',
     'sits between 1.896 and\n1.905 across the ensemble',
     'sits between 1.896 and\n1.805 across the ensemble'),
    ('단일선택 tabun 범위 오기 (1.876 대 1.976)',
     'between 1.876 and 1.884 across the',
     'between 1.776 and 1.884 across the'),
    ('단일선택 순열 p 상한 오기 (0.028 대 0.038)',
     '0.025 and 0.028\nin the single-conformer runs',
     '0.025 and 0.038\nin the single-conformer runs'),
    ('시드 수 오기 (six 대 four)',
     'under six seeds and rebuilding',
     'under four seeds and rebuilding'),
    ('구조좌표 수를 메타열 포함 시절로 (36 대 37)',
     'value of 0.026 sits at the 24.3rd',
     'value of 0.029 sits at the 24.3rd'),
    ('구조공간 명목 p 를 앙상블 이전으로 (0.026 대 0.031)',
     '0.026 sits at the 24.3rd percentile',
     '0.031 sits at the 24.3rd percentile'),
    ('유효검정수를 구판으로 (5.00 대 7.00)',
     '5.00 effective tests of eighteen',
     '7.00 effective tests of eighteen'),
]

CASES += SUPP_CASES_CONF
CASES += SUPP_CASES_STALE
SUPP_CASES += SUPP_CASES_DFT
SUPP_CASES += SUPP_CASES_EXTRA


def run():
    r = subprocess.run([sys.executable, 'verify_consistency.py'],
                       capture_output=True, text=True)
    return r.returncode


def main():
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
            if a in good:
                tgt, src = TEX, good
            elif a in supp_good:
                tgt, src = SUPP, supp_good
            else:
                print(f'  SKIP {desc}: 원문을 찾지 못함 ({a[:40]!r})')
                bad += 1
                continue
            open(tgt, 'w').write(src.replace(a, b, 1))
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
                if a not in pbak:
                    print(f'  SKIP {desc}: SI 원문을 찾지 못함 ({a[:40]!r})')
                    bad += 1
                    continue
                open(SUPP, 'w').write(pbak.replace(a, b, 1))
                if run() == 0:
                    print(f'  FAIL {desc}: 주입했는데 통과했다')
                    bad += 1
                else:
                    print(f'  ok   {desc}: 잡았다 (SI)')
                open(SUPP, 'w').write(pbak)
        finally:
            open(SUPP, 'w').write(pbak)

        sbak = open(SIDECAR).read()
        try:
            for desc, a, b in SIDECAR_CASES:
                if a not in sbak:
                    print(f'  SKIP {desc}: 원문을 찾지 못함')
                    bad += 1
                    continue
                open(SIDECAR, 'w').write(sbak.replace(a, b, 1))
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
    print(f'[test_verify] 시험 {len(CASES)+len(SUPP_CASES)+len(ANNOT_CASES)+len(SIDECAR_CASES)}건, 실패 {bad}건')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
