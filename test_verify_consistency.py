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
     '(Spearman $r=+0.004$',
     'independent (Spearman $r=-0.022$'),
    # --- 롤백 C 헤드라인 (2026-10-13 추가). 이 수치들이 주장을 떠받친다.
    ('BH 보정값 오기 (0.159 대 0.059)',
     'adjusted value of 0.159',
     'adjusted value of 0.059'),
    ('무작위 3종 부분공간 비율 오기 (36.3 대 6.3)',
     'and 36.3\\% of the subspaces reach',
     'and 6.3\\% of the subspaces reach'),
    ('주성분 3개 절단 p 오기 (0.0044 대 0.0440)',
     'from $p=0.028$ to\n$p=0.0044$',
     'from $p=0.028$ to\n$p=0.0440$'),
    ('Mantel p 를 쌍독립 가정 p 로 바꿔 씀 (0.012 대 0.038)',
     '$\\rho = +0.539$, Mantel $p = 0.012$',
     '$\\rho = +0.539$, Mantel $p = 0.038$'),
    ('omega 상관계수 오기 (-0.899 대 -0.959)',
     'r(E_{\\mathrm{LUMO}}, \\omega) =\n-0.899',
     'r(E_{\\mathrm{LUMO}}, \\omega) =\n-0.959'),
    ('GA 코사인 오기 (0.764 대 0.564)',
     'tabun, at $+0.764$ against sarin',
     'tabun, at $+0.564$ against sarin'),
    ('백분위 오기 — 구조 관측값을 더 특별해 보이게 (25.8 대 4.2)',
     'value of 0.028 sits at the 25.8th percentile',
     'value of 0.028 sits at the 4.2th percentile'),
    ('주성분3 백분위 오기 (2.1 대 0.4)',
     'value of 0.0044 sits at the 2.1st',
     'value of 0.0044 sits at the 0.4st'),
    ('주 패널 유효매칭 수 오기 (zero 대 three)',
     'twenty-eight-compound panel and zero on the primary one',
     'twenty-eight-compound panel and three on the primary one'),
    ('28종 유효매칭 수 오기 (three 대 five)',
     'against\nthree on the twenty-eight-compound panel',
     'against\nfive on the twenty-eight-compound panel'),
    ('WY 조정 p 오기 (0.075 대 0.045) — 유의/비유의가 뒤집힌다',
     '& 0.402 & 0.113 & 0.075', '& 0.402 & 0.113 & 0.045'),
    ('부트스트랩 0건 비율 오기 (58.8 대 98.8)',
     'exactly zero\nin 58.8\\% of replicates',
     'exactly zero\nin 98.8\\% of replicates'),
    ('부트스트랩 차 양수 복제수 오기 (1000 대 600)',
     'positive in 1000 of 1000 replicates',
     'positive in 600 of 1000 replicates'),
    ('GB 회복률 오기 (37.8 대 7.8) — 가장 약한 지점을 숨김',
     'as its nearest compound in 37.8\\% of replicates',
     'as its nearest compound in 7.8\\% of replicates'),
    ('35종 귀무기대 오기 (9.63 대 6.63)',
     'expectation of 9.63 set by the 21:3',
     'expectation of 6.63 set by the 21:3'),
    ('최소검출효과 오기 (0.22 대 0.42) — 검정력 주장이 뒤집힌다',
     'moved 0.22 typical nearest-neighbour',
     'moved 0.42 typical nearest-neighbour'),
    ('실제 결손 중위 오기 (0.463 대 0.163) — 검출 가능 폭보다 작아진다',
     'carry has a median of\n0.463',
     'carry has a median of\n0.163'),
    ('사린 결손 오기 (0.117 대 0.317)',
     'only at sarin, where it is 0.117',
     'only at sarin, where it is 0.317'),
    ('CEES-HD 거리 오기 (0.57 대 1.57)',
     'sulfide at 0.57 from sulfur mustard',
     'sulfide at 1.57 from sulfur mustard'),
    ('DFP 변위 거리 오기 (1.07 대 2.07) — 변위 논증이 무너진다',
     'fluorophosphate is present, at 1.07, 0.99 and 1.17',
     'fluorophosphate is present, at 2.07, 0.99 and 1.17'),
    ('사전지정 q 오기 (0.042 대 0.142) — 유의/비유의가 뒤집힌다',
     'would survive at $q=0.042$',
     'would survive at $q=0.142$'),
    ('백색화 통제 수치 오기 (16.9 대 36.9)',
     'from 36.3\\% to 16.9\\%',
     'from 36.3\\% to 36.9\\%'),
    ('meff 오기 (5.00 대 15.00)',
     'estimates 5.00 effective tests',
     'estimates 15.00 effective tests'),
    ('표 1 의 28종 전자공간 행이 구판 (0.275 대 0.297)',
     '& electronic  &  3 & 0.275 &',
     '& electronic  &  3 & 0.297 &'),
    # --- 네 번째 구성 (2026-10-20). 이 논문에서 다중검정 보정을 통과하는
    # 유일한 구성이므로 수치 하나가 틀리면 주장이 통째로 바뀐다.
    ('표 1 네 번째 블록 구조공간 행이 틀림 (0.013 대 0.113)',
     '35 cpd., incumbent 3 & structural  & 38 & \\textbf{0.013}',
     '35 cpd., incumbent 3 & structural  & 38 & \\textbf{0.113}'),
    ('표 1 네 번째 블록 병합공간 행이 틀림 (0.005 대 0.055)',
     '+ pesticides 2 & merged      & 41 & \\textbf{0.005}',
     '+ pesticides 2 & merged      & 41 & \\textbf{0.055}'),
    ('네 번째 구성 BH q 를 비유의로 오기 (0.0006 대 0.0600) — 결론이 뒤집힌다',
     'Benjamini-Hochberg value of 0.0006 over the same',
     'Benjamini-Hochberg value of 0.0600 over the same'),
    ('표 2 네 번째 열의 사전지정 값 오기 (0.015 대 0.150) — 유의성이 뒤집힌다',
     'Statistic specified in advance (3 tests) & 0.386 & 0.092 & \\textbf{0.042} & \\textbf{0.015}',
     'Statistic specified in advance (3 tests) & 0.386 & 0.092 & \\textbf{0.042} & \\textbf{0.150}'),
    ('표 2 네 번째 열의 meff 오기 (5.00 대 15.00)',
     '$m_{\\mathrm{eff}}$ of 18 (Li-Ji)     & 6.00  & 5.00  & 5.00  & 5.00',
     '$m_{\\mathrm{eff}}$ of 18 (Li-Ji)     & 6.00  & 5.00  & 5.00  & 15.00'),
    ('서술자 제거 수 오기 (13 대 23) — 51 과 38 이 안 맞게 된다',
     'removes 13 of the 51',
     'removes 23 of the 51'),
    ('초록의 네 번째 구성 보정값 오기 (0.0006 대 0.0060)',
     'produces an adjusted\nvalue of 0.0006.',
     'produces an adjusted\nvalue of 0.0060.'),
    ('결론의 네 번째 구성 보정값 오기 (0.0006 대 0.0006 유지 실패)',
     'into one with an adjusted value of\n0.0006.',
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
     'lowers that value to 0.076',
     'lowers that value to 0.176'),
    # --- 표 3 (2026-10-17 추가). 이 표는 그때까지 검증기에 등록되지 않아
    # 11행 중 9행이 구판 값인 채로 144건 검사를 통과했다. 동료 검토에서
    # 지적됐다. 반올림으로 설명되지 않던 두 행과 산문·표 불일치를 시험한다.
    ('표 3 PS 행이 구판 (2.530 대 2.70) — 반올림으로 설명 안 됨',
     'PS & DMMP & 2.530 & no',
     'PS & DMMP & 2.70 & no'),
    ('표 3 CG 행이 구판 (1.600 대 1.64)',
     'CG & DCM & 1.600 & no',
     'CG & DCM & 1.64 & no'),
    ('표 3 GA 행이 임계값 2 를 넘지 않게 보이도록 오기',
     'GA & DMMP & 2.003 & no',
     'GA & DMMP & 1.903 & no'),
    ('산문이 표 3 과 다른 소만 거리를 싣는다 (1.417 대 1.40)',
     'soman (1.417) and',
     'soman (1.40) and'),
    ('전자공간 유효매칭 개수 주장이 산문에만 있고 거리는 구판 (1.407 대 1.335)',
     'dichloromethane at 1.407,',
     'dichloromethane at 1.335,'),
    ('포스겐 경쟁 거리를 다른 공간의 값으로 섞어 씀 (1.12 대 1.14)',
     'ethyl acetate in the merged space (1.12 against',
     'ethyl acetate in the merged space (1.14 against'),
    # 2026-10-23. 역방향 최근접 작용제 부트스트랩.
    ('DPGME 의 G계열 비율이 구판이다 (46.7 대 61.2)',
     'the\nG series in 46.7\\%', 'the\nG series in 61.2\\%'),
    ('DMMP 의 사린 최근접 비율 오기 (93.3 대 88.1)',
     'DMMP is nearest sarin in\n93.3\\% of replicates',
     'DMMP is nearest sarin in\n88.1\\% of replicates'),
    ('DCM 의 G계열 비율 오기 (1.4 대 12.0)',
     'Dichloromethane is nearest a G-series agent in 1.4\\%',
     'Dichloromethane is nearest a G-series agent in 12.0\\%'),
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
SIDECAR_CASES = [
    ('커버레터가 본문과 다른 Spearman 을 싣는다',
     '$r = +0.004$ (95\\% CI $-0.132$ to $+0.143$)',
     '$r = -0.022$ (95\\% CI $-0.156$ to $+0.111$)'),
]

# 도판 주석 결함은 outputs/figure_annotations.json 에 주입한다.
ANNOT = 'outputs/figure_annotations.json'
ANNOT_CASES = [
    ('도판 주석이 데이터가 반증하는 주장을 한다',
     'Fig_toxicity', 'panel_a',
     'dichloromethane, a designated simulant, is the least acutely toxic '
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
     'Twenty-eight, incumbent three & 3/11 & 0.025 & 0.009 & 0.159 & 0.132',
     'Twenty-eight, incumbent three & 3/11 & 0.025 & 0.009 & 0.159 & 0.032'),
    ('블록균등에서 유효매칭이 변한 것처럼 오기 (0/11 대 2/11)',
     'Thirty-five, incumbent three & 0/11',
     'Thirty-five, incumbent three & 2/11'),
]


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
