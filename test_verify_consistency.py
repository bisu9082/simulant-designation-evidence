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
     'independent across the 207 compounds',
     'independent across the 218 compounds'),
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
    ('도판 캡션의 각주 유래 개수 (5 대 9)',
     'Five of the nineteen values come from the',
     'Nine of the nineteen values come from the'),
    ('적용영역 내 최대오차를 전체 최대로 오기 (2.68 대 3.42)',
     'wrong by 2.68 log units. The dotted line',
     'wrong by 3.42 log units. The dotted line'),
    ('충실도 비율 오기 (3.1 대 3.6)',
     'trichloroethylene, is 3.1 times closer',
     'trichloroethylene, is 3.6 times closer'),
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
