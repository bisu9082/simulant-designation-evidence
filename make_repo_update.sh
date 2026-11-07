#!/bin/sh
# 저장소 갱신본 묶기 (2026-10-24 신설).
#
# 종전에는 바뀐 파일을 손으로 골라 담았고, 그래서 목록이 어긋났다
# (삭제한 Fig_panel 이 계속 들어가고 새 Fig_plocal 은 빠졌다).
# 저장소에 있어야 할 것을 통째로 담는다. 풀어 덮은 뒤 git add -A 로 올린다.
set -e
cd "$(dirname "$0")"
Z="simulant_repo_update_$(date +%Y%m%d).zip"
rm -f "$Z"

zip -q -X "$Z" *.py *.sh
# 2026-11-07 (10차 R5-B4). README 와 requirements.txt 가 빠져 있었다.
# README 의 첫 명령이 pip install -r requirements.txt 인데 그 파일이
# 패키지에 없었다.
zip -q -X "$Z" README.md requirements.txt
# 2026-11-07. 서술자 캐시 목록이 불완전했다. match_by_panel_ascl3.py 와
# 서술자 관문이 읽는 세 파일이 빠져 있어, 푼 패키지에서 관문이 조용히
# 건너뛰고도 '전 수치열 재현' 을 찍었다.
zip -q -X "$Z" descriptor_cache.csv descriptor_cache_35.csv \
    descriptor_cache_xtb.csv descriptor_cache_xtb_35.csv \
    descriptor_cache_xtb_36.csv descriptor_cache_ascl3.csv \
    descriptor_cache_xtb_ascl3.csv \
    inchikeys.csv panel_all.csv panel_new.csv panel_ascl3.csv
# 2026-11-07 (10차 R5-B4). deltaE/ 는 원고가 'the deposited Delta E_00
# table' 이라고 가리키는 측정 입력인데 패키지에 한 건도 없었다. 푼 패키지
# 에서 deltaE_premise.py 가 첫 줄에서 죽는다.
zip -q -X "$Z" deltaE/*.csv
# 2026-11-04 (J-C). ACS본이 정본. JHM 평행본은 archive/jhm_v1/ 동결.
zip -q -X "$Z" paper_acs/main.tex paper_acs/supporting.tex \
    paper_acs/cover_letter.tex paper_acs/refs.bib paper_acs/tab_*.tex
zip -q -X "$Z" paper_acs/Fig_*.png paper_acs/toc_acs.png \
    paper_acs/main.pdf paper_acs/supporting.pdf paper_acs/cover_letter.pdf
# outputs/ 중 저장소에 들어갈 것만. _crop_* 는 도판 검수용 조각이고
# _superseded/ 는 폐기본이다. 2026-10-24 에 둘 다 패키지로 새어 나갔다.
zip -q -X -r "$Z" outputs -i 'outputs/*.json' 'outputs/*.csv' 'outputs/*.png' \
    'outputs/*.pdf' 'outputs/*.tex' \
    -x 'outputs/_crop_*' 'outputs/_superseded/*'
# 2026-11-06 (9차). 종전에는 sources/ 전체를 넣었다. 그 안에는 투고 전 내부
# 자기심사 메모(GATE8 채점·교정계획, T1/T2 작업결과, 원논문 Methods 직독)가
# 들어 있어서 공개 저장소에 함께 올라갔다. 원고가 예치를 약속한 것은 파싱한
# 규제값과 그 출처뿐이므로 1차 사료만 넣는다. 새 1차 사료를 받으면 여기에
# 확장자나 파일명을 더할 것. 와일드카드를 열어 두면 메모가 다시 샌다.
zip -q -X "$Z" sources/compiled_aegls_27jul2018.pdf \
    sources/registers_20261106.csv \
    sources/designations_20261106.csv \
    sources/dye_panel_20261106.csv

n=$(unzip -l "$Z" | tail -1 | awk '{print $2}')
printf '[repo] %s  (%s개 파일, %s)\n' "$Z" "$n" "$(du -h "$Z" | cut -f1)"
