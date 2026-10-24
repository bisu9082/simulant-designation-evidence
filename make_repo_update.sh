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
zip -q -X "$Z" descriptor_cache.csv descriptor_cache_35.csv \
    descriptor_cache_xtb.csv descriptor_cache_xtb_35.csv \
    inchikeys.csv panel_all.csv panel_new.csv
zip -q -X "$Z" paper/main.tex paper/supplementary.tex paper/cover_letter.tex \
    paper/highlights.tex paper/abstract_300.tex paper/refs.bib
zip -q -X "$Z" paper/Fig_*.png paper/main.pdf paper/supplementary.pdf \
    paper/cover_letter.pdf paper/highlights.pdf
# outputs/ 중 저장소에 들어갈 것만. _crop_* 는 도판 검수용 조각이고
# _superseded/ 는 폐기본이다. 2026-10-24 에 둘 다 패키지로 새어 나갔다.
zip -q -X -r "$Z" outputs -i 'outputs/*.json' 'outputs/*.csv' 'outputs/*.png' \
    'outputs/*.pdf' 'outputs/*.tex' \
    -x 'outputs/_crop_*' 'outputs/_superseded/*'
zip -q -X -r "$Z" sources

n=$(unzip -l "$Z" | tail -1 | awk '{print $2}')
printf '[repo] %s  (%s개 파일, %s)\n' "$Z" "$n" "$(du -h "$Z" | cut -f1)"
