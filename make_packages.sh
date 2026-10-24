#!/bin/sh
# 투고 패키지 묶기 (2026-10-23 신설).
#
# 종전에는 zip 을 손으로 만들었고, 그래서 삭제된 도판(Fig_panel, Fig_qsar_failure)
# 이 남고 새 도판(Fig_plocal)이 빠진 사본이 나갔다. 원고가 실제로 \includegraphics
# 하는 파일만 골라 담고, 담은 목록을 출력한다.
set -e
cd "$(dirname "$0")"
D=$(date +%Y%m%d)
OV="simulant_overleaf_${D}.zip"
RP="simulant_repo_${D}.zip"

FIGS=$(grep -ho 'includegraphics\[[^]]*\]{[^}]*}' paper/main.tex paper/supplementary.tex \
       | sed 's/.*{\(.*\)}/\1/' | sort -u)
for f in $FIGS; do
  [ -f "paper/$f" ] || { echo "원고가 부르는 도판이 없다: $f"; exit 1; }
done

rm -f "$OV" "$RP"
( cd paper && zip -q -X "../$OV" main.tex supplementary.tex cover_letter.tex \
    highlights.tex refs.bib $FIGS )
# 컴파일 순서 안내. main 과 supplementary 가 xr 로 서로를 참조하므로 한 번씩만
# 돌리면 양쪽에 ?? 가 남는다. 받는 사람이 그걸 모르면 투고본이 깨져 보인다.
if [ -f README_overleaf.txt ]; then
  zip -q -X "$OV" README_overleaf.txt
fi
printf '[pkg] %s\n' "$OV"
unzip -l "$OV" | sed -n '4,$p' | head -20

# 저장소 사본: 코드·원고·산출물 전량 (submission/ 이 이미 자기검증을 통과한 것)
zip -q -X -r "$RP" submission
printf '[pkg] %s  (%s개 항목)\n' "$RP" "$(unzip -l "$RP" | tail -1 | awk '{print $2}')"
