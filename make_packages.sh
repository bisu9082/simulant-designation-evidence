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
printf '[pkg] %s\n' "$OV"
unzip -l "$OV" | sed -n '4,$p' | head -20

# 저장소 사본: 코드·원고·산출물 전량 (submission/ 이 이미 자기검증을 통과한 것)
zip -q -X -r "$RP" submission
printf '[pkg] %s  (%s개 항목)\n' "$RP" "$(unzip -l "$RP" | tail -1 | awk '{print $2}')"
