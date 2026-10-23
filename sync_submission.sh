#!/bin/sh
# 투고 패키지 동기화 (2026-10-13 신설).
#
# submission/ 은 손으로 채워 왔고, 그래서 작업본과 어긋난 사본이 들어간 적이 있다.
# 원고·도판·코드·산출물을 한 번에 복사하고, 복사 뒤 코드 사본으로 검증기를
# 다시 돌려 패키지 자체가 자기 정합적인지 확인한다.
set -e
cd "$(dirname "$0")"

mkdir -p submission/code submission/outputs

cp -p paper/main.tex paper/refs.bib paper/cover_letter.tex \
      paper/highlights.tex paper/supplementary.tex submission/
for f in main.pdf cover_letter.pdf highlights.pdf supplementary.pdf; do
  [ -f "paper/$f" ] && cp -p "paper/$f" "submission/$f"
done
[ -f paper/main.pdf ] && cp -p paper/main.pdf submission/simulant_manuscript.pdf

for f in Fig_main.png Fig_spaces.png Fig_toxicity.png Fig_plocal.png \
         Fig_h1.png Fig_mde.png graphical_abstract.png graphical_abstract.pdf; do
  cp -p "outputs/$f" "submission/$f"
done

# 코드: 분석·도판·검증 스크립트와 빌드 스크립트 전량
for f in *.py *.sh; do
  case "$f" in sync_submission.sh) continue ;; esac
  cp -p "$f" submission/code/
done
cp -p descriptor_cache*.csv submission/code/ 2>/dev/null || true

# 산출물: 본문이 인용하는 수치의 출처
for f in experiment_summary.json permutation_sidedness.json rollback_c.json \
         multiplicity_c.json match_by_panel.json deltaE_premise.json \
         figure_annotations.json sanity_checks.json; do
  [ -f "outputs/$f" ] && cp -p "outputs/$f" submission/outputs/
done
for f in outputs/*.csv; do
  [ -f "$f" ] && cp -p "$f" submission/outputs/
done

n_code=$(ls submission/code | wc -l)
n_out=$(ls submission/outputs | wc -l)
printf '[sync] code %s개, outputs %s개 복사\n' "$n_code" "$n_out"
