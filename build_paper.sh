#!/bin/sh
# 원고 빌드 정본 (2026-10-13 신설).
#
# 종전에는 도판을 outputs/ 에 만들고 paper/ 로 손으로 복사했다. 그래서 분석을
# 다시 돌려도 원고에 박히는 PNG 는 옛것인 경우가 생겼고, GATE 8 3차에서 도판과
# 본문이 어긋난 지적이 나온 경로가 이것이다. 복사와 컴파일을 한 곳에 묶는다.
#
# 사용: sh build_paper.sh [--figs]
#   --figs 를 주면 도판을 다시 생성한다. 분석을 돌린 뒤에는 반드시 붙일 것.
set -e
cd "$(dirname "$0")"

if [ "$1" = "--figs" ]; then
  for s in make_fig_main.py make_fig_toxicity.py make_fig_h1.py \
           make_fig_electronic.py make_fig_mde.py make_ga.py; do
    echo "[fig] $s"
    python3 "$s" > /dev/null
  done
fi

for f in Fig_main.png Fig_spaces.png Fig_toxicity.png Fig_plocal.png \
         Fig_h1.png Fig_mde.png; do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "paper/$f"
  else
    echo "[경고] outputs/$f 없음 — 원고의 도판이 갱신되지 않는다" >&2
  fi
done

cd paper
# 본문과 Supplementary 가 xr 로 서로를 참조하므로 교대로 두 번씩 돌린다.
# 한 번만 돌리면 상대의 .aux 가 아직 없어 ?? 가 박힌다.
for pass in 1 2; do
  pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode supplementary.tex > /dev/null 2>&1 || true
done
bibtex main > /dev/null 2>&1 || true
bibtex supplementary > /dev/null 2>&1 || true
for pass in 1 2; do
  pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1 || true
  pdflatex -interaction=nonstopmode supplementary.tex > /dev/null 2>&1 || true
done
for aux in cover_letter highlights; do
  pdflatex -interaction=nonstopmode "$aux.tex" > /dev/null 2>&1 || true
done

for d in main supplementary; do
  printf '%-14s 쪽수 %s  미정의참조 %s  미정의인용 %s  Overfull %s\n' "$d" \
    "$(pdfinfo $d.pdf 2>/dev/null | awk '/^Pages/{print $2}')" \
    "$(grep -c 'Reference.*undefined' $d.log || true)" \
    "$(grep -c 'Citation.*undefined' $d.log || true)" \
    "$(grep -c 'Overfull' $d.log || true)"
done
