#!/bin/sh
# 투고 패키지 동기화 (2026-10-13 신설).
#
# submission/ 은 손으로 채워 왔고, 그래서 작업본과 어긋난 사본이 들어간 적이 있다.
# 원고·도판·코드·산출물을 한 번에 복사하고, 복사 뒤 코드 사본으로 검증기를
# 다시 돌려 패키지 자체가 자기 정합적인지 확인한다.
set -e
cd "$(dirname "$0")"

# 2026-11-07 (11차). code/ 와 outputs/ 를 비우고 다시 채운다. 종전에는
# 덮어쓰기만 했으므로, 루트에서 사라진 스크립트의 사본이 패키지에 계속
# 남았다. make_ga.py 를 archive/jhm_v1/ 로 옮겼는데 submission/code/ 에는
# 그대로 있었고, 그런 사본은 아무도 갱신하지 않는다.
rm -rf submission/code submission/outputs
mkdir -p submission/code submission/outputs

# 2026-11-04 (결정 J-C). 투고 패키지는 ACS Sensors 판본이다. JHM 평행본은
# archive/jhm_v1/ 로 동결했다 (사유와 되살리는 절차는 그쪽 FROZEN.md).
cp -p paper_acs/main.tex paper_acs/refs.bib paper_acs/cover_letter.tex \
      paper_acs/supporting.tex submission/
for f in acs-main.bib acs-supporting.bib; do
  [ -f "paper_acs/$f" ] && cp -p "paper_acs/$f" "submission/$f"
done
for f in main.pdf cover_letter.pdf supporting.pdf; do
  [ -f "paper_acs/$f" ] && cp -p "paper_acs/$f" "submission/$f"
done
[ -f paper_acs/main.pdf ] && cp -p paper_acs/main.pdf \
    submission/simulant_manuscript.pdf

for f in Fig_main.png Fig_plocal.png Fig_h1.png Fig_mde.png \
         Fig_talk_response.png Fig_talk_tox.png toc_acs.png \
         graphical_abstract.png graphical_abstract.pdf; do
  [ -f "outputs/$f" ] && cp -p "outputs/$f" "submission/$f"
done

# 코드: 분석·도판·검증 스크립트와 빌드 스크립트 전량
for f in *.py *.sh; do
  case "$f" in sync_submission.sh) continue ;; esac
  cp -p "$f" submission/code/
done
cp -p descriptor_cache*.csv submission/code/ 2>/dev/null || true
# 2026-11-07 (11차). README.md 와 requirements.txt 가 빠져 있었다. README 는
# 재현 지시서이고 검증기가 그 안의 재현 약속 수치를 대조하므로(키 rm_*),
# 없으면 패키지 자기검증이 '검사 누락 6건' 으로 떨어진다. 공개 zip 은
# 10차에 이미 둘을 담게 고쳤는데 이쪽만 빠져 있었다.
cp -p README.md requirements.txt submission/code/

# 산출물: 본문이 인용하는 수치의 출처. 2026-10-23 이전에는 여덟 개만 손으로
# 골라 복사했고, 그 목록이 갱신되지 않아 h1_summary.json 을 비롯해 헤드라인
# 수치의 출처 열세 개가 패키지에서 빠져 있었다. 전량 복사로 바꾼다.
for f in outputs/*.json; do
  [ -f "$f" ] && cp -p "$f" submission/outputs/
done
for f in outputs/*.csv; do
  [ -f "$f" ] && cp -p "$f" submission/outputs/
done

n_code=$(ls submission/code | wc -l)
n_out=$(ls submission/outputs | wc -l)
printf '[sync] code %s개, outputs %s개 복사\n' "$n_code" "$n_out"

# 패키지 자기 정합성. 이 절은 이 스크립트의 머리말이 처음부터 약속했으나
# 실제로는 없었다 (2026-10-23 적발). 복사된 코드 사본으로 검증기를 돌려
# 패키지 안의 원고와 산출물이 서로 맞는지 확인한다.
( cd submission/code \
  && ln -sfn ../../outputs outputs 2>/dev/null \
  && mkdir -p paper_acs \
  && cp -p ../main.tex ../supporting.tex ../cover_letter.tex paper_acs/ \
  && cp -p ../../paper_acs/tab_*.tex paper_acs/ \
  && cp -p ../Fig_*.png ../toc_acs.png paper_acs/ 2>/dev/null \
  && python3 verify_consistency.py ) || {
  echo '[sync] 패키지 자기검증 실패'; exit 1; }
rm -rf submission/code/paper_acs submission/code/outputs
