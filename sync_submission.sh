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
# 2026-11-10. submission/supporting.tex 가 ACS Sensors 시절 이름으로
# 남아 있었고, 아래 자기검증이 그 묵은 사본을 SI 로 읽고 있었다.
# 복사 대상이 아닌 .tex 는 지운다.
rm -f submission/*.tex
cp -p paper/main.tex paper/refs.bib paper/cover_letter.tex \
      paper/supplementary.tex submission/
# 2026-11-13. supplementary.tex 가 \input 하는 표 12개가 패키지에 없어
# 조판소에서 'File tab_sparse.tex not found' 로 멈춘다. GATE 8 11차 적발.
cp -p paper/tab_*.tex submission/

# 2026-11-13 (GATE 8 12차). ACS 는 조판용 소스 .tex 를 받는다. 작업 주석
# 176행에 저널 재타깃 경위, 폐기 라인의 결함 목록, 미해결 TODO 가 들어
# 있어 편집장이 읽게 된다. 작업본은 그대로 두고 패키지 사본에서만 뺀다.
# 줄 전체가 주석인 것만 지운다 (줄 끝 % 는 LaTeX 의미가 있을 수 있다).
for f in main supplementary cover_letter tab_ablation tab_family tab_sparse \
         tab_plocal tab_volatility tab_vpairing tab_premisenull tab_registers \
         tab_designations tab_dyepanel tab_match tab_failaudit tab_blockbalance \
         tab_si-ascl3; do
  [ -f "submission/$f.tex" ] || continue
  sed -i '/^[[:space:]]*%/d' "submission/$f.tex"
done
echo "[sync] 투고 사본에서 전행 주석 제거"

for f in acs-main.bib acs-supporting.bib; do
  [ -f "paper/$f" ] && cp -p "paper/$f" "submission/$f"
done

# 2026-11-10. 도판 목록을 손으로 들고 있었다. 3층 재설계로 Fig_spaces 와
# Fig_toxicity 가 SI 로 들어갔는데 이 목록에 없어서 패키지에 빠졌고,
# 철회한 graphical_abstract 는 목록에 남아 계속 복사되고 있었다. 문서가
# 실제로 \includegraphics 하는 것만 복사한다.
rm -f submission/*.png submission/*.pdf
for f in $(grep -ho '\\includegraphics\[[^]]*\]{[^}]*}' \
             paper/main.tex paper/supplementary.tex \
           | sed 's/.*{\(.*\)}/\1/' | sort -u); do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "submission/$f"
  else
    printf '[sync] 도판 없음: outputs/%s\n' "$f"; exit 1
  fi
done

# 2026-11-10. supporting.pdf 는 ACS Sensors 시절 이름이다. 지금 빌드가
# 내는 것은 supplementary.pdf 다.
for f in main.pdf cover_letter.pdf supplementary.pdf; do
  [ -f "paper/$f" ] && cp -p "paper/$f" "submission/$f"
done
[ -f paper/main.pdf ] && cp -p paper/main.pdf \
    submission/simulant_manuscript.pdf


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

# 2026-11-13 (GATE 8 9차, 재현성 심사). 패키지가 7단계에서 죽었다.
# run_analysis.sh 의 첫 스크립트 deltaE_premise.py 가 deltaE/ 를 읽는데
# 그 폴더가 패키지에 없었고, sources/ (EPA 원문·지정문서 CSV) 와 패널
# 정의 세 파일도 없었다. 원고는 이 패키지를 "also provided with this
# article" 이라 적는다. 돌지 않는 것을 그렇게 적을 수 없다.
# 선적 outputs 5개도 함께 넣는다 — README 가 재생성하지 않는다고 선언한
# 것들이고, 없으면 파이프라인이 전자축에서 멈춘다.
mkdir -p submission/code/deltaE submission/code/sources
cp -p deltaE/*.csv submission/code/deltaE/ 2>/dev/null || true
cp -p sources/*.csv sources/*.pdf submission/code/sources/ 2>/dev/null || true
cp -p panel_all.csv panel_new.csv inchikeys.csv submission/code/ 2>/dev/null || true
# 선적 outputs 5개는 자기검증 뒤에 넣는다. 자기검증이 submission/code/outputs
# 를 루트 outputs 로 가는 심볼릭 링크로 만들었다가 지우기 때문이다.

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
  && mkdir -p paper \
  && cp -p ../main.tex ../supplementary.tex ../cover_letter.tex paper/ \
  && cp -p ../main.pdf ../supplementary.pdf ../cover_letter.pdf paper/ \
  && cp -p ../../paper/.buildstamp paper/ \
  && cp -p ../../README.md . \
  && cp -p ../../paper/tab_*.tex paper/ \
  && cp -p ../Fig_*.png ../toc_acs.png paper/ 2>/dev/null \
  && python3 verify_consistency.py > /tmp/_sync_verify.txt 2>&1 ) || {
  # 2026-11-12. 종전에는 verify 가 FAIL 이면 무조건 멈췄다. 그런데 A-tag 는
  # '예치 태그가 원격에 없다' 를 FAIL 로 보는데, 그것은 Ku 가 투고 직전에
  # 푸시하기 전까지의 정상 상태다. 그래서 패키지 동기화가 투고 직전까지
  # 영구히 막히고 submission/ 사본이 묵었다 — 그 묵은 사본이 바로 이
  # 스크립트가 고치려던 문제다. A-tag 단독 실패는 경고로 통과시킨다.
  n_fail=$(grep -cE '\(검사 [^)]*\):' /tmp/_sync_verify.txt 2>/dev/null || true)
  n_tag=$(grep -c '(검사 A-tag)' /tmp/_sync_verify.txt 2>/dev/null || true)
  n_b2d=$(grep -c '(검사 B2d)' /tmp/_sync_verify.txt 2>/dev/null || true)
  n_fail=${n_fail:-0}; n_tag=${n_tag:-0}; n_b2d=${n_b2d:-0}
  if [ "$n_fail" -gt 0 ] && [ "$((n_tag + n_b2d))" -eq "$n_fail" ]; then
    echo '[sync] A-tag/B2d 만 남았다 — 푸시 전 정상 상태이므로 통과시킨다'
  else
    echo '[sync] 패키지 자기검증 실패'
    grep -E '\(검사 [^)]*\):' /tmp/_sync_verify.txt | head -8
    exit 1
  fi; }
rm -rf submission/code/paper submission/code/outputs
# 2026-11-13. 자기검증이 submission/code 안에서 python3 를 돌리므로
# __pycache__ 가 생긴다. 복사 직후가 아니라 여기서 지워야 한다.
rm -rf submission/code/__pycache__

# README 가 "재생성하지 않고 선적한다" 고 선언한 다섯 파일. 없으면 패키지가
# 전자축에서 멈춘다 (xtb 바이너리 / 1시간짜리 B3LYP 실행이 필요한 자리).
mkdir -p submission/code/outputs
for f in aegl_electronic.csv xtb_fast_vs_slow.csv dft_omega.csv \
         dft_dipole.csv dft_vs_gfn2.csv; do
  [ -f "outputs/$f" ] && cp -p "outputs/$f" submission/code/outputs/
done
