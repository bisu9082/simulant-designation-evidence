#!/bin/sh
# Analytical Chemistry (achemso) 판본 빌드.
#
# 본문과 Supplementary 가 xr 로 서로를 가리킨다. 상호 의존이므로 한쪽만
# 돌리면 상대 라벨이 undefined 로 남는다. 두 문서를 번갈아 돌려야 번호가
# 안정된다. 그 순서를 여기에 고정한다.
#
# 2026-11-08 (라인 전환). 이 스크립트를 다시 썼다. 종전 판은 pdflatex
# 출력을 전부 >/dev/null 로 버리고 `|| true` 로 실패를 삼켰다. 그래서
# 빌드가 깨져도 쪽수 한 줄만 찍고 성공처럼 끝났고, 한 번은 그 때문에
# 구판 length_audit.json 을 읽고 분량을 잘못 보고했다. ACS 판본 빌드
# 스크립트에 넣어 둔 오류 검사(치명 오류·끊어진 참조·PDF 부재·빌드
# 지문)를 이쪽으로 가져온다.
#
# 쓰는 법:  sh build_paper.sh [--figs]
set -e
cd "$(dirname "$0")"
ROOT=$(pwd)

# --figs 를 주면 도판 생성기까지 다시 돌린다.
# 목록에 없는 파일은 중단 사유다. [ -f ] 로 조용히 건너뛰면 목록이
# 틀렸다는 사실이 드러나지 않는다 (ACS 쪽에서 실제로 그렇게 새어
# 나갔다 — make_fig_plocal.py 라는 없는 파일이 목록에 있었다).
if [ "$1" = '--figs' ]; then
  # Fig_main.png / Fig_spaces.png / Fig_plocal.png 은 셋 다
  # make_fig_main.py 가 쓴다. 생성기 수가 도판 수보다 적은 것이 정상이다.
  for s in make_fig_main.py make_fig_h1.py make_fig_mde.py \
           make_fig_toxicity.py; do
    if [ ! -f "$s" ]; then
      echo "[fig] $s 없음 — 목록이 저장소와 어긋났다"; exit 1
    fi
    printf '[fig] %-26s ' "$s"
    if python3 "$s" >"$ROOT/outputs/_log_$s.txt" 2>&1; then
      echo ok
    else
      echo FAIL
      tail -20 "$ROOT/outputs/_log_$s.txt"
      exit 1
    fi
  done
  # Graphical abstract. JHM 은 TOC 도판이 아니라 GA 를 요구한다.
  if [ -f make_toc_acs.py ]; then
    printf '[fig] %-26s ' make_toc_acs.py
    if python3 make_toc_acs.py >"$ROOT/outputs/_log_make_toc_acs.txt" 2>&1; then
      echo ok
    else
      echo FAIL; tail -20 "$ROOT/outputs/_log_make_toc_acs.txt"; exit 1
    fi
  fi
fi

# toc_acs.png 은 ACS 의 TOC 그래픽이다 (3.25 x 1.75 in 실척).
# 종전 JHM 빌드는 TOC 를 paper/ 로 복사하지 않아 본문이 참조할 수 없었다.
for f in Fig_main.png Fig_spaces.png Fig_plocal.png Fig_h1.png \
         Fig_mde.png Fig_toxicity.png toc_acs.png; do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "paper/$f"
  else
    echo "[중단] outputs/$f 없음 — sh build_paper.sh --figs 를 먼저 돌릴 것" >&2
    exit 1
  fi
done

# 생성 표 동기화. tab_*.tex 는 make_si_tables.py 가 outputs/ 에 만든다.
# 손으로 유지했던 네 표가 반복적으로 구판 수치를 남긴 것이 전파 실패의
# 근원이었다. 복사를 빌드에 묶어 시차를 없앤다.
for f in tab_ablation.tex tab_family.tex tab_blockbalance.tex \
         tab_sparse.tex tab_plocal.tex tab_volatility.tex \
         tab_vpairing.tex tab_premisenull.tex tab_registers.tex \
         tab_designations.tex tab_dyepanel.tex tab_failaudit.tex; do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "paper/$f"
  else
    echo "[중단] outputs/$f 없음 — make_si_tables.py 를 먼저 돌릴 것" >&2
    exit 1
  fi
done

cd paper

# 로그를 남긴다. 버리지 않는다.
pass() {
    pdflatex -interaction=nonstopmode "$1" >"$1.buildout" 2>&1 || true
}

echo '[1/6] main (1차 — SI 라벨 없이)'
pass main
echo '[2/6] main 참고문헌'
bibtex main >main.bibout 2>&1 || true
echo '[3/6] supplementary (1차 — main.aux 를 읽는다)'
pass supplementary
bibtex supplementary >supplementary.bibout 2>&1 || true
pass supplementary
echo '[4/6] main (2차 — supplementary.aux 를 읽는다)'
pass main
pass main
echo '[5/6] supplementary (2차 — 확정된 main.aux 를 읽는다)'
pass supplementary
pass supplementary
echo '[6/6] cover letter'
pass cover_letter
pass cover_letter

# 빌드 지문. 검사 A-pdf 가 이 파일과 현재 .tex 의 해시를 비교한다.
# mtime 으로 비교했더니 주입시험이 원고를 스냅샷으로 되돌릴 때 내용은
# 같고 mtime 만 새로워져서 거짓 실패가 났다. 내용으로 본다.
for f in main supplementary cover_letter; do
    sha256sum "$f.tex"
done > .buildstamp

fail=0
for f in main supplementary cover_letter; do
    log="$f.log"
    if grep -q '^!' "$log" 2>/dev/null; then
        echo "  [$f] 치명 오류:"; grep -A2 '^!' "$log" | head -12; fail=1
    fi
    if grep -qi 'undefined' "$log" 2>/dev/null; then
        echo "  [$f] 끊어진 참조/인용:"
        grep -i 'undefined' "$log" | head -8; fail=1
    fi
done

echo
echo '--- 산출 ---'
for f in main supplementary cover_letter; do
    if [ -f "$f.pdf" ]; then
        printf '  %-16s %s쪽  Overfull %s\n' "$f.pdf" \
            "$(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}')" \
            "$(grep -c 'Overfull' "$f.log" || true)"
    else
        echo "  $f.pdf 없음"; fail=1
    fi
done

# 커버레터가 분량을 숫자로 적고 검증기가 그 숫자를
# outputs/length_audit.json 과 대조한다. 그 파일은 빌드된 PDF 를 읽어야
# 만들어지므로 여기서 갱신한다. 빌드가 성공한 뒤에만 돈다.
if [ "$fail" -eq 0 ]; then
    ( cd "$ROOT" && python3 wordcount.py >outputs/_log_wordcount.txt 2>&1 ) \
        || echo '  [주의] wordcount.py 실패 — length_audit.json 갱신 안 됨'
fi

[ "$fail" -eq 0 ] && echo '[빌드] 오류 없음' || { echo '[빌드] 위 항목 확인'; exit 1; }
