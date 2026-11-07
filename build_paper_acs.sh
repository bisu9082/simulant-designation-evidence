#!/bin/sh
# ACS Sensors 판본 빌드 (2026-11-02 신설).
#
# 본문과 SI 가 xr 로 서로를 가리킨다 (본문 -> S-, SI -> M-). 상호 의존이므로
# 한쪽만 돌리면 상대 라벨이 undefined 로 남는다. 두 문서를 번갈아 두 바퀴
# 돌려야 번호가 안정된다. 그 순서를 여기에 고정한다.
#
# 쓰는 법:  sh build_paper_acs.sh
set -e
cd "$(dirname "$0")"
ROOT=$(pwd)

# 도판 동기화 (2026-11-02 신설). 이 단계가 없었다. build_paper.sh 는 JHM
# 판본(paper/)에만 도판을 복사하므로, 파이프라인을 다시 돌려도 paper_acs/ 의
# png 는 구판으로 남았다. verify_consistency 의 ACS 도판 시차 검사가 이것을
# 잡아냈다. --figs 를 주면 생성기까지 다시 돌린다.
if [ "$1" = '--figs' ]; then
  # 2026-11-07 (11차). 이 목록에 둘이 빠져 있었다.
  #   make_fig_talk_extra.py  -> outputs/Fig_talk_tox.png = 원고의 Figure 3.
  #     아래 복사 목록에는 있는데 생성기가 없었으므로, 빈 저장소를 클론한
  #     사람은 Figure 3 을 만들 수 없고 경고 한 줄만 본다.
  #   make_fig_toxicity.py    -> outputs/fig_label_counts.json.
  #     검증기 REQUIRED_OUTPUTS 에 있어서, 없으면 사전 점검에서 멈춘다.
  #     파이프라인도 빌드도 부르지 않았다.
  for s in make_fig_main.py make_fig_h1.py make_fig_mde.py \
           make_fig_talk_response.py \
           make_fig_talk_extra.py make_fig_toxicity.py make_toc_acs.py; do
    # 2026-11-07. 종전 목록에 make_fig_plocal.py 가 있었는데 그런 파일은
    # 없다. Fig_plocal.png 은 make_fig_main.py 가 쓴다. [ -f ] 로 조용히
    # 건너뛰고 있었으므로, 목록이 틀렸다는 사실이 드러나지 않았다. 이제
    # 없는 파일은 중단 사유다.
    if [ ! -f "$s" ]; then
      echo "[fig] $s 없음 — 목록이 저장소와 어긋났다"; exit 1
    fi
    printf '[fig] %-28s ' "$s"
    if python3 "$s" >/dev/null 2>&1; then echo ok; else echo FAIL; exit 1; fi
  done
fi

for f in Fig_main.png Fig_plocal.png Fig_h1.png Fig_mde.png \
         Fig_talk_response.png Fig_talk_tox.png toc_acs.png toc_acs.pdf; do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "paper_acs/$f"
  else
    echo "[경고] outputs/$f 없음 — ACS 원고의 도판이 갱신되지 않는다" >&2
  fi
done


# 생성 표 동기화 (2026-11-03 신설). tab_*.tex 는 make_si_tables.py 가
# outputs/ 에 만든다. 손으로 유지했던 네 표가 반복적으로 구판 수치를
# 남긴 것이 전파 실패의 근원이었다. 복사를 빌드에 묶어 시차를 없앤다.
for f in tab_ablation.tex tab_family.tex tab_blockbalance.tex \
         tab_sparse.tex tab_plocal.tex tab_volatility.tex \
         tab_vpairing.tex tab_premisenull.tex tab_registers.tex \
         tab_designations.tex tab_dyepanel.tex; do
  if [ -f "outputs/$f" ]; then
    cp -p "outputs/$f" "paper_acs/$f"
  else
    echo "[경고] outputs/$f 없음 — make_si_tables.py 를 먼저 돌릴 것" >&2
    exit 1
  fi
done

cd paper_acs

pass() {
    pdflatex -interaction=nonstopmode "$1" >/dev/null 2>&1 || true
}

echo '[1/6] main (1차 — SI 라벨 없이)'
pass main
echo '[2/6] main 참고문헌'
bibtex main >/dev/null 2>&1 || true
echo '[3/6] supporting (1차 — main.aux 를 읽는다)'
pass supporting
bibtex supporting >/dev/null 2>&1 || true
pass supporting
echo '[4/6] main (2차 — supporting.aux 를 읽는다)'
pass main
pass main
echo '[5/6] supporting (2차 — 확정된 main.aux 를 읽는다)'
pass supporting
pass supporting
echo '[6/7] main (최종)'
pdflatex -interaction=nonstopmode main > main_build.log 2>&1 || true

# 2026-11-07 (11차). 커버레터를 빌드 목록에 넣는다.
#
# 왜. 이 스크립트는 본문과 SI 만 빌드했고 커버레터는 아예 부르지 않았다.
# 그래서 paper_acs/cover_letter.pdf 가 11-03 자로 멈춰 있었다 — 11차 교정
# 전체(제목 교체, 분량 사유 단락 신설)가 반영되지 않은 판본이다. 그대로
# 올리면 철회한 헤드라인 프레임의 커버레터를 투고하게 된다. 검증기는
# cover_letter.tex 를 읽으므로 글은 맞았고, 올라가는 PDF 만 틀렸다.
# 아래 시차 검사(검사 A)가 같은 자리를 다시 막는다.
echo '[7/7] cover letter'
pass cover_letter
pass cover_letter

fail=0
for f in main supporting cover_letter; do
    log="$f.log"
    if grep -q '^!' "$log" 2>/dev/null; then
        echo "  [$f] 치명 오류:"; grep -A2 '^!' "$log" | head -12; fail=1
    fi
    if grep -qi 'undefined' "$log" 2>/dev/null; then
        echo "  [$f] 끊어진 참조/인용:"
        grep -i 'undefined' "$log" | head -8; fail=1
    fi
    if grep -q 'ignored' "$log" 2>/dev/null; then
        echo "  [$f] achemso 가 절 제목을 삼켰다 (journal 코드 확인):"
        grep 'ignored' "$log" | head -4; fail=1
    fi
done

echo
echo '--- 산출 ---'
for f in main supporting cover_letter; do
    if [ -f "$f.pdf" ]; then
        printf '  %-12s %s쪽\n' "$f.pdf" \
            "$(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}')"
    else
        echo "  $f.pdf 없음"; fail=1
    fi
done

# 2026-11-07. 커버레터가 분량을 숫자로 적고 검증기가 그 숫자를
# outputs/length_audit.json 과 대조한다. 그 파일은 빌드된 PDF 를 읽어야
# 만들어지므로 여기서 갱신한다. 빌드가 성공한 뒤에만 돈다.
if [ "$fail" -eq 0 ]; then
    # 이 지점의 cwd 는 paper_acs/ 다. 상대경로에 의지하지 않는다.
    ( cd "$ROOT" && python3 wordcount.py >/dev/null 2>&1 ) \
        || echo '  [주의] wordcount.py 실패 — length_audit.json 갱신 안 됨'
fi

[ "$fail" -eq 0 ] && echo '[빌드] 오류 없음' || { echo '[빌드] 위 항목 확인'; exit 1; }
