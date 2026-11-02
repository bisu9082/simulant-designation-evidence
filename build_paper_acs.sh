#!/bin/sh
# ACS Sensors 판본 빌드 (2026-11-02 신설).
#
# 본문과 SI 가 xr 로 서로를 가리킨다 (본문 -> S-, SI -> M-). 상호 의존이므로
# 한쪽만 돌리면 상대 라벨이 undefined 로 남는다. 두 문서를 번갈아 두 바퀴
# 돌려야 번호가 안정된다. 그 순서를 여기에 고정한다.
#
# 쓰는 법:  sh build_paper_acs.sh
set -e
cd "$(dirname "$0")/paper_acs"

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
echo '[6/6] main (최종)'
pdflatex -interaction=nonstopmode main > main_build.log 2>&1 || true

fail=0
for f in main supporting; do
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
for f in main supporting; do
    if [ -f "$f.pdf" ]; then
        printf '  %-12s %s쪽\n' "$f.pdf" \
            "$(pdfinfo "$f.pdf" | awk '/^Pages/{print $2}')"
    else
        echo "  $f.pdf 없음"; fail=1
    fi
done

[ "$fail" -eq 0 ] && echo '[빌드] 오류 없음' || { echo '[빌드] 위 항목 확인'; exit 1; }
