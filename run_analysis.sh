#!/bin/sh
# 분석 파이프라인 정본 (2026-10-13 신설).
#
# 왜 필요한가. 이 프로젝트의 재발 결함은 계산 오류가 아니라 전파 실패였다.
# 서술자 정의를 바꾸면 그것을 읽는 하위 산출물을 전부 다시 만들어야 하는데,
# 순서가 머릿속에만 있었으므로 매번 몇 단계가 빠졌다. GATE 8 3차에서 지적된
# 도판-본문 불일치, 검증기가 들고 있던 낡은 p 값, outputs/ 와 paper/ 의 도판
# 시차가 모두 이 경로다. 순서를 파일로 고정한다.
#
# 2026-10-15 보강. 저장소를 빈 outputs/ 로 클론해 돌려 보니 두 가지가 빠져
# 있었다. match_by_panel.json 은 어떤 스크립트도 만들지 않는 고아 파일이었고
# (표 3·도판 2a·검증 사실 12건의 출처), 독성축 체인 전체가 이 목록에 없어
# 도판 생성과 검증이 FileNotFoundError 로 멈췄다. 둘 다 아래에 넣었다.
#
# 의존 순서 (패널·순열 축):
#   analysis_main        정본 공간 정의 -> experiment_summary.json, distance_table
#   perm_sidedness       단측/양측 p (verify_consistency 의 사실표가 읽는다)
#   sanity_checks        순열 귀무분포 .npy (도판 (d) 가 읽는다)
#   analyze_electronic   전자공간 단독 분석
#   normalize_nn         정규화 최근접 거리표
#   perm_statistics      패널 2종 x 통계량 6종
#   deltaE_premise       측정 배열 전제 검정
#   deltaE_power         Fisher 통합 + 검정력
#   rollback_c           방향 C 전용 (omega 종속, 차원 통제, Mantel, 중심화 바닥)
#   match_by_panel       패널 3구성 x 작용제 11종 최근접 모사체와 유효매칭
#   mde                  최소검출효과 — 음성 결과가 검정력 부족이 아님을 보인다
#   blockbalance         구조38 대 전자3 블록 가중 민감도 (SI res:blockbalance)
#   stat_supplement(2)   보정 체계 5종, 백색화 통제, WY minP, 서술자 부트스트랩
#
# 의존 순서 (독성축, Supplementary):
#   parse_aegl           EPA 원문(sources/) 파싱 -> toxicity_aegl.csv
#   build_tox_trainset   훈련집합 조립
#   validate_identity    화학적 동일성 검사 -> identity_validation.json.
#                        cwa_registry 가 이 파일에서 제외 레코드를 읽으므로
#                        qsar_tox 보다 반드시 먼저 와야 한다. FAIL 2건은
#                        원고가 보고한 제외 레코드이므로 비정상 종료를 허용한다.
#   qsar_tox             랜덤포레스트 + 적용영역 -> qsar_tox_summary.json, qsar_oof
#   h1_decoupling        충실도-안전 탈동조 -> h1_summary.json, h1_candidates
#   h1_screen_sensitivity 타당성심사 민감도 -> h1_screen_sensitivity.json
#   minority_regression  꼬리 회귀 진단 (qsar_oof 를 읽는다)
#   validate_aegl_parse  이중 판독 대조
#
# xtb_fast 는 이 목록에 없다. AEGL 230종에 GFN2-xTB 를 다시 돌리는 단계이고
# xtb 바이너리가 필요하다. 산출물 aegl_electronic.csv 와 xtb_fast_vs_slow.csv
# 는 저장소에 포함해 두었으므로 하위 단계가 그대로 돌아간다.
#
# 사용: sh run_analysis.sh && sh build_paper.sh --figs && python3 verify_consistency.py
set -e
cd "$(dirname "$0")"

run_strict() {
  for s in "$@"; do
    printf '[run] %-22s ' "$s"
    if python3 "$s" > "outputs/_log_${s%.py}.txt" 2>&1; then
      echo ok
    else
      echo FAIL
      tail -20 "outputs/_log_${s%.py}.txt"
      exit 1
    fi
  done
}

# 패널·순열 축
run_strict analysis_main.py perm_sidedness.py sanity_checks.py \
           analyze_electronic.py normalize_nn.py perm_statistics.py \
           deltaE_premise.py deltaE_power.py rollback_c.py match_by_panel.py \
           mde.py blockbalance.py pesticide_panel.py pesticide_family.py stat_supplement.py stat_supplement2.py nearest_agent_boot.py null_match_expectation.py

# 독성축. parse -> trainset -> 동일성검사 -> QSAR 순서를 지켜야 한다.
run_strict parse_aegl.py build_tox_trainset.py cite_facts.py

# 동일성 검사는 보고서다. FAIL 2건(Pentaborane 구조 불일치, 다성분 레코드)은
# 원고가 제외했다고 적은 그 두 건이므로 비정상 종료로 파이프라인을 세우지
# 않는다. 다만 cwa_registry 가 이 산출물을 읽으므로 건너뛸 수는 없다.
printf '[run] %-22s ' validate_identity.py
if python3 validate_identity.py > outputs/_log_validate_identity.txt 2>&1; then
  echo 'ok (FAIL 0건)'
else
  echo "ok ($(grep -c '^   FAIL' outputs/_log_validate_identity.txt)건 보고 — 원고가 제외한 레코드)"
fi
[ -f outputs/identity_validation.json ] || { echo '  identity_validation.json 생성 실패'; exit 1; }

run_strict qsar_tox.py h1_decoupling.py h1_screen_sensitivity.py \
           minority_regression.py validate_aegl_parse.py

# 다중검정 보정은 rollback_c 의 결과에서 파생된다.
python3 - <<'PY'
import json
import numpy as np
r = json.load(open('outputs/rollback_c.json'))
S = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
     'n_closer', 'n_valid']
rows = []
for lab, sp in r['panels'].items():
    ps = np.array([sp[s]['statistics'][k]['p'] for s in sp for k in S])
    o = np.sort(ps)
    m = len(o)
    bh = np.minimum.accumulate((o * m / np.arange(1, m + 1))[::-1])[::-1]
    rows.append((lab, m, float(o[0]), float(min(1, o[0] * m)), float(bh[0]),
                 int((bh < 0.05).sum())))
    print(f'  [BH] {lab:<38} 최소p {o[0]:.4f}  Bonferroni '
          f'{min(1, o[0]*m):.3f}  최소q {bh[0]:.3f}  q<0.05 {int((bh<0.05).sum())}건')
json.dump(rows, open('outputs/multiplicity_c.json', 'w'), indent=1)
PY

echo '[run] 완료. 다음: sh build_paper.sh --figs && python3 verify_consistency.py'
