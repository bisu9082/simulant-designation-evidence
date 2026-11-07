#!/usr/bin/env python3
"""지정 실패의 불변성 감사.

왜 이 스크립트가 있는가.

10·11·12·13차 네 차례 모두 헤드라인이 같은 이유로 무너졌다. 긍정형
주장을 만들려고 고른 통계량이 선언되지 않은 분석자 선택의 함수였다.
13차 에디터 게이트가 안 1(자유도 범위를 헤드라인으로)에도 같은 진단을
내렸다. 범위 0~6 과 6~30 은 감사집합에 무엇을 넣는가라는, 역시
선언되지 않은 선택의 함수다.

그런데 예치물에는 그 선택에 **불변인** 양성 명제가 이미 있었다.
감사된 모든 분석 구성에서 유효 매칭이 한 번도 성립하지 않는 작용제가
있고, 성립하지 않는 이유도 같다. 이 스크립트는 그것을 센다.

세는 방식.

감사 구성은 두 산출물에서 온다.
  outputs/match_rule_audit.json  패널×지정집합×필터규칙×비교군판독 18셀
  outputs/sparse_indicators.json 서술자 인코딩 7변형
두 파일의 작용제별 판정 스키마가 다르므로 아래에서 각각 읽는다.

판정 기준은 원 감사와 같다. match_rule_audit 는 within_threshold 와
beats_comparison_set 이 모두 참일 때 유효 매칭이고, sparse_indicators
는 match 플래그가 그 둘을 이미 합친 값이다.

보고하는 것.
  n_cells          작용제가 등장한 감사 셀 수
  n_match          그중 유효 매칭이 성립한 셀 수
  invariant_fail   n_match == 0 (모든 구성에서 실패)
  blockers         실패 셀에서 기록된 차단자 (더 가까운 화합물)
  d_norm_range     정규화 거리의 범위. 역치는 THRESHOLD.
  fail_mode        threshold_only / blocker_only / both

[C-3] 이 스크립트는 기존 예치 산출물만 읽는다. 새로 계산하는 값은
집계뿐이고 거리나 매칭을 다시 판정하지 않는다.
"""
import json
import os

THRESHOLD = 2.0          # match_rule_audit.py:28 과 같은 값
OUT = 'outputs/designation_failure_audit.json'


def load_cells():
    """[(셀이름, {작용제: 판정dict})] 를 돌려준다."""
    cells = []

    p = 'outputs/match_rule_audit.json'
    if os.path.exists(p):
        d = json.load(open(p))
        for blk, arms in d.items():
            if blk == 'panel_swap':      # 패널 교체는 구성이 아니라 별 감사
                continue
            for arm, v in arms.items():
                per = {}
                for a in v.get('agents', []):
                    ok = bool(a.get('within_threshold')) and \
                         bool(a.get('beats_comparison_set'))
                    per[a['agent']] = dict(
                        match=ok, d_norm=a.get('d_norm'),
                        nearest=a.get('nearest_simulant'),
                        blocker=a.get('blocker'),
                        within_threshold=a.get('within_threshold'))
                cells.append((f'match_rule/{blk}/{arm}', per))
    else:
        print(f'  [수집 불가: {p} 없음]')

    p = 'outputs/sparse_indicators.json'
    if os.path.exists(p):
        d = json.load(open(p))
        for v in d.get('variants', []):
            per = {}
            for nm, a in v.get('agents', {}).items():
                dn = a.get('d_norm')
                per[nm] = dict(
                    match=bool(a.get('match')), d_norm=dn,
                    nearest=a.get('nearest'), blocker=None,
                    within_threshold=(None if dn is None
                                      else bool(dn <= THRESHOLD)))
            cells.append((f"sparse/{v['variant']}", per))
    else:
        print(f'  [수집 불가: {p} 없음]')

    return cells


def main():
    cells = load_cells()
    if not cells:
        print('[수집 불가: 감사 셀이 없다]')
        return

    agents = {}
    for name, per in cells:
        for nm, a in per.items():
            agents.setdefault(nm, []).append((name, a))

    report = dict(generated_by='designation_failure_audit.py',
                  threshold=THRESHOLD,
                  n_cells=len(cells),
                  cell_names=[n for n, _ in cells],
                  agents={})

    print('=' * 78)
    print('[지정 실패 감사] 감사된 모든 분석 구성에서의 작용제별 판정')
    print('=' * 78)
    print(f'  감사 셀 {len(cells)}개')
    print()
    print(f'  {"작용제":<22}{"매칭/셀":>9}  {"d_norm":^15} 차단자')

    for nm in sorted(agents, key=lambda x: (
            sum(1 for _, a in agents[x] if a['match']), x)):
        rec = agents[nm]
        nmat = sum(1 for _, a in rec if a['match'])
        ds = [a['d_norm'] for _, a in rec if a['d_norm'] is not None]
        bl = sorted({a['blocker'] for _, a in rec
                     if not a['match'] and a.get('blocker')})
        # 실패 양태. 역치를 넘어서 실패한 셀과 차단자 때문에 실패한 셀.
        n_over = sum(1 for _, a in rec
                     if not a['match'] and a.get('within_threshold') is False)
        n_blk = sum(1 for _, a in rec
                    if not a['match'] and a.get('blocker'))
        nfail = len(rec) - nmat
        if nfail == 0:
            mode = 'never_fails'
        elif n_over == nfail:
            mode = 'threshold_only'
        elif n_blk == nfail and n_over == 0:
            mode = 'blocker_only'
        else:
            mode = 'both'

        report['agents'][nm] = dict(
            n_cells=len(rec), n_match=nmat,
            invariant_fail=bool(nmat == 0),
            invariant_match=bool(nmat == len(rec)),
            d_norm_min=(min(ds) if ds else None),
            d_norm_max=(max(ds) if ds else None),
            n_cells_over_threshold=n_over,
            n_cells_blocked=n_blk,
            fail_mode=mode,
            blockers=bl,
            nearest_designated=sorted({a['nearest'] for _, a in rec
                                       if a.get('nearest')}))
        rng = f'{min(ds):.3f}-{max(ds):.3f}' if ds else '[수집 불가]'
        mark = ' <-' if nmat == 0 else '   '
        print(f'  {nm:<22}{nmat:>4}/{len(rec):<4}{mark}{rng:^15} '
              f'{", ".join(bl) if bl else "(기록 없음)"}')

    inv = sorted(k for k, v in report['agents'].items()
                 if v['invariant_fail'])
    report['invariant_failures'] = inv
    report['n_agents'] = len(report['agents'])
    report['n_invariant_failures'] = len(inv)

    print()
    print(f'  감사된 {len(cells)}개 구성 전부에서 유효 매칭이 성립하지 '
          f'않는 작용제: {len(inv)}/{len(report["agents"])}종')
    for nm in inv:
        v = report['agents'][nm]
        print(f'    {nm:<22} {v["fail_mode"]:<15} '
              f'{", ".join(v["blockers"]) if v["blockers"] else ""}')

    with open(OUT, 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print(f'\n[done] {OUT}')


if __name__ == '__main__':
    main()
