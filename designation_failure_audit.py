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
                # sparse_indicators 는 차단자를 계산하지 않는다. 종전에
                # blocker=None 을 넣어서 '차단자 없음' 과 구별되지
                # 않았고, 그 결과 fail_mode 가 틀린 라벨을 냈다
                # (2026-11-07 에디터 게이트 2회차). 미기록을 명시한다.
                per[nm] = dict(
                    match=bool(a.get('match')), d_norm=dn,
                    nearest=a.get('nearest'),
                    blocker=None, blocker_recorded=False,
                    within_threshold=(None if dn is None
                                      else bool(dn <= THRESHOLD)))
            cells.append((f"sparse/{v['variant']}", per))
    else:
        print(f'  [수집 불가: {p} 없음]')

    return cells


def _threshold_free(rec):
    """역치를 뺀 판정. (통과, 실패, 미정) 을 셈한다.

    분모는 n_match 와 같은 전체 셀 수여야 한다. 통과 + 실패 + 미정 =
    len(rec) 를 항등식으로 지킨다. 2026-11-09 게이트 3회차 교정.
    """
    ok = und = 0
    for _, a in rec:
        if a['match']:
            ok += 1                       # 두 조건 다 통과했으므로 비교군도 통과
            continue
        if a.get('blocker_recorded', True):
            # match_rule 셀. beats_comparison_set 이 blocker 로 표현돼 있다.
            if a.get('blocker') is None:
                ok += 1                   # 비교군은 통과, 역치에서만 실패
        else:
            # sparse 셀. 역치가 구속이 아니었다면 비교군에서 실패한 것.
            if a.get('within_threshold') is False:
                und += 1                  # 역치 초과라 비교군 조건을 알 수 없다
    return dict(n_match_no_threshold=ok,
                n_fail_no_threshold=len(rec) - ok - und,
                n_undetermined_no_threshold=und)


def verdict_vectors(cells):
    """25구성을 '같은 판정을 내는 묶음' 으로 접는다.

    왜. 2026-11-07 에디터 게이트 2회차가 25구성이 서로 다른 분석
    선택이 아니라고 적발했다. sparse/baseline 은
    match_rule/panel28_incumbent3/deterministic/wide 와 같은 셀이고,
    panel28 과 panel35+lit8 에서 wide 와 narrow 는 비교군이 동일하며,
    sparse 7변형은 희소성 절단값 하나를 중첩 샘플링한 것이다.
    25 를 분모로 쓰면 중복 가중치가 대비를 만든다.
    """
    groups = {}
    for name, per in cells:
        key = tuple(sorted(nm for nm, a in per.items() if a['match']))
        groups.setdefault(key, []).append(name)
    return groups


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
        # 차단자가 기록되지 않는 셀(sparse 7변형)을 따로 센다. 종전에는
        # 이 셀들이 '차단자 없음' 과 섞여서 else 분기로 떨어졌고,
        # n_over=0 인 작용제에 'both' 가 붙었다. 2026-11-07 교정.
        n_unrec = sum(1 for _, a in rec
                      if not a['match'] and not a.get('blocker_recorded',
                                                      True))
        nfail = len(rec) - nmat
        n_judged = nfail - n_unrec
        if nfail == 0:
            mode = 'never_fails'
        elif n_over == nfail:
            mode = 'threshold_in_every_failing_cell'
        elif n_over == 0 and n_blk == n_judged and n_judged > 0:
            mode = ('blocker_only' if n_unrec == 0
                    else 'blocker_only_where_recorded')
        elif n_over > 0 and n_blk > 0:
            mode = 'mixed'
        else:
            mode = 'undetermined'

        report['agents'][nm] = dict(
            n_cells=len(rec), n_match=nmat,
            invariant_fail=bool(nmat == 0),
            invariant_match=bool(nmat == len(rec)),
            d_norm_min=(min(ds) if ds else None),
            d_norm_max=(max(ds) if ds else None),
            n_cells_over_threshold=n_over,
            n_cells_blocked=n_blk,
            n_cells_blocker_unrecorded=n_unrec,
            fail_mode=mode,
            # [자유도: 역치] 역치 2.0 을 제거하고 '비교군에 더 가까운
            # 것이 없다' 만 요구하면 어떻게 되는가. 본문 652-655 가
            # 이미 28종 패널에서 3 -> 5 로 보고하는 축인데 당초 감사에는
            # 들어 있지 않았다. 2026-11-07 게이트 2회차 지적.
            #
            # [교정 게이트 3회차] 첫 구현이 틀렸다. 차단자가 기록된
            # 18셀만 분모로 삼아 25셀 분모의 n_match 와 나란히 냈고,
            # 그래서 11종 중 6종에서 단조성이 깨졌다 (PS 가 1 -> 0).
            # 제약을 제거하면 매칭은 줄 수 없다. 셀마다 비교군 조건만
            # 따로 판정한다.
            #   match_rule 셀: beats_comparison_set 이 그 조건이다.
            #   sparse 셀: match 는 두 조건의 연접이므로
            #     match=True            -> 비교군 조건 통과
            #     match=False, 역치 이내 -> 비교군 조건에서 실패
            #     match=False, 역치 초과 -> 비교군 조건 미정
            **_threshold_free(rec),
            blockers=bl,
            nearest_designated=sorted({a['nearest'] for _, a in rec
                                       if a.get('nearest')}))
        rng = f'{min(ds):.3f}-{max(ds):.3f}' if ds else '[수집 불가]'
        mark = ' <-' if nmat == 0 else '   '
        print(f'  {nm:<22}{nmat:>4}/{len(rec):<4}{mark}{rng:^15} '
              f'{", ".join(bl) if bl else "(기록 없음)"}')

    inv = sorted(k for k, v in report['agents'].items()
                 if v['invariant_fail'])
    # 역치를 자유도로 넣으면 남는 것. 역치 제거만으로 매칭이 되는
    # 작용제는 '모든 구성에서 실패' 가 아니다.
    inv_th = sorted(k for k in inv
                    if report['agents'][k]['n_match_no_threshold'] == 0)
    report['invariant_failures'] = inv
    report['invariant_failures_threshold_in_audit'] = inv_th
    report['n_agents'] = len(report['agents'])
    report['n_invariant_failures'] = len(inv)
    report['n_invariant_failures_threshold_in_audit'] = len(inv_th)

    groups = verdict_vectors(cells)
    report['n_distinct_verdict_vectors'] = len(groups)
    report['verdict_vectors'] = [
        dict(matched=list(k), n_cells=len(v), cells=sorted(v))
        for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1]))]

    print()
    print(f'  감사 셀 {len(cells)}개는 서로 다른 판정벡터 '
          f'{len(groups)}개로 접힌다. 25 는 독립 선택의 수가 아니다.')
    for g in report['verdict_vectors']:
        print(f'    셀 {g["n_cells"]:>2}개 -> 매칭 '
              f'{{{", ".join(g["matched"]) if g["matched"] else "없음"}}}')

    print()
    print(f'  역치 2.0 을 고정한 채로: 모든 구성에서 실패 '
          f'{len(inv)}/{len(report["agents"])}종')
    for nm in inv:
        v = report['agents'][nm]
        extra = ('' if v['n_match_no_threshold'] == 0 else
                 f'  <- 역치 제거하면 {v["n_match_no_threshold"]}셀에서 매칭')
        print(f'    {nm:<22} {v["fail_mode"]:<32} '
              f'{", ".join(v["blockers"]) if v["blockers"] else ""}{extra}')
    print()
    print(f'  역치를 자유도로 넣으면: {len(inv_th)}종 '
          f'({", ".join(inv_th)})')

    with open(OUT, 'w') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    print(f'\n[done] {OUT}')


if __name__ == '__main__':
    main()
