#!/usr/bin/env python3
"""물리적 타당성 심사(4규칙)에 대한 H1 결론의 민감도 (2026-10-23 신설).

왜 필요한가. S0.3 절은 걸러낸 레코드를 되살려도 결론이 바뀌지 않는다고
서술한다. 종전에는 그 수치가 스크립트 밖에서 한 번 계산된 값이었고, 충실도
공간에서 omega 를 뺀 뒤로는 더 이상 맞지 않았다. 재현 가능한 산출물로 만든다.

방법. h1_decoupling.py 를 H1_KEEP_UNPHYSICAL=1 로 한 번 더 돌린다. 그 경로는
네 규칙으로 걸러낼 레코드를 풀에 그대로 두고 같은 시드·같은 복제수로 같은
검정을 수행하며, 정식 산출물은 건드리지 않는다.

출력: outputs/h1_screen_sensitivity.json
"""
import json
import os
import subprocess
import sys

RAW = 'outputs/h1_unphysical_kept_raw.json'


def main():
    base = json.load(open('outputs/h1_summary.json'))
    env = dict(os.environ, H1_KEEP_UNPHYSICAL='1')
    if os.path.exists(RAW):
        os.remove(RAW)
    r = subprocess.run([sys.executable, 'h1_decoupling.py'], env=env,
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(RAW):
        sys.stderr.write(r.stdout[-2000:] + r.stderr[-2000:])
        raise SystemExit('[sens] 민감도 경로 실패')
    kept = json.load(open(RAW))
    os.remove(RAW)

    n_screened = sum(base['pool_screen'].values())
    out = {'screen_rules': base['pool_screen'],
           'n_records_screened': n_screened,
           'screened': {'n_candidates': base['n_candidates'],
                        'spearman_r': round(base['spearman_r'], 4),
                        'spearman_ci': [round(v, 3)
                                        for v in base['spearman_ci']],
                        'verdict': base['verdict']},
           'kept': {'n_candidates': kept['n_candidates'],
                    'spearman_r': round(kept['spearman_r'], 4),
                    'spearman_ci': [round(v, 3) for v in kept['spearman_ci']],
                    'verdict': kept['verdict']},
           'delta_spearman': round(kept['spearman_r'] - base['spearman_r'], 4),
           'verdict_unchanged': kept['verdict'] == base['verdict']}
    with open('outputs/h1_screen_sensitivity.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(f"[sens] 심사 적용 {out['screened']['n_candidates']}종 "
          f"Spearman {out['screened']['spearman_r']:+.4f}  ->  "
          f"미적용 {out['kept']['n_candidates']}종 "
          f"{out['kept']['spearman_r']:+.4f}  "
          f"(차 {out['delta_spearman']:+.4f}, 판정 "
          f"{'불변' if out['verdict_unchanged'] else '변동'})")
    print('[sens] outputs/h1_screen_sensitivity.json')


if __name__ == '__main__':
    main()
