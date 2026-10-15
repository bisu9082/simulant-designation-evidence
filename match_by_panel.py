#!/usr/bin/env python3
"""outputs/match_by_panel.json 생성기 (2026-10-15 신설).

왜 필요한가. 이 파일은 표 3, 도판 2(a), verify_consistency 의 사실 12건
(패널별 유효매칭 수, DFP·CEES·DCNP 변위 거리, VX 분리거리)의 출처인데
어떤 스크립트도 만들지 않는 고아 파일이었다. permutation_sidedness.json 과
같은 경로이고, 저장소를 클론해 run_analysis.sh 를 돌리면 도판 생성과 검증이
FileNotFoundError 로 멈추는 것으로 드러났다. 생성기를 둬서 그 경로를 막는다.

정의는 기존 코드에서 그대로 가져온다. 공간과 거리는 rollback_c.zdist
(분산0 제거 -> |r|>0.95 상관필터 -> 표준화 -> 유클리드), 매칭 규칙은
normalize_nn 의 논리곱 (정규화거리 <= 2 AND 모든 비모사체보다 가까움),
정규화 분모는 패널 전체 최근접거리의 중앙값이다.
"""
import json

import numpy as np
import pandas as pd

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick

CASES = [
    ('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
    ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
    ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
     ['simulant', 'literature_simulant']),
]
THRESHOLD = 2.0


def merged_space(df):
    """본문 2.3절의 병합공간. omega 는 두 궤도에너지의 정확한 함수이므로 제외."""
    n = len(df)
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == n and c not in META_NUM]
    return [c for c in num if c != OMEGA_DERIVED]


def zdist(df, cols):
    keep, _ = pick(df, cols, PROTECTED)
    X = df[keep].to_numpy(float)
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    return D, len(keep)


def one_panel(path, sim_cls):
    df = pd.read_csv(path)
    D, nk = zdist(df, merged_space(df))
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))

    kl = df.compound_class.tolist()
    name = df.compound_name.tolist()
    i_agent = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_cls]
    # 비모사체 = 작용제도 모사체도 아닌 나머지. 매칭 2조건 중 둘째의 비교군.
    i_other = [i for i, k in enumerate(kl)
               if k != 'real_agent' and k not in sim_cls]

    agents = []
    for i in i_agent:
        js = min(i_sim, key=lambda j: D[i, j])
        d_sim = D[i, js]
        d_oth = min(D[i, j] for j in i_other)
        d_norm = d_sim / nn_med
        agents.append({
            'agent': name[i],
            'nearest_simulant': name[js],
            'd_norm': round(float(d_norm), 3),
            'match': bool(d_norm <= THRESHOLD and d_sim < d_oth),
        })
    return {
        'n_dims': nk,
        'nn_median': round(nn_med, 3),
        'n_match': int(sum(a['match'] for a in agents)),
        'agents': agents,
    }


def main():
    out = {label: one_panel(path, cls) for label, path, cls in CASES}
    with open('outputs/match_by_panel.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    for label, v in out.items():
        print(f'  {label:<30} {v["n_dims"]}차원  '
              f'전형최근접 {v["nn_median"]:.3f}  유효매칭 {v["n_match"]}/11')
    return out


if __name__ == '__main__':
    main()
