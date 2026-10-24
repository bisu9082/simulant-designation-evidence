#!/usr/bin/env python3
"""서술자 공간 정의 정본 (2026-10-13 분리).

왜 분리하는가. 이 정의들은 지금까지 analysis_main.py 안에 있었고, 여섯 개
스크립트가 거기서 임포트했다. analysis_main.py 는 모듈이 아니라 스크립트여서
임포트하는 것만으로 argparse 가 돌고 전체 분석이 재실행되며
outputs/experiment_summary.json 을 덮어쓴다. 이 세션에서만 두 번, 어떤 스크립트를
임포트했다는 이유만으로 산출물이 구판으로 되돌아갔다. 정의만 여기 두고
analysis_main.py 는 이 모듈을 읽는 스크립트로 남긴다. 이 파일은 부작용이 없다.
"""

# 전자 서술자 전체. 표 보고용 집합이며, 구조공간에서 제외할 대상을 정의한다.
ELECTRONIC_28 = ['dipole_moment_D', 'E_LUMO_eV', 'E_HOMO_eV',
                 'electrophilicity_index_eV']

# Parr 친전자성 지수 omega = mu^2/2eta 는 HOMO/LUMO 의 결정론적 함수다
# (정의식 대비 최대 편차 3.6e-15 eV). 독립 좌표가 아니므로 어떤 거리공간에도
# 넣지 않는다. 상관필터(|r|>0.95)는 r(E_LUMO, omega) = -0.899 이므로 못 잡는다.
# 2026-10-13 GATE 8 3차 M3/C1 지적 반영.
OMEGA_DERIVED = 'electrophilicity_index_eV'
ELECTRONIC_DIST = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]

# 인 국소 서술자. 비인 화합물에서 정의되지 않으므로 전 패널 거리에 쓸 수 없다.
P_LOCAL = ['q_P', 'q_phosphoryl_O', 'wiberg_PO', 'fukui_f_plus_P']

META_NUM = ['n_conformers', '_n_conf_xtb']           # 실행 메타데이터, 화학 아님

PROTECTED = ['MR_crippen', 'logP_crippen', 'TPSA', 'HBA', 'lser_B_proxy',
             'lser_S_proxy'] + ELECTRONIC_28

R_CUT = 0.95


_DEP_CACHE = {}

DEP_TOL = 1e-10          # 정확 종속 판정 허용오차
DEP_MAX_TERMS = 2        # 몇 개 좌표의 조합까지 볼 것인가

# 종속관계가 발견되면 그 관계에 묶인 좌표 중 어느 것을 버릴지 정해야 한다.
# omega 선례를 그대로 따른다. 파생 합성량을 버리고 원시량을 남긴다. 그러지
# 않으면 MR_crippen(몰굴절) 같은 해석 가능한 원시 서술자가 빠지고 그것으로부터
# 만들어진 lser_S_proxy 가 남는다. 아래 열들을 먼저 검사해 먼저 떨어뜨린다.
DERIVED_FIRST = ('lser_A_proxy', 'lser_B_proxy', 'lser_S_proxy', 'lser_L_proxy',
                 'is_organophosphorus', 'n_P', 'n_PdO', 'heavy_atoms')


def exact_dependents(df, cols, tol=DEP_TOL, max_terms=DEP_MAX_TERMS):
    """다른 좌표 max_terms 개의 정확한 아핀함수인 좌표를 찾는다.

    2026-10-24 신설. 본 논문은 Parr omega 가 두 궤도에너지의 정확한 함수인데
    |r|>0.95 상관필터가 못 잡는다는 것을 방법론적 기여로 내세운다. 그 교훈을
    선언만 하고 omega 한 열에만 손으로 적용하고 있었다. 심사에서 구조공간의
    lser_S_proxy = (TPSA + MR_crippen)/100 이 같은 결함으로 남아 있다는 지적이
    나왔고, 28종 패널에서 실제로 유지되고 있었다(35종에서는 상관필터가 우연히
    걸러낸다). 규칙으로 만든다.

    상관필터가 못 잡는 것은 두 열 이상의 조합인 경우다. 한 열과의 정확한
    비례관계는 |r| = 1 이므로 상관필터가 먼저 떨어뜨린다. 그래서 이 검사는
    상관필터 뒤에 돌리며, 남는 것은 다중 좌표 조합뿐이다.

    omega 처럼 비선형인 종속은 이 선형 검사로 잡히지 않는다. 그쪽은 정의식으로
    따로 확인했고 ELECTRONIC_DIST 가 제외한다. 두 경로를 모두 유지한다.
    """
    import itertools
    import numpy as np
    n = len(df)
    # pick() 은 루프 안에서 여러 번 불린다. 이 검사는 열 수의 세제곱에 비례해
    # lstsq 를 돌려 한 번에 7초가 걸리므로, 캐시가 없으면 파이프라인이 몇십 분
    # 느려진다 (2026-10-24 측정).
    key = (tuple(cols), n,
           float(np.nansum(df[list(cols)].to_numpy(float))) if len(cols) else 0.0)
    if key in _DEP_CACHE:
        return _DEP_CACHE[key]
    keep, dropped = list(cols), []
    order = ([c for c in DERIVED_FIRST if c in cols]
             + [c for c in cols if c not in DERIVED_FIRST])
    for t in order:
        others = [c for c in keep if c != t]
        y = df[t].to_numpy(float)
        hit = None
        for k in range(1, max_terms + 1):
            if hit:
                break
            for combo in itertools.combinations(others, k):
                X = np.c_[df[list(combo)].to_numpy(float), np.ones(n)]
                coef, *_ = np.linalg.lstsq(X, y, rcond=None)
                if np.abs(X @ coef - y).max() < tol:
                    hit = combo
                    break
        if hit:
            keep.remove(t)
            dropped.append((t, tuple(hit)))
    _DEP_CACHE[key] = (keep, dropped)
    return keep, dropped


def pick(df, cols, protected):
    """분산 0 제거 → |r|>R_CUT 쌍에서 비보호 열 제거 → 정확 종속 좌표 제거."""
    cols = [c for c in cols if df[c].std(ddof=1) > 1e-12]
    C = df[cols].corr().abs().to_numpy()
    keep, dropped = list(cols), []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = cols[i], cols[j]
            if a not in keep or b not in keep or C[i, j] <= R_CUT:
                continue
            drop = b if (a in protected or b not in protected) else a
            if drop in protected and (a if drop == b else b) in protected:
                drop = b                                   # 둘 다 보호면 뒤를 버림
            keep.remove(drop); dropped.append((drop, a if drop == b else b,
                                               float(C[i, j])))
    keep, dep = exact_dependents(df, keep)
    dropped.extend(dep)
    return keep, dropped
