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

# 2026-10-27. conformer_w_max 가 빠져 있었다. 2026-10-25 에 앙상블 평균을
# 도입하면서 보고용으로 남긴 최대 Boltzmann 가중치인데, META_NUM 을 갱신하지
# 않아 모든 거리에 좌표로 들어갔다. 본문은 '실행 메타데이터 제외'라고
# 선언하고 있었다. 화학량이 아니고 임베딩 시드의 함수이기도 하다.
META_NUM = ['n_conformers', '_n_conf_xtb', 'conformer_w_max']  # 실행 메타데이터

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


def _pick_legacy(df, cols, protected):
    """2026-11-02 이전 상관필터. 입력 열 순서의 함수다. 비교 목적으로만 남긴다.

    쌍을 (i, j) 이중 루프로 훑으며 같은 상관군에서 '뒤쪽' 열을 버리므로,
    어느 좌표가 살아남는지가 CSV 의 열 배치로 결정된다. GATE 8 7차 R4 지적.
    order_invariance.py 가 이 함수와 pick() 의 분포를 나란히 보고한다.
    """
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


def pick(df, cols, protected):
    """분산 0 제거 → |r|>R_CUT 상관필터 → 정확 종속 좌표 제거.

    상관필터는 2026-11-02 에 교체했다. 이전 판은 쌍을 입력 열 순서로 훑으며
    뒤쪽을 버렸으므로 유지되는 좌표 집합이 CSV 의 열 배치의 함수였고, 그
    배치를 바꾸면 28종 패널의 유효매칭이 4건과 3건 사이를 오갔다. 결과가
    분석자가 고른 적 없는 변수의 함수이면 결과가 아니다.

    교체한 규칙은 이 논문이 만든 것이 아니라 Kuhn 의 caret findCorrelation
    규칙이다. |r| 이 가장 큰 쌍을 먼저 보고, 그 두 좌표 중 나머지 유지 좌표에
    대한 평균 |r| 이 큰 쪽(더 중복인 쪽)을 버린다. 하나를 버리면 상관행렬을
    다시 계산하고 역치 초과 쌍이 없을 때까지 반복한다.

    결정성을 위해 세 가지를 고정한다.
      (1) 입력 열을 이름순으로 정렬해 시작한다. 호출자가 어떤 순서로 넘겨도
          같은 집합에서 출발한다.
      (2) 보호 좌표(PROTECTED)가 한쪽에만 있으면 비보호 쪽을 버린다. 평균 |r|
          비교보다 앞선다. 본문이 해석하는 좌표를 중복 때문에 잃지 않기 위한
          것이고, 이전 판에도 있던 규칙이다.
      (3) 평균 |r| 차가 1e-12 안이면 이름이 사전순 뒤인 쪽을 버린다. 완전
          중복 열(|r| = 1)에서만 걸리는 경로다.

    (1) 때문에 반환되는 keep 은 항상 이름순이며, 뒤에 붙는 exact_dependents
    의 검사 순서도 그래서 결정적이 된다.
    """
    import numpy as np
    cols = sorted(c for c in cols if df[c].std(ddof=1) > 1e-12)
    keep, dropped = list(cols), []
    while len(keep) > 1:
        C = df[keep].corr().abs()
        A = C.to_numpy(float).copy()
        np.fill_diagonal(A, np.nan)
        iu = np.triu_indices(len(keep), 1)
        over = [(float(A[i, j]), keep[i], keep[j])
                for i, j in zip(*iu) if A[i, j] > R_CUT]
        if not over:
            break
        # 가장 강한 쌍. r 동률이면 이름으로 깨서 결정적으로 만든다.
        r, a, b = max(over, key=lambda t: (t[0], t[1], t[2]))
        pa, pb = a in protected, b in protected
        if pa != pb:
            drop = b if pa else a                      # 비보호 쪽을 버린다
        else:
            ma = float(np.nanmean(A[keep.index(a)]))
            mb = float(np.nanmean(A[keep.index(b)]))
            if abs(ma - mb) > 1e-12:
                drop = a if ma > mb else b             # 더 중복인 쪽을 버린다
            else:
                drop = max(a, b)                       # 이름 사전순 뒤
        keep.remove(drop)
        dropped.append((drop, b if drop == a else a, r))
    keep, dep = exact_dependents(df, keep)
    dropped.extend(dep)
    return keep, dropped
