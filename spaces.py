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


def pick(df, cols, protected):
    """분산 0 제거 → |r|>R_CUT 쌍에서 비보호 열을 떨어뜨린다."""
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
    return keep, dropped
