#!/usr/bin/env python3
"""친전자성 지수를 Parr 정의로 교정한다 (2026-10-12).

GATE 8 4차 지적 (리뷰어 B M3, 리뷰어 C 신규 발견, 내가 재확인).

본문 2.2절과 2.7절은 omega = mu^2 / 2*eta 를 Parr 등의 정의로 인용하고,
eta 를 "half the HOMO-LUMO gap" 이라고 명시한다. 그런데 두 코드 경로가
모두 eta 에 간극 '전체' 를 넣는다.

    xtb_fast.py:76              mu, eta = (ho+lu)/2.0, (lu-ho)
    xtb_descriptors.py:190      eta = E_LUMO - E_HOMO

사린으로 대조: HOMO -11.9449, LUMO -2.6193, mu -7.2821, gap 9.3256.
    저장값          2.8432  (= mu^2 / (2*gap))
    Parr 정의       5.6864  (= mu^2 / (2*(gap/2)) = mu^2 / gap)
즉 보고된 모든 omega 가 정확히 절반이다. 프로젝트 금지사항의
"실제로 사용하지 않은 방법을 사용했다고 서술" 에 해당한다.

영향 범위. omega 는 표준화 후 거리에 들어가므로 상수배는 거리를 바꾸지
않는다. 아래에서 그것을 수치로 확인한다. 바뀌는 것은 Table 3 의 omega 열,
Fig 1(f), 3.3절 본문 수치, 2.7절의 NaCN 서술(545 -> 1090 eV) 뿐이다.

이 스크립트는 코드 두 곳을 고치고, 캐시된 값을 HOMO/LUMO 에서 재계산한다
(xtb 재실행 불필요 — omega 는 궤도 에너지의 결정론적 함수다).
"""
import numpy as np
import pandas as pd


def omega_parr(homo, lumo):
    """Parr 정의. eta = (LUMO - HOMO)/2 이므로 omega = mu^2 / (LUMO - HOMO)."""
    mu = (homo + lumo) / 2.0
    gap = lumo - homo
    return np.where(gap > 0, mu * mu / gap, np.nan)


def patch_code():
    changed = []
    # xtb_fast.py
    p = open('xtb_fast.py').read()
    old = '            mu, eta = (ho + lu) / 2.0, (lu - ho)'
    new = ('            # Parr 정의: eta 는 간극의 절반이다 (2026-10-12 교정).\n'
           '            # 종전에는 간극 전체를 넣어 omega 가 전부 절반이었다.\n'
           '            mu, eta = (ho + lu) / 2.0, (lu - ho) / 2.0')
    if old in p:
        open('xtb_fast.py', 'w').write(p.replace(old, new))
        changed.append('xtb_fast.py')

    # /home/claude/xtb_descriptors.py (정밀 프로토콜)
    import os
    xd = '/home/claude/xtb_descriptors.py'
    if os.path.exists(xd):
        q = open(xd).read()
        o2 = "        eta = out['E_LUMO_eV'] - out['E_HOMO_eV']"
        n2 = ("        # Parr 정의: eta = 간극의 절반 (2026-10-12 교정)\n"
              "        eta = (out['E_LUMO_eV'] - out['E_HOMO_eV']) / 2.0")
        if o2 in q:
            q = q.replace(o2, n2)
            q = q.replace(
                "  electrophilicity_index_eV  ω = μ² / (2η),  μ=(E_H+E_L)/2,  η=E_L−E_H",
                "  electrophilicity_index_eV  ω = μ² / (2η),  μ=(E_H+E_L)/2,  η=(E_L−E_H)/2")
            # f+ 방향 표기 교정
            q = q.replace('인 원자 친전자 Fukui 지수 f(+)',
                          '인 원자 친핵성공격 감수성 Fukui 지수 f(+)')
            open(xd, 'w').write(q)
            changed.append(xd)
    return changed


def patch_cache(path, homo='E_HOMO_eV', lumo='E_LUMO_eV',
                col='electrophilicity_index_eV'):
    D = pd.read_csv(path)
    if homo not in D.columns:
        return None
    before = D[col].copy()
    D[col] = omega_parr(D[homo].to_numpy(float), D[lumo].to_numpy(float))
    ok = D[col].notna() & before.notna()
    ratio = (D.loc[ok, col] / before[ok]).replace([np.inf, -np.inf], np.nan).dropna()
    D.to_csv(path, index=False)
    return len(D), float(ratio.min()), float(ratio.max())


def distances_unchanged():
    """표준화 거리가 omega 의 상수배에 불변임을 수치로 확인한다."""
    from spaces import ELECTRONIC_28, META_NUM, PROTECTED, pick
    df = pd.read_csv('descriptor_cache_xtb.csv')
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == 28 and c not in META_NUM]
    out = {}
    for tag, cols in [('structural', [c for c in num if c not in ELECTRONIC_28]),
                      ('electronic', [c for c in num if c in ELECTRONIC_28]),
                      ('merged', num)]:
        keep, _ = pick(df, cols, PROTECTED)
        X = df[keep].to_numpy(float)
        Xz = (X - X.mean(0)) / X.std(0, ddof=1)
        out[tag] = np.sqrt(((Xz[:, None, :] - Xz[None, :, :]) ** 2).sum(-1))
    return out


if __name__ == '__main__':
    print('[1] 교정 전 거리 행렬 보존')
    before = distances_unchanged()

    print('[2] 코드 교정')
    for c in patch_code():
        print('   ', c)

    print('[3] 캐시 재계산 (HOMO/LUMO 에서 결정론적으로 산출)')
    for f in ['descriptor_cache_xtb.csv', 'outputs/aegl_electronic.csv',
              'outputs/xtb_fast_vs_slow.csv']:
        try:
            r = patch_cache(f)
            if r:
                print(f'    {f:<34} {r[0]}행  배율 {r[1]:.4f} ~ {r[2]:.4f}')
        except Exception as e:
            print(f'    {f:<34} 건너뜀 ({e})')
    # compare 파일의 _slow 열도 같은 규약으로
    try:
        M = pd.read_csv('outputs/xtb_fast_vs_slow.csv')
        if 'electrophilicity_index_eV_slow' in M.columns:
            C = pd.read_csv('descriptor_cache_xtb.csv')
            m = dict(zip(C.compound_name, C.electrophilicity_index_eV))
            M['electrophilicity_index_eV_slow'] = M.compound_name.map(m)
            M.to_csv('outputs/xtb_fast_vs_slow.csv', index=False)
            print('    xtb_fast_vs_slow.csv 의 _slow 열도 동기화')
    except Exception as e:
        print('    _slow 동기화 건너뜀', e)

    print('[4] 교정 후 거리 행렬 대조')
    after = distances_unchanged()
    allok = True
    for k in before:
        d = float(np.abs(before[k] - after[k]).max())
        # 허용오차 1e-4. omega 를 CSV 에 저장된 HOMO/LUMO(소수 5자리)에서
        # 재계산하므로, 원래의 전정밀도 omega 의 정확한 2배와 미세하게 다르다.
        # 순수한 반올림 잔차이며 거리 규모(약 1~15)에 비해 1e-6 수준이다.
        print(f'    {k:<11} 최대 차이 {d:.3e}  '
              f'{"불변(반올림 잔차 이하)" if d < 1e-4 else "*** 변함 ***"}')
        allok &= d < 1e-4
    print(f'\n[판정] 표준화 거리 {"불변 — p값·H1·후보 전부 영향 없음" if allok else "변함, 재확인 필요"}')
