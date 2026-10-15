#!/usr/bin/env python3
"""
H1 결합해제 검정 — 충실도 축과 안전성 축이 직교하는가?

후보 풀: EPA AEGL Final 중 구조·전자서술자·실측 AEGL-3 를 모두 가진 화합물.
  충실도 = 전자공간(HOMO, LUMO, omega, dipole, z표준화)에서 진본 작용제까지 최단거리
  안전성 = 실측 log10(AEGL-3 10분, mg/m3)  (클수록 안전)

H1 기각(강한 음의 상관: 가까울수록 독하다) -> 충실+안전은 구조적으로 불가 -> RQ-3 폐기
H1 유지(무상관)                              -> '가깝고 안전한' 사분면 존재 -> RQ-3 진행

이중용도 가드: 결과는 '어떤 후보를 배제할지'로만 쓴다. 구조 경보(OP 항콜린에스터라제
파마코포어)에 걸리는 후보는 상관과 무관하게 전량 제외한다 (걸러내기 전용, 서열화 아님).
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

SEED, N_BOOT = 42, 10000
rng = np.random.default_rng(SEED)
EL = ['E_HOMO_eV', 'E_LUMO_eV', 'electrophilicity_index_eV', 'dipole_moment_D']

# ---------------------------------------------- 구조 경보 (걸러내기 전용)
# Tier A = 작용제 모티프 (무조건 배제) / Tier B = 주의 / Tier C = 기록만
ALERTS = {
    # Tier A — P(=O) 에 '좋은 이탈기'가 붙은 것만이 신경작용제 모티프다.
    #          F / CN / 티오알킬 / 활성화 아릴옥시. 단순 알콕시(DMMP)는 여기 아님.
    'A_OP_leaving_F_CN': '[PX4](=[O,S])(~[#6,#7,#8])(~[#6,#7,#8])[F,$(C#N)]',
    'A_OP_thiolate': '[PX4](=[O,S])(~[#6,#7,#8])(~[#6,#7,#8])[SX2][#6]',
    'A_OP_nitroaryloxy': '[PX4](=[O,S])(~[#6,#7,#8])(~[#6,#7,#8])Oc1ccc([N+](=O)[O-])cc1',
    'A_arsenical': '[As]',
    'A_sulfur_mustard': '[SX2]([CH2][CH2][Cl,Br,I])[CH2][CH2][Cl,Br,I]',
    'A_nitrogen_mustard': '[NX3]([CH2][CH2][Cl,Br,I])[CH2][CH2][Cl,Br,I]',
    # Tier A — 혈액작용제(시안) 계열: 흡입 시 HCN 방출
    'A_cyanide_inorganic': '[#6-]#[#7]',
    'A_cyanogen_halide': '[F,Cl,Br,I]C#N',
    'A_hydrogen_cyanide': '[CX2H1]#[NX1]',
    # Tier B — 역사적 최루·수포 모티프 (알파-할로 카보닐, 할로겐화 옥심)
    'B_alpha_halo_carbonyl': '[CX3](=O)[CX4][F,Cl,Br,I]',
    'B_halo_oxime': '[CX3](=[NX2][OX2H])[F,Cl,Br,I]',
    # Tier B — 질식작용제 계보(포스겐 유도체: 클로로포메이트·아실할라이드)
    'B_chloroformate': '[OX2]([#6])[CX3](=O)[F,Cl,Br]',
    'B_acyl_halide': '[CX3](=O)[F,Cl,Br]',
    'B_carbonyl_dihalide': '[CX3](=O)([F,Cl,Br])[F,Cl,Br]',
    # Tier B — 실용 배제: 격렬 가수분해 / 강산화제
    'B_halosilane': '[Si][F,Cl,Br,I]',
    'B_sulfonyl_halide': '[SX4](=O)(=O)[F,Cl,Br]',
    'B_peroxy_acid': '[CX3](=O)[OX2][OX2H]',
    # Tier C — 단순 인산/포스폰산 에스터 (DMMP 류). 배제하지 않고 기록만 한다.
    'C_phosphoryl_ester': '[PX4](=[O,S])([OX2][#6])[OX2,CX4]',
}
PAT = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}
assert all(p is not None for p in PAT.values()), 'SMARTS 컴파일 실패'


def alerts(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, []
    hit = [k for k, p in PAT.items() if m.HasSubstructMatch(p)]
    return Chem.MolToInchiKey(m), hit


def tier(hits):
    if any(h.startswith('A_') for h in hits):
        return 'A'
    if any(h.startswith('B_') for h in hits):
        return 'B'
    return 'C' if hits else '-'


# ---------------------------------------------- 자료 결합
C = pd.read_csv('descriptor_cache_xtb.csv')
FAST28 = pd.read_csv('outputs/xtb_fast_vs_slow.csv')          # 28종 빠른 모드
A = pd.read_csv('outputs/aegl_electronic.csv')
A = A[A.ok.fillna(False)].reset_index(drop=True)
# 동정 검증 실패 레코드 제외 (2026-10-11)
from cwa_registry import CWA as CWA_CAS_MAP, IDENTITY_FAILED
_nbad = int(A.cas.astype(str).isin(IDENTITY_FAILED).sum())
A = A[~A.cas.astype(str).isin(IDENTITY_FAILED)].reset_index(drop=True)
if _nbad:
    print(f'[pool] 동정 검증 실패 {_nbad}건 제외')
A['is_cwa'] = A.cas.astype(str).isin(CWA_CAS_MAP)

# 전자 서술자가 물리적으로 성립하지 않는 레코드 분리 (2026-10-11).
# 친전자성 지수는 omega = mu^2 / 2*eta 이고 eta 는 HOMO-LUMO 간극의 절반이다.
# 간극이 0 에 가까워지면 분모가 사라져 omega 가 발산한다. 이 풀에서
# 시안화나트륨은 간극 0.097 eV 에서 omega = 545 eV 를 받는다. 다음으로 작은
# 간극은 1.31 eV 이므로 이 한 건만 떨어져 있다. 간극이 이 수준이면 닫힌껍질
# 단일참조 GFN2 기술 자체가 성립하지 않으므로, 값이 극단적이라기보다 의미가
# 없다. 표준화된 거리를 쓰는 이상 한 건이 omega 축 전체를 압축한다
# (풀 표준편차 5.67 -> 35.24, 6배).
# 2026-10-12. 종전에는 간극 하한 하나만 두었고, 그 임계값은 적발된 한 건
# (시안화나트륨)에 사후적으로 맞춘 규칙이었다. 리뷰어 C 지적에 따라 물리적
# 타당성 심사 규칙을 사전에 선언하고 전부 적용한다. 네 규칙은 모두 "이 전자
# 서술자 집합이 그 화학종에 대해 정의되는가" 라는 하나의 질문에서 나온다.
GAP_MIN_EV = 0.5
SCREEN = {}


def _screen(df):
    import re as _re
    from rdkit import Chem as _Chem
    smi = df.smiles.astype(str)
    mols = [_Chem.MolFromSmiles(s) for s in smi]
    # R1 다성분 레코드: 조각 사이의 상대 배치가 임의이므로 쌍극자가 정의되지
    #    않고, 환산 분자량도 명명된 화합물의 것이 아닐 수 있다.
    r1 = smi.str.contains(r'\.', regex=True)
    # R2 개각 화학종: omega = mu^2/2eta 가 이 형태로 정의되지 않고, 궤도 파싱이
    #    alpha/beta 중 한쪽만 읽는다.
    r2 = pd.Series([bool(m) and any(a.GetNumRadicalElectrons() for a in m.GetAtoms())
                    for m in mols], index=df.index)
    # R3 단원자 화학종(수은 증기, 카드뮴, 적린 등): 결합이 없어 분자 서술자가
    #    의미를 갖지 않는다. 중원자 1개라도 수소가 붙은 HCl·HF·NH3·H2S 는
    #    통상적인 분자이므로 제외하지 않는다.
    r3 = pd.Series([bool(m) and m.GetNumHeavyAtoms() == 1
                    and sum(a.GetTotalNumHs() for a in m.GetAtoms()) == 0
                    for m in mols], index=df.index)
    # R4 간극 하한: eta -> 0 에서 omega 가 발산하고 닫힌껍질 단일참조 기술이
    #    성립하지 않는다.
    r4 = df.gap_eV < GAP_MIN_EV
    return pd.DataFrame(dict(multicomponent=r1, open_shell=r2,
                             single_atom=r3, tiny_gap=r4), index=df.index)


_S = _screen(A)
SCREEN = {k: int(v.sum()) for k, v in _S.items()}
_bad = _S.any(axis=1)
print(f'[pool] 물리적 타당성 심사 — 사전 선언한 4규칙으로 {int(_bad.sum())}건 제외')
for k, n in SCREEN.items():
    print(f'   {k:<16} {n:>3}건')
for r in A[_bad].head(6).itertuples():
    why = ','.join(k for k in SCREEN if _S.loc[r.Index, k])
    print(f'   예: {r.cas:12s} {str(r.aegl_name)[:28]:28s} 간극 {r.gap_eV:6.3f}  '
          f'omega {r.electrophilicity_index_eV:>10,.1f}  [{why}]')
A_unscreened = A.copy()
A = A[~_bad].reset_index(drop=True)

own = FAST28.merge(C[['compound_name', 'compound_class', 'canonical_smiles']],
                   on='compound_name')
own['smiles'] = own.canonical_smiles
own['label'] = own.compound_name
own['y_tox'] = np.nan
own['src'] = 'study28'
A['label'] = A.aegl_name
A['compound_class'] = 'aegl_pool'
A['y_tox'] = np.log10(A.aegl3_10min_mgm3)
A['src'] = 'aegl'
own['is_cwa'] = own.compound_class == 'real_agent'

cols = ['label', 'compound_class', 'smiles', 'y_tox', 'src', 'is_cwa'] + EL
P = pd.concat([own[cols], A[cols]], ignore_index=True)
ik, al = zip(*[alerts(s) for s in P.smiles])
P['inchikey'], P['alerts'] = ik, ['|'.join(a) for a in al]
P['tier'] = [tier(a) for a in al]
P = P[P.inchikey.notna()].reset_index(drop=True)

# 중복 제거: 28종 레코드를 우선 보존하되 AEGL 실측 독성값을 가져온다
P['skeleton'] = P.inchikey.str[:14]
tox_by_ik = (P[P.src == 'aegl'].dropna(subset=['y_tox'])
             .groupby('skeleton').y_tox.median())
P['y_tox'] = P.apply(lambda r: r.y_tox if r.y_tox == r.y_tox
                     else tox_by_ik.get(r.skeleton, np.nan), axis=1)
# 동일 skeleton 이 작용제로 표시된 레코드가 하나라도 있으면 전체를 작용제로 본다
_cwa_sk = set(P.loc[P.is_cwa.fillna(False), 'skeleton'])
P['is_cwa'] = P.skeleton.isin(_cwa_sk)
P['prio'] = (P.src == 'study28').astype(int)
P = (P.sort_values('prio', ascending=False)
       .drop_duplicates('skeleton', keep='first').reset_index(drop=True))
print(f'[pool] 중복 제거 후 {len(P)}종  '
      f"(연구 28종 {int((P.src=='study28').sum())} / AEGL 풀 {int((P.src=='aegl').sum())})")
print(f'[pool] 실측 독성 보유 {int(P.y_tox.notna().sum())}종')

# ---------------------------------------------- 충실도 축
AG = P[P.compound_class == 'real_agent']
print(f'[fidelity] 기준 진본 작용제 {len(AG)}종: ' + ', '.join(AG.label))
Z = P[EL].to_numpy(float)
mu, sd = Z.mean(0), Z.std(0, ddof=1)
Zz = (Z - mu) / sd
Az = Zz[(P.compound_class == 'real_agent').to_numpy()]
P['fidelity_dist'] = [np.sqrt(((Az - z) ** 2).sum(1)).min() for z in Zz]
# 작용제 자신은 거리 0 이므로 후보 풀에서 제외
Cand = P[(~P.is_cwa) & P.y_tox.notna()].reset_index(drop=True)
print(f'[H1] 후보 {len(Cand)}종 (군용작용제 {int(P.is_cwa.sum())}종 전량 제외, 실측 독성 보유)')

# ---------------------------------------------- H1 검정
x, y = Cand.fidelity_dist.to_numpy(), Cand.y_tox.to_numpy()
sw = stats.shapiro(y)
pr, pp = stats.pearsonr(x, y)
sr, sp = stats.spearmanr(x, y)


def boot_ci(fn):
    v = [fn(*[a[i] for a in (x, y)])
         for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    return np.percentile(v, [2.5, 97.5])


ci_p = boot_ci(lambda a, b: stats.pearsonr(a, b)[0])
ci_s = boot_ci(lambda a, b: stats.spearmanr(a, b)[0])
print(f'\n[H1] Shapiro(y) W={sw.statistic:.4f} p={sw.pvalue:.4f}')
print(f'[H1] Pearson  r = {pr:+.4f}  p = {pp:.4g}  95% CI [{ci_p[0]:+.3f}, {ci_p[1]:+.3f}]')
print(f'[H1] Spearman r = {sr:+.4f}  p = {sp:.4g}  95% CI [{ci_s[0]:+.3f}, {ci_s[1]:+.3f}]')
print(f'[H1] 결정계수 r^2 = {pr**2:.4f}  (충실도가 설명하는 독성 분산 {100*pr**2:.1f}%)')

verdict = ('H1 기각 — 충실도와 독성이 결합돼 있다 (가까울수록 독하다). RQ-3 재검토'
           if (sr < -0.3 and sp < 0.05) else
           'H1 유지 — 충실도와 안전성이 통계적으로 분리돼 있다. RQ-3 성립 가능')
print(f'[H1] 판정: {verdict}')

# ---------------------------------------------- 사분면 (RQ-2)
fid_cut = float(np.percentile(Cand.fidelity_dist, 25))     # 상위 25% 충실도
_TX = pd.read_csv('outputs/toxicity_aegl.csv')
_dcm = _TX.loc[_TX.compound_name == 'dichloromethane', 'AEGL3_10min_mgm3']
assert len(_dcm) == 1 and _dcm.iloc[0] == _dcm.iloc[0], 'DCM AEGL-3 미확보'
DCM_AEGL3 = float(_dcm.iloc[0])            # 하드코딩 금지: 파싱 산출물에서 직접 읽음
tox_cut = float(np.log10(DCM_AEGL3))       # DCM(현행 모사체) 기준선
Cand['quadrant'] = np.where(
    (Cand.fidelity_dist <= fid_cut) & (Cand.y_tox >= tox_cut), 'faithful_and_safe',
    np.where(Cand.fidelity_dist <= fid_cut, 'faithful_but_toxic',
             np.where(Cand.y_tox >= tox_cut, 'safe_but_unfaithful', 'neither')))
Cand['alert_flag'] = Cand.tier.isin(['A', 'B'])   # Tier A/B 만 배제
print(f'\n[RQ-2] 충실도 경계 = 상위 25퍼센타일 거리 {fid_cut:.3f}')
print(f'[RQ-2] 안전 경계 = 현행 모사체 DCM 의 AEGL-3 {DCM_AEGL3:,.0f} mg/m3 (log {tox_cut:.3f})')
print(Cand.quadrant.value_counts().to_string())
FS = Cand[(Cand.quadrant == 'faithful_and_safe') & ~Cand.alert_flag] \
    .sort_values('fidelity_dist')
print(f'\n[RQ-2] 충실+안전 사분면, 구조경보 통과: {len(FS)}종')
print(FS[['label', 'fidelity_dist', 'y_tox', 'src']].head(25)
      .to_string(index=False, float_format=lambda v: f'{v:.3f}'))

# ---------------------------------------------- 구조 경보 검증
print('\n[gate] 구조 경보 검증 (연구 28종)')
S = P[P.src == 'study28']
for _, r in S[S.tier != '-'].sort_values('tier').iterrows():
    print(f'  Tier {r.tier}  {r.label:<24} {r.alerts}')
ra = S[S.compound_class == 'real_agent']
missA = ra[ra.tier != 'A']
print(f'  진본 작용제 {len(ra)}종 중 Tier A {int((ra.tier=="A").sum())}종, '
      f'미검출 {len(missA)}종: ' + (', '.join(missA.label) if len(missA) else '없음'))
fp = S[(S.compound_class == 'interferent') & S.tier.isin(['A', 'B'])]
print(f'  간섭물 Tier A/B 오경보 {len(fp)}종'
      + (f": {', '.join(fp.label)}" if len(fp) else ''))
print(f'  모사체 판정: ' + ', '.join(
    f'{r.label}=Tier {r.tier}' for _, r in S[S.compound_class == 'simulant'].iterrows()))

# 안전 기준 스윕. 본문이 쓰는 수치는 전부 여기서 산출해 JSON 에 남긴다.
# 종전에는 완화 기준의 후보 수를 도판 스크립트(make_fig_h1.py)에서 계산했고,
# 분석 산출물에 기록되지 않아 본문 수치를 재현할 경로가 없었다.
TIERS = [
    ('ten_times_safer', tox_cut + 1.0, 'DCM 대비 10배 안전'),
    ('incumbent',       tox_cut,       'DCM 와 동등 이상 안전'),
    ('relaxed_10',      tox_cut - 1.0, 'DCM 의 1/10 (10배 독성)까지 허용'),
    ('relaxed_100',     tox_cut - 2.0, 'DCM 의 1/100 (100배 독성)까지 허용'),
]
tier_out = {}
print('\n[RQ-2] 안전 기준 스윕 (충실도 경계 고정, 구조경보 Tier A/B 배제 후)')
for key, bar, desc in TIERS:
    sel = Cand[(Cand.fidelity_dist <= fid_cut) & (Cand.y_tox >= bar)
               & ~Cand.alert_flag].sort_values('fidelity_dist')
    # 현행 모사체보다 충실한가 (DCM 의 충실도 거리와 비교)
    _d = Cand[Cand.label == 'dichloromethane'].fidelity_dist
    dcm_fid = float(_d.iloc[0]) if len(_d) else None
    more_faithful = (int((sel.fidelity_dist < dcm_fid).sum())
                     if dcm_fid is not None else None)
    tier_out[key] = dict(bar_log10=round(float(bar), 4),
                         bar_mgm3=round(float(10 ** bar), 1),
                         description=desc, n=int(len(sel)),
                         n_more_faithful_than_dcm=more_faithful,
                         labels=sel.label.tolist()[:30])
    print(f'  {desc:<34} AEGL-3 >= {10**bar:>12,.0f} mg/m3  '
          f'후보 {len(sel):>2}종'
          + (f'  (그중 DCM 보다 충실 {more_faithful}종)'
             if more_faithful is not None else ''))
    if len(sel):
        print(sel[['label', 'fidelity_dist', 'y_tox', 'tier']].head(20)
              .to_string(index=False, float_format=lambda v: f'{v:.3f}'))
FS2 = Cand[(Cand.fidelity_dist <= fid_cut) & (Cand.y_tox >= tox_cut + 1.0)
           & ~Cand.alert_flag].sort_values('fidelity_dist')

json.dump({r.label: r.tier for _, r in P[P.src == 'study28'].iterrows()},
          open('outputs/h1_gate_tiers.json', 'w'), indent=2, ensure_ascii=False)
Cand.to_csv('outputs/h1_candidates.csv', index=False)
json.dump(dict(seed=SEED, bootstrap=N_BOOT, n_pool=int(len(P)),
               n_candidates=int(len(Cand)),
               pearson_r=float(pr), pearson_p=float(pp),
               pearson_ci=[float(ci_p[0]), float(ci_p[1])],
               spearman_r=float(sr), spearman_p=float(sp),
               spearman_ci=[float(ci_s[0]), float(ci_s[1])],
               r_squared=float(pr ** 2), shapiro_p=float(sw.pvalue),
               verdict=verdict, fidelity_cut=fid_cut, tox_cut=tox_cut,
               dcm_aegl3_mgm3=DCM_AEGL3,
               quadrants={k: int(v) for k, v in Cand.quadrant.value_counts().items()},
               n_faithful_safe_passing_alerts=int(len(FS)),
               safety_tiers=tier_out,
               pool_screen=SCREEN,
               structural_alerts=ALERTS),
          open('outputs/h1_summary.json', 'w'), indent=2, ensure_ascii=False)
print('\n[done] outputs/h1_candidates.csv, h1_summary.json')
