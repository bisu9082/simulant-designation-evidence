#!/usr/bin/env python3
"""
AEGL-3 (10분, mg/m3) 예측 QSAR — simulant 안전성 축 보간용.

이중용도 가드 (설계에 내장):
  G-1 화학작용제(신경작용제·수포작용제)는 학습셋에서 전량 제외한다.
      이들은 held-out 평가점으로만 쓴다.
  G-2 적용영역(applicability domain) 밖 화합물에는 수치를 부여하지 않고
      '[적용영역 밖]' 으로 표기한다. 신규 작용제 서열화 용도로 쓰일 수 없게 하는
      장치이자, 동시에 과학적으로 올바른 제한이다.
  G-3 모델은 '후보가 용매 수준으로 안전한가'를 거르는 floor test 로만 해석한다.
      독성 순위를 매기는 용도로 쓰지 않는다.

통계 기준 (프로젝트 GATE): 시드 5개 × 5-fold 반복 CV, bootstrap 10,000, 95% CI,
baseline 2종 대비.
"""
import json
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
RDLogger.DisableLog('rdApp.*')

SEEDS = [11, 22, 33, 44, 55]          # 5개 (GATE: 시드 >= 5)
N_BOOT = 10000
rng = np.random.default_rng(42)

# ---- G-1 : 학습 제외 대상 (CWA) — CAS 명부. 문자열 매칭 금지 (2026-10-11 교정)
from cwa_registry import CWA, IDENTITY_FAILED


def desc(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    d = {}
    d['MW'] = Descriptors.MolWt(m)
    d['logP'] = Crippen.MolLogP(m)
    d['MR'] = Crippen.MolMR(m)
    d['TPSA'] = rdMolDescriptors.CalcTPSA(m)
    d['HBA'] = rdMolDescriptors.CalcNumHBA(m)
    d['HBD'] = rdMolDescriptors.CalcNumHBD(m)
    d['rot'] = rdMolDescriptors.CalcNumRotatableBonds(m)
    d['rings'] = rdMolDescriptors.CalcNumRings(m)
    d['arom'] = rdMolDescriptors.CalcNumAromaticRings(m)
    d['heavy'] = m.GetNumHeavyAtoms()
    d['fcsp3'] = rdMolDescriptors.CalcFractionCSP3(m)
    d['hetero'] = rdMolDescriptors.CalcNumHeteroatoms(m)
    d['chi1v'] = rdMolDescriptors.CalcChi1v(m)
    d['chi2v'] = rdMolDescriptors.CalcChi2v(m)
    d['kappa2'] = rdMolDescriptors.CalcKappa2(m)
    d['labute'] = rdMolDescriptors.CalcLabuteASA(m)
    d['balaban'] = Descriptors.BalabanJ(m)
    d['bertz'] = Descriptors.BertzCT(m)
    for sym, key in [('F', 'nF'), ('Cl', 'nCl'), ('Br', 'nBr'), ('N', 'nN'),
                     ('O', 'nO'), ('S', 'nS'), ('P', 'nP'), ('As', 'nAs'),
                     ('Si', 'nSi')]:
        d[key] = sum(a.GetSymbol() == sym for a in m.GetAtoms())
    for key, sma in [('n_CN_nitrile', '[C]#[N]'), ('n_carbonyl', '[CX3]=[OX1]'),
                     ('n_nitro', '[N+](=O)[O-]'), ('n_epox', 'C1OC1'),
                     ('n_PdO', '[P]=[O]'), ('n_halo_C', '[#6][F,Cl,Br,I]'),
                     ('n_SH', '[SX2H]'), ('n_OH', '[OX2H]'),
                     ('n_isocyanate', 'N=C=O'), ('n_acidhal', '[CX3](=O)[F,Cl,Br]')]:
        p = Chem.MolFromSmarts(sma)
        d[key] = len(m.GetSubstructMatches(p)) if p else 0
    return d


D = pd.read_csv('outputs/aegl_trainset_raw.csv')
rows, bad = [], []
for r in D.itertuples():
    d = desc(r.smiles)
    if d is None:
        bad.append(r.cas); continue
    d.update(cas=r.cas, aegl_name=r.aegl_name, section=r.section,
             y=np.log10(r.aegl3_10min_mgm3), aegl3=r.aegl3_10min_mgm3)
    rows.append(d)
A = pd.DataFrame(rows)
FEAT = [c for c in A.columns if c not in ('cas', 'aegl_name', 'section', 'y', 'aegl3')]
A = A.replace([np.inf, -np.inf], np.nan).dropna(subset=FEAT).reset_index(drop=True)
print(f'[data] AEGL 화합물 {len(A)}건, 서술자 {len(FEAT)}개  '
      f'(RDKit 파싱 실패 {len(bad)})')

bad_id = A.cas.astype(str).isin(IDENTITY_FAILED)
if bad_id.any():
    print(f'[G-0] 동정 검증 실패로 제외 {int(bad_id.sum())}건: '
          + ', '.join(f'{r.cas} {r.aegl_name}' for r in A[bad_id].itertuples()))
A = A[~bad_id].reset_index(drop=True)
is_cwa = A.cas.astype(str).isin(CWA)
TR, HO = A[~is_cwa].reset_index(drop=True), A[is_cwa].reset_index(drop=True)
print(f'[G-1] 학습 {len(TR)}건 / CWA held-out {len(HO)}건: '
      f"{', '.join(HO.aegl_name)}")
print(f'[data] y = log10(AEGL-3 10min mg/m3)  범위 {TR.y.min():.2f} ~ {TR.y.max():.2f}'
      f'  (중앙 {TR.y.median():.2f})')

X, y = TR[FEAT].to_numpy(float), TR.y.to_numpy(float)


def cv_scores(make_model, X, y, seeds=SEEDS, k=5):
    oof_all = []
    for s in seeds:
        oof = np.empty_like(y)
        for tr, te in KFold(k, shuffle=True, random_state=s).split(X):
            sc = StandardScaler().fit(X[tr])
            m = make_model(s).fit(sc.transform(X[tr]), y[tr])
            oof[te] = m.predict(sc.transform(X[te]))
        oof_all.append(oof)
    return np.array(oof_all)


def summarise(oof_all, y, tag):
    r2 = [1 - ((o - y) ** 2).sum() / ((y - y.mean()) ** 2).sum() for o in oof_all]
    rmse = [float(np.sqrt(((o - y) ** 2).mean())) for o in oof_all]
    mae = [float(np.abs(o - y).mean()) for o in oof_all]
    pooled = oof_all.mean(0)
    res = (pooled - y) ** 2
    boot = [np.sqrt(rng.choice(res, res.size, replace=True).mean())
            for _ in range(N_BOOT)]
    ci = np.percentile(boot, [2.5, 97.5])
    print(f'  {tag:<26} R2 {np.mean(r2):.3f} +/- {np.std(r2, ddof=1):.3f}   '
          f'RMSE {np.mean(rmse):.3f} +/- {np.std(rmse, ddof=1):.3f}   '
          f'MAE {np.mean(mae):.3f}   RMSE 95% CI [{ci[0]:.3f}, {ci[1]:.3f}]')
    return dict(tag=tag, r2_mean=float(np.mean(r2)), r2_sd=float(np.std(r2, ddof=1)),
                rmse_mean=float(np.mean(rmse)), rmse_sd=float(np.std(rmse, ddof=1)),
                mae_mean=float(np.mean(mae)),
                rmse_ci95=[float(ci[0]), float(ci[1])], n_seeds=len(oof_all))


print(f'\n[CV] 시드 {SEEDS} × 5-fold, bootstrap {N_BOOT}')
ia = [FEAT.index('logP'), FEAT.index('MW')]
results = []
results.append(summarise(np.array([np.full_like(y, y.mean())] * len(SEEDS)),
                         y, 'baseline-1 평균예측'))
results.append(summarise(cv_scores(lambda s: LinearRegression(), X[:, ia], y),
                         y, 'baseline-2 logP+MW 선형'))
oof_rf = cv_scores(lambda s: RandomForestRegressor(
    n_estimators=600, min_samples_leaf=2, max_features=0.4,
    random_state=s, n_jobs=-1), X, y)
results.append(summarise(oof_rf, y, 'RandomForest (전체 서술자)'))

# out-of-fold 예측을 내보낸다. 꼬리 회귀 진단(minority_regression.py)이 이 파일을
# 읽는다. 진단 쪽에서 모델을 다시 적합하면 같은 수치가 두 경로에서 나오게 되고,
# 이 파이프라인은 이미 그 중복으로 본문 수치가 갈린 적이 있다.
pd.DataFrame(dict(cas=TR.cas, aegl_name=TR.aegl_name, section=TR.section,
                  y=y, aegl3=TR.aegl3,
                  oof_pred=oof_rf.mean(0))).to_csv(
    'outputs/qsar_oof.csv', index=False)
print(f'[oof] outputs/qsar_oof.csv  {len(TR)}행 '
      f'(시드 {len(SEEDS)}개 평균 out-of-fold 예측)')

# ------------------------------------------------- G-2 적용영역
sc_all = StandardScaler().fit(X)
Xz = sc_all.transform(X)


def knn_dist(Q, ref=Xz, k=5):
    out = []
    for q in Q:
        d = np.sort(np.sqrt(((ref - q) ** 2).sum(1)))
        out.append(d[1:k + 1].mean() if np.isclose(d[0], 0) else d[:k].mean())
    return np.array(out)


tr_knn = knn_dist(Xz)
AD_T = float(np.percentile(tr_knn, 95))
print(f'\n[G-2] 적용영역 경계 = 학습셋 자체 kNN(k=5) 거리의 95퍼센타일 = {AD_T:.3f}')

rf_full = RandomForestRegressor(n_estimators=600, min_samples_leaf=2,
                                max_features=0.4, random_state=11,
                                n_jobs=-1).fit(Xz, y)

# held-out CWA 평가
if len(HO):
    Hz = sc_all.transform(HO[FEAT].to_numpy(float))
    hp, hd = rf_full.predict(Hz), knn_dist(Hz)
    print('\n[G-1 검증] 학습에서 제외한 CWA 의 예측 (모델이 작용제를 맞히는가?)')
    for n, t, p, dd in zip(HO.aegl_name, HO.y, hp, hd):
        print(f'  {n:<34} 실측 log10 {t:+.2f}  예측 {p:+.2f}  '
              f'오차 {p-t:+.2f} log  AD거리 {dd:.2f}'
              f'{"  [적용영역 밖]" if dd > AD_T else ""}')
    # held-out 예측을 파일로 남긴다. 종전에는 콘솔에만 찍혀, 본문 3.5절의
    # 헤드라인 세 개(오차 범위, 배수, 적용영역 내 비율)를 확인하려면 코드를
    # 재실행해야 했다 (GATE 8 재현성 지적).
    _E = pd.DataFrame(dict(
        cas=HO.cas.to_numpy(), aegl_name=HO.aegl_name.to_numpy(),
        measured_log10=HO.y.to_numpy(), measured_mgm3=HO.aegl3.to_numpy(),
        predicted_log10=np.asarray(hp), predicted_mgm3=10 ** np.asarray(hp),
        error_log10=np.asarray(hp) - HO.y.to_numpy(),
        fold_too_safe=10 ** (np.asarray(hp) - HO.y.to_numpy()),
        ad_distance=np.asarray(hd),
        in_domain=np.asarray(hd) <= AD_T))
    _E.sort_values('error_log10').to_csv('outputs/qsar_heldout.csv', index=False)
    _e = _E.error_log10
    print(f'  [heldout] outputs/qsar_heldout.csv  {len(_E)}행 | '
          f'오차 {_e.min():+.2f} ~ {_e.max():+.2f} log '
          f'({10**_e.min():.1f} ~ {float("%.2g" % (10**_e.max())):,.0f}배) | '
          f'전부 안전 쪽 {bool((_e > 0).all())} | '
          f'적용영역 내 {int(_E.in_domain.sum())}/{len(_E)} '
          f'(그중 최대오차 {_E[_E.in_domain].error_log10.max():+.2f} log)')
    print(f'  -> CWA 평균 절대오차 {np.abs(hp - HO.y.to_numpy()).mean():.2f} log 단위, '
          f'적용영역 밖 {int((hd > AD_T).sum())}/{len(HO)}')

# ------------------------------------------------- 28종 예측
C = pd.read_csv('descriptor_cache_xtb.csv')
T = pd.read_csv('outputs/toxicity_aegl.csv')
qrows = []
for r in C.itertuples():
    d = desc(r.canonical_smiles)
    d.update(compound_name=r.compound_name, compound_class=r.compound_class)
    qrows.append(d)
Q = pd.DataFrame(qrows)
Qz = sc_all.transform(Q[FEAT].to_numpy(float))
Q['pred_log10'] = rf_full.predict(Qz)
Q['ad_dist'] = knn_dist(Qz)
Q['in_domain'] = Q.ad_dist <= AD_T
Q = Q.merge(T[['compound_name', 'AEGL3_10min_mgm3']], on='compound_name', how='left')
Q['pred_mgm3'] = np.where(Q.in_domain, 10 ** Q.pred_log10, np.nan)
Q['status'] = np.where(Q.AEGL3_10min_mgm3.notna(), '실측 AEGL-3',
                       np.where(Q.in_domain, '예측(적용영역 내)', '[적용영역 밖]'))
# G-3: DMMP 의 예측은 본문에서 쓰지 않는다. 적용영역 거리가 작용제들 사이에
# 있어 모델이 한두 자릿수 틀리는 구간이기 때문이다. 산출 CSV 에도 그렇게 적어
# 둔다 (종전에는 in_domain=True 로만 보여 본문과 어긋났다).
Q.loc[Q.compound_name == 'DMMP', 'status'] = '[예측 보류: 작용제 영역]'
out = Q[['compound_name', 'compound_class', 'AEGL3_10min_mgm3', 'pred_log10',
         'pred_mgm3', 'ad_dist', 'in_domain', 'status']]
out.to_csv('outputs/toxicity_predicted.csv', index=False)
print('\n[예측] 28종')
print(out.sort_values('ad_dist').to_string(index=False,
      float_format=lambda v: f'{v:.3g}'))
print(f"\n[예측] 적용영역 내 {int(Q.in_domain.sum())}/28, "
      f"밖 {int((~Q.in_domain).sum())}/28")

json.dump(dict(seeds=SEEDS, bootstrap=N_BOOT, n_train=int(len(TR)),
               n_heldout_cwa=int(len(HO)), n_features=len(FEAT),
               ad_threshold=AD_T, cv=results,
               guards=['G-1 CWA 학습 제외', 'G-2 적용영역 밖 수치 미부여',
                       'G-3 floor test 해석 한정']),
          open('outputs/qsar_tox_summary.json', 'w'), indent=2, ensure_ascii=False)
print('\n[done] outputs/toxicity_predicted.csv, qsar_tox_summary.json')
