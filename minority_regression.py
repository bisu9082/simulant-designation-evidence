#!/usr/bin/env python3
"""꼬리 회귀 진단 — 모델이 독성 극단에 닿지 못하는 이유를 수치로 보인다.

본문 3.5절이 쓰는 수치를 산출한다. 종전에는 outputs/minority_regression.json
만 존재하고 생성 스크립트가 패키지에 없었다 (GATE 8 3차 지적). 여기서 재현
경로를 복원한다.

모델을 다시 적합하지 않고 qsar_tox.py 가 내보낸 out-of-fold 예측
(outputs/qsar_oof.csv) 을 읽는다. 같은 수치를 두 경로에서 계산하면 갈라진다는
것을 이 파이프라인에서 이미 겪었다 (분자량 84.93/84.933, cyclosarin 이름 매칭).

검정하는 것
  1) 학습세트에 아세틸콜린에스터라제 억제 기작이 실제로 들어 있는가.
     '기작이 없어서 작용제를 못 맞힌다' 가설의 반증.
  2) 전체 out-of-fold 편향은 0 근처인데 독성 최하위 십분위에서만 안전 쪽으로
     쏠리는가.
  3) 억제제 부분집합 안에서도 같은 꼬리 회귀가 나타나는가.
"""
import json

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

# 아세틸콜린에스터라제 억제 약리단을 구조로 판정한다. 손으로 쓴 이름·CAS 목록을
# 쓰지 않는다. 종전 본문은 '11종(유기인 8 + 카바메이트 3)' 이라고만 적었고 어느
# 화합물인지 밝히지 않아 독자가 검증할 수 없었다. 여기서는 SMARTS 로 뽑고
# 적중 화합물을 전부 이름과 수치로 출력한다.
#   유기인 에스터: 이탈기를 가진 포스포릴/티오포스포릴 에스터
#   카바메이트 에스터: N-H 카바메이트 (옥심 카바메이트 포함)
SMARTS = {
    'organophosphorus_ester':
        '[PX4](=[OX1,SX1])([OX2,SX2,NX3])([OX2,SX2,NX3,CX4])[OX2,SX2,F,$(C#N)]',
    'carbamate_ester': '[NX3;H1,H2][CX3](=[OX1])[OX2]',
}

O = pd.read_csv('outputs/qsar_oof.csv')
O['cas'] = O.cas.astype(str).str.strip()
O['err'] = O.oof_pred - O.y          # 양수 = 실제보다 안전하게 예측
_T = pd.read_csv('outputs/aegl_trainset_raw.csv')
_T['cas'] = _T.cas.astype(str).str.strip()
O = O.merge(_T[['cas', 'smiles']], on='cas', how='left')

_mols = {r.cas: Chem.MolFromSmiles(str(r.smiles)) for r in O.itertuples()}
groups = {}
for key, sma in SMARTS.items():
    pat = Chem.MolFromSmarts(sma)
    groups[key] = [c for c, m in _mols.items()
                   if m is not None and pat is not None
                   and m.HasSubstructMatch(pat)]
ache_cas = sorted(set(sum(groups.values(), [])))
ache = O[O.cas.isin(ache_cas)]
n = len(O)

print(f'[train] 학습 {n}건')
print(f'[ache]  아세틸콜린에스터라제 억제 약리단 보유 {len(ache)}건 = '
      f'학습세트의 {100*len(ache)/n:.1f}%')
for key, cs in groups.items():
    sub = O[O.cas.isin(cs)].sort_values('aegl3')
    print(f'  {key} {len(cs)}건')
    for r in sub.itertuples():
        print(f'    {r.cas:12s} {str(r.aegl_name)[:30]:30s} '
              f'AEGL-3 {r.aegl3:>10,.3f} mg/m3  오차 {r.err:+.2f} log')

SARIN = 0.366708
worse = ache[ache.aegl3 < SARIN].sort_values('aegl3')
print(f'\n[ache]  사린({SARIN:.3f} mg/m3)보다 급성독성이 강한 억제제 {len(worse)}건: '
      + (', '.join(f'{r.aegl_name} {r.aegl3:.3f}' for r in worse.itertuples())
         or '없음'))

print(f'\n[bias]  전체 out-of-fold 편향  {O.err.mean():+.3f} log')
print(f'[bias]  억제제 부분집합         {ache.err.mean():+.3f} log')
print(f'[bias]  그 외                   {O[~O.cas.isin(ache_cas)].err.mean():+.3f} log')

cut = O.y.quantile(0.10)
low = O[O.y <= cut]
hi = O[O.y >= O.y.quantile(0.90)]
# 본문에는 분위수 경계가 아니라 그 구간에 실제로 들어간 값의 극값을 쓴다.
# 종전에는 경계값 3.6 / 7871 을 적었는데 실제 최대/최소는 3.21 / 8094 였다.
low_max = float(10 ** low.y.max())
hi_min = float(10 ** hi.y.min())
print(f'\n[tail]  독성 최상위 십분위 (구간 내 최대 AEGL-3 {low_max:.2f} mg/m3, n={len(low)}): '
      f'평균 편향 {low.err.mean():+.3f} log, 안전 쪽 {100*(low.err>0).mean():.0f}%')
print(f'[tail]  독성 최하위 십분위 (구간 내 최소 AEGL-3 {hi_min:,.0f} mg/m3, '
      f'n={len(hi)}): 평균 편향 {hi.err.mean():+.3f} log, '
      f'안전 쪽 {100*(hi.err>0).mean():.0f}%')

a = ache.sort_values('aegl3')
most, least = a.iloc[0], a.iloc[-1]
print(f'\n[tail]  억제제 부분집합 내부')
print(f'        가장 독성 강한 {most.aegl_name}: '
      f'{10**most.err:,.0f}배 안전하게 예측 ({most.err:+.3f} log)')
print(f'        가장 독성 약한 {least.aegl_name}: '
      f'{10**(-least.err):,.0f}배 독하게 예측 ({least.err:+.3f} log)')

out = dict(
    n_train=int(n),
    smarts=SMARTS,
    n_ache=int(len(ache)),
    n_by_group={k: int(len(v)) for k, v in groups.items()},
    frac_ache=float(len(ache) / n),
    ache_compounds=[dict(cas=r.cas, name=r.aegl_name,
                         aegl3_mgm3=float(r.aegl3), err_log=float(r.err))
                    for r in a.itertuples()],
    sarin_aegl3_mgm3=SARIN,
    n_ache_more_toxic_than_sarin=int(len(worse)),
    ache_more_toxic_than_sarin=[dict(name=r.aegl_name, aegl3_mgm3=float(r.aegl3))
                                for r in worse.itertuples()],
    err_all=float(O.err.mean()), err_ache=float(ache.err.mean()),
    err_rest=float(O[~O.cas.isin(ache_cas)].err.mean()),
    low10_cut_mgm3=float(10 ** cut), low10_max_mgm3=low_max,
    low10_n=int(len(low)),
    low10_err=float(low.err.mean()), low10_safe_frac=float((low.err > 0).mean()),
    high10_cut_mgm3=float(10 ** O.y.quantile(0.90)), high10_min_mgm3=hi_min,
    high10_n=int(len(hi)),
    high10_err=float(hi.err.mean()), high10_safe_frac=float((hi.err > 0).mean()),
    most_toxic_ache=dict(name=most.aegl_name, err_log=float(most.err),
                         fold_too_safe=float(10 ** most.err)),
    least_toxic_ache=dict(name=least.aegl_name, err_log=float(least.err),
                          fold_too_toxic=float(10 ** (-least.err))),
    source='outputs/qsar_oof.csv (qsar_tox.py 가 내보낸 out-of-fold 예측)',
)
json.dump(out, open('outputs/minority_regression.json', 'w'),
          indent=2, ensure_ascii=False)
print('\n[done] outputs/minority_regression.json')
