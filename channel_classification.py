#!/usr/bin/env python3
"""유효 채널 수에 대한 분류 성능 민감도 (2026-11-02 신설, T2-7).

왜 필요한가. 원고는 2차 런의 교차검증 분류 정확도 0.964 를 인용해 '어레이
자체는 문제가 아니다' 라고 적는다. 그 값은 29 염료를 선언한 어레이에서 나온
것인데, 이 저장소의 중복 채널 축약은 2차 런 자외 채널의 유효 열이 23 개임을
보인다 (Dye8 과 Dye14 가 전부 0, Dye22 가 Dye23~27 과 동일). 선언 29 와
유효 23 의 차이가 분류 결론을 바꾸는지 확인하지 않고 두 수를 같은 문장에
쓰는 것은 맞지 않다.

할 수 없는 것부터 적는다. 0.964 자체는 이 저장소의 데이터로 재현할 수 없다.
예치된 피벗은 (물질, 농도) 한 행이고 반복측정이 이미 평균되어 있어 원논문의
교차검증 설계(반복측정을 표본으로 쓰는)를 복원할 방법이 없다. 따라서
[수집 불가: 반복측정 수준 데이터] 로 적고, 인용값을 우리 수치로 대체하지
않는다 ([C-3]).

할 수 있는 것은 같은 질문의 약한 판본이다. 0 µM 을 뺀 네 농도점을 표본으로
놓으면 2차 런은 7급 x 4표본 = 28 표본이 된다. 여기에 전체교차검증(LOO)을
걸어 29 선언 채널과 23 유효 채널의 정확도를 나란히 낸다. 표본이 28 개뿐이고
농도점을 독립 표본으로 쓰므로 이 수치는 원논문의 0.964 와 같은 양이 아니다.
그 점을 산출물에 명시한다.

분류기는 두 가지를 쓴다. 하나는 이 논문의 나머지가 쓰는 것과 같은 전처리를
따르는 최근접 중심(단위정규화 패턴의 코사인), 다른 하나는 축소추정 LDA 다.
표본 28 개에 변수 29 개이므로 정칙화 없는 LDA 는 성립하지 않는다.

사용: python3 channel_classification.py
"""
import json

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

import deltaE_premise as m

CONC_NZ = [c for c in m.CONC if c != 0]
QUOTED = {'second_run_accuracy': 0.964, 'second_run_macro_f1': 0.964,
          'third_run_accuracy': 1.000,
          'reproducible_from_deposit': False,
          'note': '[수집 불가: 반복측정 수준 ΔE 가 예치물에 없다. 피벗은 '
                  '(물질, 농도) 한 행이며 반복측정 평균이다. 원논문의 '
                  '교차검증 설계를 복원할 수 없으므로 인용값을 재현하거나 '
                  '대체하지 않는다.]'}


def unit_rows(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    n[n == 0] = 1.0
    return X / n


def nearest_centroid_loo(X, y):
    """단위정규화 패턴의 코사인 최근접 중심. 본문 나머지와 같은 전처리."""
    X = unit_rows(X)
    ok = 0
    for tr, te in LeaveOneOut().split(X):
        cl = {}
        for c in np.unique(y[tr]):
            v = X[tr][y[tr] == c].mean(0)
            nv = np.linalg.norm(v)
            cl[c] = v / nv if nv else v
        s = {c: float(np.dot(X[te][0], v)) for c, v in cl.items()}
        ok += int(max(s, key=s.get) == y[te][0])
    return ok / len(y)


def lda_cv(X, y, folds=None):
    """축소추정 LDA. folds=None 이면 LOO."""
    sp = LeaveOneOut().split(X) if folds is None else StratifiedKFold(
        n_splits=folds, shuffle=True, random_state=42).split(X, y)
    ok = tot = 0
    for tr, te in sp:
        sc = StandardScaler().fit(X[tr])
        md = LinearDiscriminantAnalysis(solver='lsqr', shrinkage='auto')
        md.fit(sc.transform(X[tr]), y[tr])
        pr = md.predict(sc.transform(X[te]))
        ok += int((pr == y[te]).sum()); tot += len(te)
    return ok / tot


def cell(tag, mode):
    D = pd.read_csv(f'deltaE/{tag}.{mode}_deltaE_pivot.csv')
    declared = [c for c in D.columns if c.lower().startswith('dye')]
    effective = m.collapse_duplicate_channels(D.copy(), list(declared),
                                              f'{tag}.{mode}')
    A = D[D.Concentration.isin(CONC_NZ)]
    y = A.Agent.to_numpy()
    out = {'n_samples': int(len(A)), 'n_classes': int(len(set(y))),
           'n_declared_channels': len(declared),
           'n_effective_channels': len(effective),
           'dropped_channels': sorted(set(declared) - set(effective))}
    for nm, cols in (('declared', declared), ('effective', effective)):
        X = A[cols].to_numpy(float)
        out[nm] = {
            'nearest_centroid_loo': round(nearest_centroid_loo(X, y), 4),
            'lda_shrinkage_loo': round(lda_cv(X, y), 4),
        }
        if len(A) >= 4 * out['n_classes']:
            out[nm]['lda_shrinkage_4fold'] = round(lda_cv(X, y, folds=4), 4)
    return out


def main():
    print('=' * 78)
    print('[T2-7] 선언 채널 대 유효 채널의 분류 성능')
    print('=' * 78)
    print('\n인용값 0.964 / 1.000 은 재현 대상이 아니다:')
    print('  ' + QUOTED['note'])
    print('\n아래는 예치 데이터로 가능한 약한 판본이다. 농도점을 표본으로 쓴다.')

    rep = {'quoted': QUOTED, 'cells': {}}
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            r = cell(tag, mode)
            rep['cells'][f'{tag}.{mode}'] = r
            dd = (f"제거 {','.join(r['dropped_channels'])}"
                  if r['dropped_channels'] else '제거 없음')
            print(f"\n  {tag}.{mode}  표본 {r['n_samples']} 급 {r['n_classes']}  "
                  f"선언 {r['n_declared_channels']} -> 유효 "
                  f"{r['n_effective_channels']}  ({dd})")
            for nm in ('declared', 'effective'):
                e = r[nm]
                extra = (f"  4겹LDA {e['lda_shrinkage_4fold']:.4f}"
                         if 'lda_shrinkage_4fold' in e else '')
                print(f"     {nm:<10} 최근접중심LOO {e['nearest_centroid_loo']:.4f}"
                      f"   축소LDA LOO {e['lda_shrinkage_loo']:.4f}{extra}")

    json.dump(rep, open('outputs/channel_classification.json', 'w'), indent=1,
              ensure_ascii=False)
    print('\n-> outputs/channel_classification.json')


if __name__ == '__main__':
    main()
