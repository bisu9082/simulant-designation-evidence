#!/usr/bin/env python3
"""서술자 캐시의 3D 열을 앙상블 평균본으로 교체한다 (2026-10-25 신설).

왜 전체 캐시를 다시 만들지 않는가. xtb 전자 서술자(E_HOMO, E_LUMO, 쌍극자,
Wiberg, Fukui)는 GFN2 계산 결과이고 이 환경에는 xtb 가 없다. 그 열을 다시
만들 수 없으므로, 재현 가능한 RDKit 3D 열 10개만 제자리에서 교체한다.

교체 대상: build_descriptors.ensemble_3d() 가 내는 10열 + n_conformers +
conformer_w_max. 나머지 열은 손대지 않는다.

--check 로 부르면 바꿀 폭만 보고하고 쓰지 않는다.
"""
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

import build_descriptors as B

RDLogger.DisableLog('rdApp.*')

CACHES = ['descriptor_cache.csv', 'descriptor_cache_35.csv',
          'descriptor_cache_xtb.csv', 'descriptor_cache_xtb_35.csv']


def main():
    apply = '--check' not in sys.argv
    for path in CACHES:
        df = pd.read_csv(path)
        # input_smiles 를 쓴다. descriptors() 와 verify() 가 그 열을 쓰므로
        # canonical_smiles 로 임베딩하면 원자 순서가 달라져 같은 시드에서도
        # 다른 형태이성질체 집합이 나온다. 2026-10-25 에 그 경로로 재현
        # 검사가 계속 실패했다.
        smi_col = 'input_smiles'
        rows = []
        for smi in df[smi_col]:
            m = Chem.MolFromSmiles(smi)
            mh, cids, ff = B.embed(m)
            ens = B.ensemble_3d(mh, ff)
            ens['n_conformers'] = mh.GetNumConformers()
            rows.append(ens)
        worst = []
        for c in B.D3_NAMES:
            if c not in df.columns:
                continue
            old = df[c].to_numpy(float)
            new = np.array([r[c] for r in rows], float)
            den = np.where(np.abs(old) > 1e-12, np.abs(old), 1.0)
            worst.append((float(np.max(np.abs(new - old) / den)), c))
        worst.sort(reverse=True)
        print(f'{path}  최대 상대차 상위 3: ' +
              ', '.join(f'{c} {v:.3f}' for v, c in worst[:3]))
        if apply:
            for c in B.D3_NAMES:
                if c in df.columns:
                    df[c] = [r[c] for r in rows]
            df['n_conformers'] = [r['n_conformers'] for r in rows]
            df['conformer_w_max'] = [r['_w_max'] for r in rows]
            df.to_csv(path, index=False)
    print('[적용]' if apply else '[미적용] --check 였다')
    return 0


if __name__ == '__main__':
    sys.exit(main())
