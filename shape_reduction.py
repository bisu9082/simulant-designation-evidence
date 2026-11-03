#!/usr/bin/env python3
r"""관성텐서 파생 좌표 축소 민감도 (2026-11-03 신설, 작업 A-1).

왜 필요한가. GATE 8 8차에서 리뷰어 2 가 집은 것이고, 내가 확인했다.
보유 39좌표 중 9개가 하나의 관성텐서 3-벡터에서 나온다:

  PMI1 PMI2 (PMI3) NPR1 NPR2 Asphericity Eccentricity
  InertialShapeFactor SpherocityIndex RadiusOfGyration

이 아홉은 세 개의 주관성모멘트의 함수다. 그런데 z 표준화는 좌표마다 같은
기여를 주므로, 실효 자유도 3 인 블록이 전체 378쌍 평균 제곱거리의 23.1% 를
가져간다. 독립 좌표 3개인 전자공간은 7.7% 다. 셋 대 셋인데 세 배다.

개별 쌍에서는 더 심하다. HCN 과 시안화염소는 준선형이라 NPR1=0, NPR2=1,
Asphericity=1, Eccentricity=1 이 동시에 경계값에 붙는다. 그 결과 형태
블록이 HCN↔아세토니트릴 거리 제곱의 69.8%, CK↔아세토니트릴의 61.9% 를
차지한다. 아세토니트릴은 match_rule_audit 이 그 두 작용제의 blocker 로
기록한 화합물이다. 11종 중 2종의 매칭 차단이 준선형 분자의 형태 좌표
쌓임으로 결정된다.

논문은 이미 같은 이유로 omega 를 모든 거리에서 뺐다 — "동반 좌표의
결정론적 함수이고 선형 상관필터는 그것을 못 본다"(spaces.OMEGA_DERIVED).
같은 근거가 이 블록에 그대로 적용된다. 상관필터가 못 잡는 이유도 같다:
관계가 비(ratio)이므로 선형 상관이 0.95 를 넘지 않는다(최대 0.9364).

한 가지는 리뷰어가 틀렸다. 리뷰어는 항등식이 기계정밀도로 성립한다고
적었는데, 예치 캐시에서 NPR1xPMI3 대 PMI1 의 최대 상대오차는 6.1%,
InertialShapeFactor 는 15.3% 다. 서술자마다 독립으로 볼츠만 평균하기
때문이다. 즉 중복은 정확한 함수관계가 아니라 1~15% 여유가 있는 근사이고,
그 여유는 화학이 아니라 평균 방식의 산물이다. 그래서 축소는 더 정당하다 —
남겨 둘 이유가 되는 '독립 정보' 가 그 여유뿐이기 때문이다.

이 스크립트는 정본을 덮어쓰지 않는다. 두 공간을 나란히 돌려 예치한다.
어느 쪽을 1차로 할지는 결과를 보고 정한다(논문 자신의 방법론이 요구하는
절차다). 산출: outputs/shape_reduction.json
"""
import itertools
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from spaces import ELECTRONIC_28, META_NUM, OMEGA_DERIVED, PROTECTED, pick

ELEC3 = [c for c in ELECTRONIC_28 if c != OMEGA_DERIVED]

# 관성텐서에서 파생되는 좌표 전부. 아래 PRIMITIVE 만 남기고 나머지는
# 거리에서 뺀다. omega 와 같은 취급이다.
SHAPE_ALL = ['PMI1', 'PMI2', 'PMI3', 'NPR1', 'NPR2', 'Asphericity',
             'Eccentricity', 'InertialShapeFactor', 'SpherocityIndex',
             'RadiusOfGyration']
# 원시량. omega 선례의 규칙("파생 합성량을 버리고 원시량을 남긴다")을 따른다.
SHAPE_PRIMITIVE = ['PMI1', 'PMI2', 'PMI3']
SHAPE_DERIVED = [c for c in SHAPE_ALL if c not in SHAPE_PRIMITIVE]

PANELS = [('panel28_incumbent3', 'descriptor_cache_xtb.csv', ['simulant']),
          ('panel35_incumbent3', 'descriptor_cache_xtb_35.csv', ['simulant']),
          ('panel35_incumbent_plus_lit8', 'descriptor_cache_xtb_35.csv',
           ['simulant', 'literature_simulant'])]
STATS = ['mean_diff', 'median_diff', 'mean_logratio', 'mean_dsim',
         'n_closer', 'n_valid']


def candidates(df):
    num = [c for c in df.select_dtypes('number').columns
           if df[c].notna().sum() == len(df) and c not in META_NUM]
    return [c for c in num if c not in
            ('q_P', 'q_phosphoryl_O', 'wiberg_PO', 'fukui_f_plus_P',
             'gap_eV')]


def zmat(df, keep):
    X = df[keep].to_numpy(float)
    return (X - X.mean(0)) / X.std(0, ddof=1)


def dist(df, keep):
    Z = zmat(df, keep)
    return np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))


def identity_audit(df):
    """항등식이 캐시에서 얼마나 깨져 있는가. 리뷰어 주장 대 실측."""
    out = {}
    for name, a, b in [
            ('NPR1*PMI3 vs PMI1', df.NPR1 * df.PMI3, df.PMI1),
            ('NPR2*PMI3 vs PMI2', df.NPR2 * df.PMI3, df.PMI2),
            ('PMI2/(PMI1*PMI3) vs ISF',
             df.PMI2 / (df.PMI1 * df.PMI3), df.InertialShapeFactor)]:
        m = np.isfinite(a) & np.isfinite(b) & (b.abs() > 1e-12)
        out[name] = float(((a[m] - b[m]).abs() / b[m].abs()).max())
    out['note'] = ('볼츠만 평균을 서술자마다 독립으로 적용하므로 항등식이 '
                   '깨진다. 기계정밀도가 아니다.')
    return out


def block_share(df, keep):
    """형태 블록과 전자 블록이 평균 제곱거리의 몇 %를 가져가는가."""
    Z = zmat(df, keep)
    sq = (Z[:, None, :] - Z[None, :, :]) ** 2
    tot = sq.sum(-1)
    iu = np.triu_indices(len(df), 1)
    def share(cols):
        i = [keep.index(c) for c in cols if c in keep]
        if not i:
            return 0.0
        return float(100 * sq[:, :, i].sum(-1)[iu].sum() / tot[iu].sum())
    nm = list(df.compound_name)
    pair = {}
    for a, b in [('hydrogen cyanide (AC)', 'acetonitrile'),
                 ('cyanogen chloride (CK)', 'acetonitrile')]:
        if a in nm and b in nm:
            i, j = nm.index(a), nm.index(b)
            ii = [keep.index(c) for c in SHAPE_ALL if c in keep]
            pair[f'{a} | {b}'] = float(
                100 * sq[i, j, ii].sum() / tot[i, j]) if ii else 0.0
    return {'shape_pct': share(SHAPE_ALL), 'electronic_pct': share(ELEC3),
            'n_shape_kept': len([c for c in keep if c in SHAPE_ALL]),
            'n_elec_kept': len([c for c in keep if c in ELEC3]),
            'near_linear_pairs_shape_pct': pair}


def exact_perm(D, kl, sim_classes):
    """rollback_c 와 같은 귀무·같은 통계량."""
    D = D.copy().astype(float)
    np.fill_diagonal(D, np.inf)
    nn_med = float(np.median(D.min(1)))
    i_real = [i for i, k in enumerate(kl) if k == 'real_agent']
    i_pool = [i for i, k in enumerate(kl) if k != 'real_agent']
    i_sim = [i for i, k in enumerate(kl) if k in sim_classes]
    pos = {g: j for j, g in enumerate(i_pool)}
    true_combo = tuple(sorted(pos[s] for s in i_sim))
    ALL = list(itertools.combinations(range(len(i_pool)), len(true_combo)))
    Ds = D[np.ix_(i_real, i_pool)]
    C = np.asarray(ALL, dtype=np.int32)
    null = {s: np.empty(len(ALL)) for s in STATS}
    npool = len(i_pool)
    for a0 in range(0, len(C), 8000):
        cb = C[a0:a0 + 8000]
        b = len(cb)
        m = np.zeros((b, npool), bool)
        m[np.arange(b)[:, None], cb] = True
        A3 = Ds[:, None, :]
        ds = np.where(m[None, :, :], A3, np.inf).min(2)
        di = np.where(~m[None, :, :], A3, np.inf).min(2)
        null['mean_diff'][a0:a0 + b] = (ds - di).mean(0)
        null['median_diff'][a0:a0 + b] = np.median(ds - di, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            lr = np.log(ds / di)
        fin = np.isfinite(lr)
        null['mean_logratio'][a0:a0 + b] = np.where(
            fin.any(0), np.nansum(np.where(fin, lr, 0.0), 0)
            / np.maximum(fin.sum(0), 1), np.nan)
        null['mean_dsim'][a0:a0 + b] = ds.mean(0)
        null['n_closer'][a0:a0 + b] = -(ds < di).sum(0)
        null['n_valid'][a0:a0 + b] = -((ds / nn_med <= 2.0) & (ds <= di)).sum(0)
    ti = ALL.index(true_combo)
    return {s: dict(observed=float(null[s][ti]),
                    p=float((null[s] <= null[s][ti]).sum() / len(ALL)),
                    rank=int((null[s] <= null[s][ti]).sum()))
            for s in STATS}, len(ALL)


def mantel_family(df, keep):
    """12셀 Mantel 전수열거 + BH. sparse_indicators 와 같은 절차."""
    import deltaE_premise as m
    idx = {nm: i for i, nm in enumerate(df.compound_name)}
    st = [c for c in keep if c not in ELECTRONIC_28]
    el = [c for c in keep if c in ELEC3]
    SP = {'structural': st, 'electronic3': el, 'merged': st + el}
    DIST = {k: dist(df, v) for k, v in SP.items()}
    ps, cells = [], []
    for tag in ['d1', 'd2', 'd3']:
        for mode in ['uv', 'normal']:
            try:
                D, dy = m.load(tag, mode)
            except Exception:
                continue
            P = m.consensus(m.centred(m.raw_patterns(D, dy)))
            ags = [a for a in sorted(P) if a in m.MAP and m.MAP[a] in idx]
            if len(ags) < 4:
                continue
            k = len(ags)
            pr = list(itertools.combinations(range(k), 2))
            dis = np.array([1 - float(np.dot(P[ags[i]], P[ags[j]]))
                            for i, j in pr])
            perms = list(itertools.permutations(range(k)))
            for sp in SP:
                M = np.array([[DIST[sp][idx[m.MAP[a]], idx[m.MAP[b]]]
                               for b in ags] for a in ags])
                d = np.array([M[i, j] for i, j in pr])
                r = float(stats.spearmanr(d, dis).statistic)
                cnt = 0
                for pm in perms:
                    dp = np.array([M[pm[i], pm[j]] for i, j in pr])
                    rp = stats.spearmanr(dp, dis).statistic
                    if rp == rp and abs(rp) >= abs(r):
                        cnt += 1
                p = cnt / len(perms)
                ps.append(p)
                cells.append({'cell': f'{tag}_{mode}', 'space': sp,
                              'spearman': r, 'p_mantel': p})
    o = np.sort(np.array(ps))
    K = len(o)
    bh = np.minimum.accumulate((o * K / np.arange(1, K + 1))[::-1])[::-1]
    return {'n_cells': int(K), 'min_p': float(o[0]), 'min_q': float(bh[0]),
            'n_sig_q': int((bh < 0.05).sum()),
            'max_abs_rho': float(max(abs(c['spearman']) for c in cells)),
            'cells': cells}


def run(variant, drop):
    res = {}
    for label, cache, sim_cls in PANELS:
        df = pd.read_csv(cache)
        cand = [c for c in candidates(df) if c not in set(drop)]
        keep, dropped = pick(df, cand, PROTECTED)
        st = [c for c in keep if c not in ELECTRONIC_28]
        el = [c for c in keep if c in ELEC3]
        SP = {'structural': st, 'electronic3': el, 'merged': st + el}
        panel = {'n_dims': {k: len(v) for k, v in SP.items()},
                 'kept_shape': [c for c in keep if c in SHAPE_ALL],
                 'block_share_merged': block_share(df, st + el)}
        for sp, cols in SP.items():
            r, nperm = exact_perm(dist(df, cols), list(df.compound_class),
                                  sim_cls)
            panel[sp] = {'n_permutations': nperm, 'statistics': r}
        if label == 'panel28_incumbent3':
            panel['mantel'] = mantel_family(df, keep)
            panel['identity_audit'] = identity_audit(df)
        res[label] = panel
    return res


def main():
    print('=' * 78)
    print('[A-1] 관성텐서 파생 좌표 축소 — 정본을 덮지 않고 나란히 예치')
    print('=' * 78)
    out = {'generated_by': 'shape_reduction.py',
           'shape_all': SHAPE_ALL, 'shape_primitive': SHAPE_PRIMITIVE,
           'shape_derived_dropped': SHAPE_DERIVED,
           'rationale': ('omega 선례와 같다. 아홉 좌표가 관성텐서 3-벡터의 '
                         '함수이고 관계가 비(ratio)이므로 선형 상관필터가 '
                         '못 잡는다. 원시량 PMI1/2/3 만 남긴다.'),
           'variants': {}}
    for variant, drop in [('full', []), ('reduced', SHAPE_DERIVED)]:
        print(f'\n--- {variant} ---')
        r = run(variant, drop)
        out['variants'][variant] = r
        for label, p in r.items():
            md = p['merged']['statistics']['mean_diff']
            nv = p['merged']['statistics']['n_valid']
            bs = p['block_share_merged']
            print(f"  {label:<30} dims {p['n_dims']}  "
                  f"형태보유 {bs['n_shape_kept']}")
            print(f"     병합 mean_diff p={md['p']:.4f}  "
                  f"유효매칭 {int(-nv['observed'])} (p={nv['p']:.4f})  "
                  f"형태기여 {bs['shape_pct']:.1f}% 대 전자 "
                  f"{bs['electronic_pct']:.1f}%")
            if 'mantel' in p:
                mt = p['mantel']
                print(f"     Mantel {mt['n_cells']}셀 최소 p={mt['min_p']:.4f} "
                      f"최소 q={mt['min_q']:.4f} "
                      f"유의 {mt['n_sig_q']}  |rho|max={mt['max_abs_rho']:.3f}")
    os.makedirs('outputs', exist_ok=True)
    with open('outputs/shape_reduction.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print('\n-> outputs/shape_reduction.json')


if __name__ == '__main__':
    main()
