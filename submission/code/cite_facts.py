#!/usr/bin/env python3
"""외부 문헌에서 인용한 수치의 단일 출처 (2026-10-20 신설).

무엇이고 무엇이 아닌가. 이 스크립트는 계산하지 않는다. 본문이 다른 논문에서
가져다 쓴 수치를, 출처 DOI 와 조회일과 함께 한 곳에 적어 두는 것이다.
그러면 본문이 그 수치와 어긋날 때 verify_consistency.py 가 잡는다.

왜 필요한가. 표 3 사고의 원인은 계산 오류가 아니라 '본문 수치가 어떤
산출물과도 대조되지 않는 상태'였다. 외부 문헌에서 옮겨 적은 수치는 바로 그
상태에 있다. 아래 값은 2026-10-20 Semantic Scholar API 로 받은 각 논문의
출판사 초록에서 읽은 것이며, 본문 문장을 고칠 때는 원문을 다시 확인해야 한다.

Lavoie 2011 은 출판사가 초록을 비공개 처리했으므로 수치를 싣지 않는다.
본문도 그 논문에서 수치를 인용하지 않는다.
"""
import json

RETRIEVED = '2026-10-20'
SOURCE = 'Semantic Scholar Graph API, publisher abstract field'

FACTS = {
    'james2018': {
        'doi': '10.1002/jat.3527',
        'journal': 'Journal of Applied Toxicology',
        'volume': 38, 'number': 1, 'pages': '113-121', 'year': 2018,
        'n_simulants_identified': 22,
        'n_unsuitable': 17,
        'n_remaining': 5,
        'quote_basis': ('Twenty-two simulants were identified, of which 17 were '
                        'determined unsuitable for use in HVTs. The remaining '
                        'simulants (n = 5) were further scrutinized.'),
    },
    'james2021': {
        'doi': '10.3390/ijerph18168681',
        'journal': 'International Journal of Environmental Research and Public Health',
        'volume': 18, 'number': 16, 'pages': '8681', 'year': 2021,
        'n_chemicals_identified': 171,
        'n_discounted': 78,
        'n_suitable': 39,
        'n_possible': 54,
        'quote_basis': ('Of the 171 unique chemicals identified, 78 were '
                        'discounted for the risk they could pose to human '
                        'volunteers, 39 were deemed suitable for use, and a '
                        'further 54 were considered to be possible simulants.'),
    },
    'lavoie2011': {
        'doi': '10.1016/j.jhazmat.2011.07.077',
        'journal': 'Journal of Hazardous Materials',
        'volume': 194, 'pages': '85-91', 'year': 2011,
        'abstract_status': '[수집 불가: 출판사가 초록 비공개 처리]',
        'note': ('본문은 이 논문에서 수치를 인용하지 않는다. 방법(제약 유사성 '
                 '검색의 모사체 식별 적용)과 목적(가깝지만 덜 독성인 화합물) '
                 '진술만 쓰며, 그 범위는 제목과 S2 TLDR 로 확인했다.'),
    },
}


def main():
    out = {'retrieved': RETRIEVED, 'source': SOURCE, 'papers': FACTS}
    with open('outputs/cite_facts.json', 'w') as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(f'[cite_facts] outputs/cite_facts.json 기록 ({len(FACTS)}편, '
          f'조회 {RETRIEVED})')


if __name__ == '__main__':
    main()
