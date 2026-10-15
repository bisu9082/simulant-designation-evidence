#!/usr/bin/env python3
"""군용 화학작용제 CAS 명부 및 동정 검증 실패 목록 (2026-10-11 신설).

종전 제외 규칙은 문자열 매칭이었다:
    CWA_NAME = ('agent g', 'agent v', 'mustard', 'lewisite', 'sarin', 'soman', 'tabun')
AEGL 표의 'GF Agent' 는 'agent g' 에 걸리지 않아 신경작용제가 QSAR 학습세트에
그대로 들어갔다. 포스겐 옥심, 아다미사이트, 비소계 작용제도 마찬가지였다.
여기서는 CAS 로 명시한다.
"""

# AEGL Compiled Table 에서 군용 화학작용제로 식별되는 전 항목.
# 28종 패널의 real_agent 11종을 포함한다.
CWA = {
    '77-81-6':    'Agent GA (Tabun)',
    '107-44-8':   'Agent GB (Sarin)',
    '96-64-0':    'Agent GD (Soman)',
    '329-99-7':   'GF Agent (Cyclosarin)',
    '50782-69-9': 'Agent VX',
    '505-60-2':   'Sulfur Mustard (HD)',
    '51-75-2':    'Nitrogen Mustard-2 (HN-2)',
    '538-07-8':   'Nitrogen Mustard-1 (HN-1)',
    '555-77-1':   'Nitrogen Mustard-3 (HN-3)',
    '541-25-3':   'Lewisite 1 (L)',
    '593-89-5':   'Methyldichloroarsine (MD)',
    '712-48-1':   'Diphenylchloroarsine (DA)',
    '578-94-9':   'Adamsite (DM)',
    '1794-86-1':  'Phosgene oxime (CX)',
    '74-90-8':    'Hydrogen cyanide (AC)',
    '506-77-4':   'Cyanogen chloride (CK)',
    '75-44-5':    'Phosgene (CG)',
    '76-06-2':    'Chloropicrin (PS)',
    # 2026-10-11 추가. 파서 교정으로 'AEGL 3  0.52  0.17  0.086  NR  NR' 처럼
    # 긴 지속시간이 NR 인 행이 복구되면서 학습세트에 새로 편입된 군용작용제다.
    # 종전에는 행 자체가 폐기돼 보이지 않았다. 셋 모두 화학무기금지협약 등재
    # 작용제이므로 학습에서 제외하고 held-out 평가점으로만 쓴다.
    '6581-06-2':  'Agent BZ (3-quinuclidinyl benzilate)',
    '598-14-1':   'Ethyldichloroarsine (ED)',
    '696-28-6':   'Phenyldichloroarsine (PD)',
}

# 원표에 AEGL-3 가 있으나 단일 분자 구조가 없어 구조 기반 모델에 쓸 수 없는 항목.
# Ricin 은 단백질 독소이면서 협약 등재 작용제이므로 제외 사유가 중첩된다.
NO_SINGLE_STRUCTURE = {
    '9009-86-3':  'Ricin — 단백질 독소, 단일 분자 구조 없음 (협약 등재 작용제)',
    '1330-20-7':  'Xylenes — 이성질체 혼합물, 단일 구조 없음',
}

# validate_identity.py 가 적발한 동정 오류. 저장 구조가 이름과 다른 화합물이다.
#   19624-22-7 Pentaborane            -> EPN (유기인계 살충제) 구조로 해석됨
#   27137-85-5 Trichlorophenylsilane  -> PhSiCl3 + 2 HCl 다성분 레코드,
#                                        환산 분자량 284.47 (실제 211.55)
# 두 건 모두 신뢰할 구조를 자동으로 재지정할 수 없어 분석에서 제외한다.
#
# 2026-10-11. 이 목록을 손으로 유지했을 때 검증기 FAIL 집합과 갈라졌다
# (검증기는 107-07-3 을 FAIL, 목록은 27137-85-5 를 제외해 교집합이 1건이었다).
# 원인은 검증기 쪽 결함 두 건이었고 고쳤다. 재발을 막기 위해 이제 목록을
# 검증기 출력에서 읽고, 아래 문서화된 집합과 다르면 즉시 실패한다.
_DOCUMENTED = {
    '19624-22-7': 'Pentaborane — 저장 구조가 EPN',
    '27137-85-5': 'Trichlorophenylsilane — 다성분 레코드, 환산 분자량 오류',
}


def _load_identity_failed():
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'outputs', 'identity_validation.json')
    if not os.path.exists(path):
        raise RuntimeError(
            'outputs/identity_validation.json 없음. '
            'validate_identity.py 를 먼저 실행해야 한다.')
    d = json.load(open(path))
    found = {f['cas']: f.get('missing') for f in d.get('failures', [])}
    if set(found) != set(_DOCUMENTED):
        raise RuntimeError(
            '동정 검증 FAIL 집합이 문서화된 제외목록과 다르다.\n'
            f'  검증기: {sorted(found)}\n'
            f'  문서화: {sorted(_DOCUMENTED)}\n'
            '어느 쪽이 옳은지 확인한 뒤 양쪽을 맞춰야 한다.')
    return dict(_DOCUMENTED)


IDENTITY_FAILED = _load_identity_failed()


def is_cwa(cas):
    return str(cas).strip() in CWA


def identity_ok(cas):
    return str(cas).strip() not in IDENTITY_FAILED
