# simulant — 원논문 Methods 직독으로 드러난 사실 오류 3건 (2026-11-02)

Ku 가 올린 세 편의 자기 논문 PDF 를 직접 읽어 d1/d2/d3 의 출처를 확정하고,
원고의 사실 오류 세 건을 바로잡았다. 셋 다 "수행하지 않은 측정을 서술한"
유형이므로 [C-3] 위반에 해당한다. 커밋 d1a74d6.

## 데이터셋과 원논문의 대응 (확정)

d1 = Kim, Kang, Shin, Kim, Yoo, Kim, Lee, Lee, *J. Hazard. Mater.* **501**
(2026) 140703, DOI 10.1016/j.jhazmat.2025.140703 (bib: hazref).
작용제 GA, GB, GD, GF, VX + DMMP. 원자료 1,740행 = 29염료 × 6화합물 ×
5농도 × 2조명.

d2 = Kim, Yoo, Kang, Kim, Shin, Cha, Churchill, Kim, Lee, *Discover Sensors*
**1**(1) (2025) 29, DOI 10.1007/s44397-025-00032-y (bib: multiplexref, 신설).
작용제 AC, CK, CG, PS, HD, HN, L. Springer 페이지에서 저자·권·논문번호
실측했고 저자에 Kang 이 있어 자기인용 감사를 통과한다.

d3 = Shin, Kang, Kim, Kim, Yoo, Jang, Kim, Anslyn, Lee, *Sens. Actuators B*
**457** (2026) 139706, DOI 10.1016/j.snb.2026.139706 (bib: group3ref).
작용제 A-230, A-232, A-234, A-242.

## 염료 번호 — 세 편이 동일하다

Ku 는 "세 논문이 각각 다이 순번이 있다" 고 했으나, 세 논문의 Table 1 은
29종이 **같은 번호, 같은 순서**다. 1 Anthracene, 2 Pyrene, 3 Allura Red AC,
4 Quinine, 5 Rhodamine B, 6 Methyl Blue, 7 Fluorescein, 8 Methyl Orange,
9 Nile Red, 10 Safranin O, 11 Toluidine Blue O, 12 Nile Blue A, 13 Cresol
Red, 14 Eriochrome Black T, 15 Quinoline Yellow, 16 Bromophenol Blue,
17 Neutral Red, 18 Eosin Y, 19 Indigo Carmine, 20 Phenol Red, 21 Crystal
Violet, 22 L-Glutathione (환원형), 23 2,7-Dibromofluorene, 24 Ethyl
viologen dibromide, 25 4,4'-Dihydroxybenzophenone, 26 2,4-DNPH,
27 2,5-Diaminobenzene-1,4-dithiol 2염산염, 28 2-Hydrazinobenzothiazole,
29 Rhodamine 6G.

유일한 불일치는 **염료 15 Quinoline Yellow 의 용매**로, d1 과 d2 는 DMSO,
d3 는 DW 로 적었다. 셋 중 하나가 오기이며 Ku 확인이 필요하다.
[확인 필요: 염료 15 용매]

## 교정한 오류 3건

### 1. 농도 단위 — ppm 이 아니라 µM

세 원논문 모두 아세토니트릴 용액 **0, 10, 50, 100, 500 µM** 이다.
ACS 판본은 `10~ppm` / `500~ppm` 으로 적었고, JHM 보충자료는 "nominal
**vapour** concentrations ... (arbitrary units as reported)" 로 유보해
두었다. 기상이 아니라 용액이고 단위도 확인되었으므로 양쪽을 교정했다.

### 2. 반복측정은 존재하지 않는다

ACS 2.1 이 "Replicates of one compound correlate at +0.345 to +0.542" 로
적었다. `deltaE_premise.py` 의 within 항은 `ks = sorted(cen[a])`, 즉 같은
화합물의 서로 다른 두 **농도** 사이 코사인이다. n_within = 36 = 6화합물 ×
C(4,2)농도쌍 이 이를 확인한다. 원자료는 (light, dye, agent, conc) 조합당
한 행뿐이므로 반복측정이 아예 없다.

d2 만 "independent triplicate measurements (n=3) for the selected optimal
dye-analyte combinations" 을 SI Fig. S4 에 보고했고, 그 자료는 우리가 쓰는
피벗에 들어 있지 않다. d1·d3 은 반복측정 언급 자체가 없다.

따라서 이 검정이 세우는 명제는 "반복 일치가 물질간 일치를 넘는다" 가
아니라 "한 화합물이 농도를 가로질러 자기 자신과 일치하는 정도가 다른
화합물과 일치하는 정도를 넘는다" 다. 측정 반복성은 이 자료로 세울 수 없다.

### 3. d2 의 출처 오인용

JHM 보충자료가 AC/CG/CK/HD/HN/L/PS 세트의 출처를 `ondeviceref`
(Bae, Kang, Kim, Jang, Lee, YOLOv8, *Chemometrics Intell. Lab. Syst.*
**261** (2025) 105365) 로 적었으나, 그 논문은 이 데이터의 출처가 아니다.
실제 출처는 d2 = multiplexref 다.

## ACS Experimental Section 신설 (939단어 / 예산 1,100)

원논문에서 직접 가져온 어레이 측정 상세:
29종 100 mM 원액(EtOH/MeOH/DMSO/DW/THF/Toluene, 초음파 용해) → MeCN 에
1:1000 희석하여 작업농도 100 µM, 잔류 원액용매 0.1% (v/v) 미만 →
2.0 mL 염료용액에 MeCN 연속희석 분석물 투입 → 무배양 즉시 촬영 →
가시광 및 365 nm UV 이중조명, 암상자 내 스마트폰 카메라 고정,
22 ± 3 °C / 35 ± 5% RH → ROI 중앙 20 × 10 px, 5지점(중앙·상·하·좌·우)
RGB 평균 → CIE L*a*b* → 동일 염료·동일 채널의 0 µM 영상 대비 CIEDE2000.
작용제는 전량 화생방연구소(OPCW 지정시설) 제조·취급, GC-MS 순도는
d2 기준 루이사이트 67% 외 전부 95% 초과.

## 검증기 보강

- **사전점검 중복행 검사 신설.** 피벗에 (Agent, Concentration) 이 두 번
  나오면 반복측정이 실제로 존재한다는 뜻이므로 멈춘다. 중복 주입으로
  작동을 확인했다.
- **SUPERSEDED 금지항목 2종 신설.** ppm 재발과 반복측정 주장 재발.
  AEGL 의 정당한 ppm 은 면제 패턴으로 통과시킨다.
- **ACS 전용 주입시험 5건 신설.** 111 → 116건, 실패 0.
- **주입시험의 자기오염 결함 수정.** SIGTERM 으로 중단되면 주입값이 본문에
  남았다. 실제로 이번에 `paper/main.tex` 에 0.0221 자리의 0.0321 이 남아
  있었고 검증기가 그것을 잡았다. atexit 와 시그널 핸들러로 어느 경로로
  죽든 원본을 되돌리게 했다.
- conc_gb_lo/hi 앵커를 단위째 µM 으로 교정해 ppm 회귀가 다시 잡히게 했다.

## 상태

본문 29쪽, 보충 22쪽, ACS 18쪽 전부 오류 없이 빌드. 검증 402건 PASS,
주입시험 116건 실패 0. elsarticle 은 컨테이너 재생성으로 소실되어 CTAN 에서
TEXMFLOCAL 에 재설치했다.

## 남은 일

- 염료 15 용매 불일치 확인 (d1·d2 DMSO vs d3 DW)
- TOC 도판 재작도 (현재 JHM 프레임판)
- SI 를 ACS 형식으로 전환
- ACS Sensors 사전문의 발송 (초안 완료, 편집장 이름 필요)
- f+(P) DFT 재실행 (화합물 20에서 중단된 상태)
- 태그 v1.0 이 아직 aec3951 에 있다
- 발표용 Beamer 패키지가 저장소 밖에 있어 컨테이너 재생성으로 소실.
  `make_fig_talk_*.py` 는 저장소에 있으므로 재생성 가능.
