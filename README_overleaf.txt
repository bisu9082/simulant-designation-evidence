Whether Chemical Warfare Agent Simulant Designations Look Supported
Depends on Which Compounds Are in the Comparison

Journal of Hazardous Materials 투고본 (2026-10-24, 공저자 1차 교정 반영)
저장소: https://github.com/bisu9082/simulant-designation-evidence


담긴 파일
---------
main.tex            본문 (elsarticle, 1p/times)
supplementary.tex   부속자료
cover_letter.tex    커버레터
highlights.tex      하이라이트
refs.bib            참고문헌 84건
Fig_main.png        본문 Figure 1
Fig_spaces.png      본문 Figure 2
Fig_toxicity.png    본문 Figure 3
Fig_h1.png          부속 Figure S3
Fig_plocal.png      부속 Figure S1
Fig_mde.png         부속 Figure S2


컴파일 순서 (중요)
------------------
main.tex 와 supplementary.tex 는 xr 패키지로 서로를 참조한다.
  main.tex          \externaldocument[S-]{supplementary}
  supplementary.tex \externaldocument[M-]{main}

상호참조이므로 한 번씩만 돌리면 양쪽 모두 ?? 가 남는다. 순서는 이렇다.

  1. supplementary.tex 컴파일   (M- 참조가 ?? 로 나온다. 정상)
  2. main.tex 컴파일            (S- 참조가 풀린다)
  3. supplementary.tex 재컴파일 (M- 참조가 풀린다)
  4. main.tex 재컴파일

오버리프에서는 Menu > Main document 를 바꿔 가며 각각 컴파일한다.
4회를 돌린 뒤 로그에 'There were undefined references' 가 남아 있으면
한 번 더 돌린다. 정상 산출물은 본문 25쪽, 부속 16쪽이다.


저자 (6인)
----------
Ku Kang / Jin Yoo          공동 제1저자        CBRN Defense Research Institute
Jeongyun Kim                                   CBRN Defense Research Institute
Miae Kim                                       CBRN Defense Research Institute
Doo-Hee Lee                공동교신            Korea Military Academy
Yoon Jeong Jang            공동교신            CBRN Defense Research Institute

1쪽 각주에는 교신저자 두 명의 주소만 찍힌다. elsarticle 이 모든 ead 를 한
각주로 모으기 때문에 비교신저자까지 실으면 전원이 교신인 것처럼 읽힌다.
비교신저자 네 명의 이메일과 전원의 ORCID 는 main.tex 머리의 주석 블록에
있으니 투고 시스템에 그대로 입력한다. Miae Kim 은 ORCID 없음(확인 완료).


별도 제출물
-----------
graphical_abstract.pdf 는 본문이 includegraphics 하지 않으므로 이 묶음에
없다. 투고 시스템에 Graphical Abstract 로 따로 올린다.


투고 전 남은 확인 사항
----------------------
- 커버레터 수신인이 'Dear Editors,' 이다. JHM 은 주제별 복수 Editor-in-Chief
  체제이므로 담당 에디터 확인 후 'Dear Dr. <Name>,' 로 바꾸거나 그대로 둔다.
- 문헌 모사체 5종 중 DCNP 만 출처가 붙어 있다. DFP, DCP, CEES,
  methyl salicylate 는 아직 인용 없이 소개된다. 공저자에게 출처를 받으면
  2.1 절에 같은 방식으로 추가한다.
- arsenic trichloride(루이사이트 모사체)는 패널에 넣지 않고 Limitations 에서
  '시험하지 않은 경로'로 인용했다. 패널에 넣으려면 전체 재실행이 필요하다.
