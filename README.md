# Simulant designation evidence

Analysis code, descriptor tables and permutation results for:

> **Whether Chemical Warfare Agent Simulant Designations Look Supported Depends
> on Which Compounds Are in the Comparison**
> Ku Kang, Yoon Jeong Jang.
> CBRN Defense Research Institute, Ministry of National Defense, Republic of Korea.
> Submitted to *Journal of Hazardous Materials*.

The paper asks whether three compounds in routine use as chemical warfare agent
simulants (dimethyl methylphosphonate, dichloromethane, dipropylene glycol
monomethyl ether) sit measurably closer to eleven authentic agents than
compounds drawn at random, under an exact permutation test over every possible
designation. The answer depends on which compounds are in the comparison, and
that dependence is the result.

## What is here

| Path | Contents |
|---|---|
| `*.py`, `*.sh` | Analysis, figure and verification scripts |
| `sources/` | EPA compiled AEGL document, the primary source for every toxicity value used here |
| `descriptor_cache_xtb*.csv` | Structural (RDKit) and electronic (GFN2-xTB) descriptors |
| `panel_all.csv`, `panel_new.csv` | Compound panels and class labels |
| `inchikeys.csv` | Chemical identity resolution input |
| `deltaE/` | Measured 29-dye colorimetric array response, three compound sets, two imaging channels |
| `outputs/` | Permutation results, distance tables, parsed AEGL values, figure annotations |
| `paper/` | Manuscript, Supplementary Material, cover letter, figures, PDFs |

The layout is flat rather than split into `src/` and `data/` because the scripts
address their inputs and `outputs/` by relative path. A clone runs as it stands.

Two files in `outputs/` are shipped rather than regenerated, because
regenerating them requires the `xtb` binary to recompute GFN2-xTB descriptors
for the roughly 230 AEGL training compounds:

- `aegl_electronic.csv`
- `xtb_fast_vs_slow.csv`

`run_analysis.sh` therefore does not call `xtb_fast.py`. Everything downstream
of those two files is regenerated from the sources in this repository,
including the parse of the EPA document.

## Reproducing

```sh
pip install -r requirements.txt
sh run_analysis.sh                      # analysis, in dependency order
sh build_paper.sh --figs                # figures, then LaTeX
python3 verify_consistency.py           # manuscript against outputs
python3 test_verify_consistency.py      # does the verifier catch injected defects
python3 ai_style_scan.py                # prose style scan
```

Run from a clone with the two shipped files above in `outputs/` and the rest of
`outputs/` empty, this reproduces the reported values exactly: smallest nominal
`p` of 0.0250 on the twenty-eight-compound panel, 0.1512 on the primary
thirty-five-compound panel, 0.0393 on the panel with eight simulants, and
Benjamini-Hochberg minima of 0.159, 0.617 and 0.212.

`run_analysis.sh` fixes the dependency order. Two constraints in it are not
obvious and were found by running a cold clone rather than by reading the code:

- `validate_identity.py` must precede `qsar_tox.py`, because `cwa_registry`
  reads the excluded records from `outputs/identity_validation.json`.
- `validate_identity.py` exits non-zero when it finds records whose stored
  structure belongs to a different compound. It currently reports two, which
  are the two the manuscript says were excluded. The pipeline treats that exit
  status as a report rather than a failure.

Running the scripts in another order leaves stale downstream outputs, which is
the failure mode this project hit repeatedly before the order was written down.

The exact permutation tests are exhaustive rather than sampled:
C(17,3) = 680, C(24,3) = 2024 and C(24,8) = 735,471 designations, enumerated in
full. The remaining analyses are deterministic given the seed `0xC0DE`.

## Verification

`verify_consistency.py` reads every numerical claim in `paper/main.tex`,
`paper/supplementary.tex`, the cover letter, the highlights and the figure
annotations, and checks each against the generator output it came from. It
currently performs 137 context checks, all passing. `test_verify_consistency.py`
injects 34 historical and synthetic defects into the manuscript or the
Supplementary Material and confirms the verifier catches all 34.

Both are part of the deliverable. The manuscript reports four audits that each
changed its numbers, and these scripts are what closed the propagation paths
those audits exposed.

## Data provenance

Structural descriptors are computed with RDKit from canonical SMILES. Electronic
descriptors are computed with GFN2-xTB in the gas phase. Chemical identities are
resolved against PubChem and against the identifiers printed in the regulatory
source documents; no property value is taken from a prediction service.

Acute inhalation values are AEGL-3 figures parsed from EPA primary sources, with
the extraction validated over the whole table by a second reader that shares
neither the text-extraction engine nor the algorithm.

Measured colorimetric array response is reported from earlier work by this
laboratory, cited in the manuscript.

`sources/compiled_aegls_27jul2018.pdf` is a United States Environmental
Protection Agency publication and is included so that the parse, and the
parsing errors the manuscript reports finding and correcting, can be checked
against the document itself.

## Dual-use

The optimisation direction is fixed throughout and stated in the manuscript:
detection-relevant fidelity is maximised while phosphorus electrophilicity,
acetylcholinesterase reactivity, alkylating capacity and predicted mammalian
toxicity are minimised. No analysis here ranks compounds by potency. The
structural filter excludes candidates and never orders them.

This repository contains computed molecular descriptors for known compounds,
along with regulator-published toxicity values. It contains no synthesis route,
no procedure, and no predictive model of toxicity for structures that do not
already exist. The toxicity model included here is deliberately not extended to
the organophosphorus class; the manuscript explains why.

## Citation

Please cite the paper. For the code and data, cite this repository at the
tagged release `v1.0`, which is the state the manuscript refers to.
