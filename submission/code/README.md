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
| `descriptor_cache_xtb*.csv` | Structural (RDKit) and electronic (GFN2-xTB) descriptors |
| `panel_all.csv`, `panel_new.csv` | Compound panels and class labels |
| `inchikeys.csv` | Chemical identity resolution input |
| `deltaE/` | Measured 29-dye colorimetric array response, three compound sets, two imaging channels |
| `outputs/` | Permutation results, distance tables, parsed AEGL values, figure annotations |
| `paper_acs/` | Manuscript, Supporting Information, cover letter, figures, PDFs (ACS submission, canonical) |
| `archive/jhm_v1/` | Frozen parallel version for a different journal; not maintained |
| `sources/` | Primary sources: the parsed EPA AEGL document, the collected register and designation records, the array component list |

The layout is flat rather than split into `src/` and `data/` because the scripts
address their inputs and `outputs/` by relative path. A clone runs as it stands.

Five files in `outputs/` are shipped rather than regenerated, because
regenerating them needs either the `xtb` binary or a B3LYP/def2-SVPD run of
about an hour:

- `aegl_electronic.csv` and `xtb_fast_vs_slow.csv` (GFN2-xTB over the
  roughly 230 AEGL training compounds)
- `dft_omega.csv`, `dft_dipole.csv` and `dft_vs_gfn2.csv`
  (B3LYP/def2-SVPD over the twenty-eight-compound panel)

`run_analysis.sh` therefore calls neither `xtb_fast.py` nor `dft_omega.py`;
both commands are written out in its comments. Everything downstream of those
five files is regenerated from the sources in this repository, including the
parse of the EPA document.

The version of the `xtb` binary that produced `descriptor_cache_xtb*.csv` was
not recorded, so that axis is not version-pinned. `requirements.txt` says so.

One step in `run_analysis.sh` needs the network and is slow:
`build_tox_trainset.py` resolves about 3,900 compound names against the PubChem
REST service at four requests a second, so a cold clone spends roughly twenty
minutes there. It is the only step that leaves the machine. Nothing else in the
pipeline does.

## Reproducing

```sh
pip install -r requirements.txt
sh run_analysis.sh                      # analysis, in dependency order
sh build_paper_acs.sh --figs            # figures, then LaTeX
python3 verify_consistency.py           # manuscript against outputs
python3 test_verify_consistency.py      # does the verifier catch injected defects
python3 ai_style_scan.py                # prose style scan
```

Run from a clone with the five shipped files above in `outputs/` and the rest
of `outputs/` empty, this reproduces the reported values exactly: smallest
nominal `p` of 0.0250 on the twenty-eight-compound panel, 0.1917 on the primary
thirty-five-compound panel and 0.0386 on the panel with eight simulants, with
Benjamini-Hochberg minima of 0.1632, 0.6299 and 0.2122. These six numbers are
checked against `outputs/multiplicity_c.json` by `verify_consistency.py`, which
reads this file, so a stale promise here fails the build rather than misleading
a reproducer. An earlier version of this paragraph said "the two shipped files"
and gave four of the six wrongly.

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

The designation and Mantel tests are exhaustive rather than sampled:
C(17,3) = 680, C(24,3) = 2024 and C(24,8) = 735,471 designations and 6! = 720
or 4!^3 = 13,824 relabellings, enumerated in full. Five analyses are sampled
instead, and each is seeded and declared where it is used: the Mantel power
simulation, the null match expectation, the premise label-permutation null in
four of its six cells, the random-subspace control and the descriptor
resampling. The rest is deterministic given the seed `0xC0DE`.

## Verification

`verify_consistency.py` reads every numerical claim in `paper_acs/main.tex`,
`paper_acs/supporting.tex`, the cover letter and the figure annotations, and
checks each against the generator output it came from. `test_verify_consistency.py`
injects historical and synthetic defects into the manuscript or the Supporting
Information and confirms the verifier catches each one. Both print their
current counts; we do not repeat the counts here, because a number in a README
is one more thing that goes stale.

The descriptor reproduction gate has a stated scope and prints it on every run.
`build_descriptors.py --verify` recomputes the RDKit-derived columns of the
five caches the analysis reads and compares them to the deposited values; that
is 53 columns in each, out of 53 to 62 numeric columns depending on the cache.
The nine columns it does not recompute are the eight semiempirical electronic
coordinates and their derivatives, which need the `xtb` binary whose version
was not recorded, together with the count of conformers that binary used. The gate names those eight rather than counting them as
checked, fails if any other column falls outside that declared list, and fails
if a cache file is missing instead of skipping it. The eight are checked
separately by the B3LYP/def2-SVPD recomputation in the Supporting Information.

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

## Licence

Code (`*.py`, `*.sh`) is under the MIT licence; data, tables, manuscript
sources and this file are under CC BY 4.0. Two classes of file in this deposit
were produced elsewhere and are not ours to license: the reproduced EPA AEGL
document under `sources/`, and the colorimetric Delta E_00 tables under
`deltaE/`, which were measured and published in the three source studies cited
in the manuscript. `LICENSE` states this in full.

The state the manuscript reports is the tag `acs-submission-20261107`. Branches
move; cite the tag.

## Citation

Please cite the paper. For the code and data, cite this repository at the
tagged release `v1.0`, which is the state the manuscript refers to.
