# MI 702 — Tn-seq homework

Analysis of a transposon insertion sequencing screen in Group A *Streptococcus*,
for MI 702 Applied Infection and Immunity.

## For students

**[▶ Open the homework in Google Colab](https://colab.research.google.com/github/nkorotkova-uky/mi702-tnseq/blob/main/notebooks/MI702_Tnseq_homework.ipynb)**

That link opens your own private copy. Run the cells in order, starting with the
setup cell. Nothing you do affects this repository or any other student's work.

If you break something, close the tab without saving and open the link again.

## What this is

A saturated *S. pyogenes* M1T1 5448 Krmit (mariner) transposon library was exposed
to human group IIA secreted phospholipase A₂ (hGIIA), a bactericidal host enzyme.
Four treated and four untreated replicates were sequenced. The homework asks which
genes change susceptibility to the enzyme — and how far the data supports the
answer.

| | |
|---|---|
| Published study | van Hensbergen et al. 2018, *PLoS Pathog* 14(10):e1007348 ([PMC6201954](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6201954/)) |
| Reads | ENA [PRJEB27626](https://www.ebi.ac.uk/ena/browser/view/PRJEB27626), run ERR2683524 |
| Reference genome | GenBank CP008776 / RefSeq NZ_CP008776 |
| Library | Le Breton et al. 2015, *Sci Rep* ([PMC4440532](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4440532/)) |

## Repository layout

```
notebooks/MI702_Tnseq_homework.ipynb   the assignment
mi702tnseq/                            analysis package the notebook imports
  pipeline.py    demultiplex, trim, map insertions, count per gene
  stats.py       fold change, BH correction, permutation and t-tests
  plots.py       volcano plot, insertion tracks
  data.py        release-asset download and loaders
pyproject.toml
```

The analysis code lives in the package rather than the notebook on purpose: the
notebook holds the questions and the student's answers, so an accidental edit
cannot break the pipeline, and `pip install` pins a tagged version.

## Data

Datasets are **release assets**, not repository files.

| Asset | Size | Contents |
|---|---|---|
| `mi702-tnseq-core.zip` | 3.9 MB | genome, annotation, TA-site list, barcodes, per-gene tables, insertion counts |
| `mi702-tnseq-reads.zip` | 156 MB | multiplexed, demultiplexed and trimmed FASTQ |

The notebook fetches `core` in the setup cell and `reads` only when Part B needs
them, so most of the assignment never downloads the large file.

**The instructor answer key is deliberately not in this repository or in either
release asset.**

## Dataset caveats

The reads are a **300,000-per-sample subset** of a run that originally had about
31.6 million reads per sample. The subset reproduces the published effect sizes
well (Spearman ρ ≈ 0.75 genome-wide, ≈ 0.88 across the published hits) but cannot
reproduce their statistical significance. One published hit, *fabT*, has no reads
at this depth. These limitations are intentional and are the subject of several
homework questions.

Per-sample insertion-site saturation in the subset (3.8–6.1%) is a property **of
the subset**, not of the library, which was far denser.

## License and citation

The published per-gene tables are CC-BY from PLoS Pathogens; the reads are public
in ENA. Please retain the citation to van Hensbergen et al. 2018 in any derived
teaching material.
