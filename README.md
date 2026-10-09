# ie_finder_site

A reduced, self-contained copy of the MGE_finder discovery workflow. The integrase HMMs, the search scripts, and the conda environment are in this repository. `UPSTREAM` records which commit of [MGE_finder](https://github.com/rljech13/MGE_finder) the copy was taken from. A checkout of that repository is not required.

The pipeline reads a genome assembly in FASTA and writes attL/attR, the integrase, and the element span in the coordinates of that assembly.

## What an element is

A confident element is a phage-type tyrosine integrase next to an opposite-strand tRNA, with a direct repeat of that tRNA 3′ end. The integrase is longer than 300 aa, attL does not overlap a CDS, and the repeat is either an exact run of at least 14 bp or a BLAST hit of at least 17 bp. Alignment gaps in the repeat are allowed. Setting `attl_reject_gapped: true` rejects them, as the published finder did. Among the BLAST hits anchored at the tRNA 3′ end, attL is the one with the highest bitscore, so a long but poor match far from the tRNA does not win over a near-perfect repeat. Setting `attl_select_by: length` takes the longest hit instead, as the published finder did. Those thresholds are in `search_params.yaml`. Cohort deduplication from the main repository is not part of this copy.

## Requirements

Linux. The first `./run.sh` builds the conda environment named in `envs/IE_finder.yaml` (`IE_finder_site`): Snakemake, Python 3.10, Biopython, pandas, PyYAML, Prodigal, HMMER, Aragorn, and BLAST+. Later runs reuse that environment. The first build needs a network connection.

If `conda` is already installed, the script uses it. Otherwise it installs Miniforge into `~/miniforge3` (`MINIFORGE_ROOT` changes that directory) and creates the environment there.

## Install

```bash
git clone https://github.com/rljech13/ie_finder_site.git
cd ie_finder_site
```

## Run

One assembly, or a directory of assemblies. The sample name is the file name without its last suffix, so name the file after the strain.

```bash
./run.sh TTHB27c.fasta outdir
./run.sh genomes_dir outdir
```

Accepted suffixes in a directory are `.fna`, `.fa`, `.fasta`, and `.fsa`. Two files that share a stem are an error.

`outdir/` receives three files per assembly:

| file | role |
|---|---|
| `TTHB27c.ie.gff3` | features on the input contigs, for a browser track or a download |
| `TTHB27c.ie.gbk` | the same contigs and features in GenBank |
| `TTHB27c.ie.report.txt` | one report: counts, why a candidate was rejected, the published coordinates, and any step log that was not empty |

Coordinates are 1-based and inclusive. They refer to the contigs of the FASTA you passed in. The first column of the GFF3 is that FASTA record id.

Each element has four features:

| GFF3 type | GenBank key | meaning |
|---|---|---|
| `mobile_genetic_element` | `mobile_element` | span from attR to attL |
| `attachment_site` with `Note=attR` | `misc_feature` / `attR` | the tRNA that was the integration site |
| `attachment_site` with `Note=attL` | `misc_feature` / `attL` | the direct repeat of the tRNA 3′ end |
| `CDS` | `CDS` | the tyrosine integrase |

`confidence=confident` passed the filter. A run with no confident element still writes a valid GFF3 header and an empty GenBank, plus a report that says so.

Intermediate tables, BLAST output, and the cut-out island sequences are written to a temporary directory and removed when the run succeeds. Do not publish those island GenBank files: their coordinates are local to the extracted sequence, which may have been reverse-complemented. Only `*.ie.gff3` and `*.ie.gbk` use the assembly coordinates.

## Options

```bash
ANNOTATE_ALL=1 ./run.sh TTHB27c.fasta outdir
KEEP_WORK=1 ./run.sh TTHB27c.fasta outdir
SEARCH_PARAMS=strict.yaml ./run.sh TTHB27c.fasta outdir
./run.sh TTHB27c.fasta outdir --cores 8
```

No activation step. `./run.sh` puts `IE_finder_site` on `PATH` for that run. After a change to `envs/IE_finder.yaml`, delete the environment (`conda env remove -n IE_finder_site`) and run again. The script rebuilds it.

`ANNOTATE_ALL=1` also writes candidates that have an attL coordinate but failed a later filter. In the GFF3 and the report their `confidence` is `candidate`.

`KEEP_WORK=1` leaves the temporary directory. The path is printed only when a run fails or when this flag is set. Use it to see the raw tables behind the report.

`SEARCH_PARAMS=strict.yaml` reads search thresholds from that file on top of `search_params.yaml`. Keys set in it win, and keys left out keep their values from `search_params.yaml`. An unknown key or a bad value stops the run before any step starts. For example, a file with these two lines brings back the published attL rules:

```yaml
attl_select_by: length
attl_reject_gapped: true
```

Arguments after the output directory are passed to Snakemake.

## The report

`strain.ie.report.txt` is the log for that genome. It has five blocks:

- `search_params` — the thresholds this run used, after any `SEARCH_PARAMS` override
- `counts` — integrase HMM hits, integrase–tRNA pairs, attL BLAST hits, how many passed
- `candidates` — one line each: `confident` or `rejected` plus the first failing reason (`exact_run_too_short`, a CDS overlap, an integrase shorter than 300 aa, a gapped alignment)
- `published_coordinates` and `audit` — the full tables, so a rejected locus can be reconstructed without the temporary files
- `logs` — stdout of the steps that printed something. Empty logs are omitted

If the run crashes, the report is still written from whichever steps finished.
