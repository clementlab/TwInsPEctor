# TwInsPEctor

TwInsPEctor is a tool for analyzing twin prime editing outcomes from next-generation sequencing reads, including data processing and visualization.

It utilizes CRISPResso2 for alignment of reads to a wiltype-edited compound reference.

Reads are categorized into one of the following eight allele types complete with detailed visualizations:
- Perfect TPE - complete programmed edit without indels.
- TPE Indel -  complete programmed edit with indels.
- Left Flap - at least N consecutive programmed bases starting from the left but not from the right.
- Right Flap - at least N consecutive programmed bases starting from the right but not from the left.
- Imperfect TPE - incomplete programmed edit (neither or both flaps).
- Imperfect WT - incomplete wildtype sequence and none of the programmed edit.
- WT Indel - complete wildtype sequence with indels and none of the programmed edit.
- WT - complete wildtype sequence without indels and none of the programmed edit.

## Features

- Command-line interface for streamlined analysis
- Data processing utilities for twin prime editing experiments
- Visualization tools for editing outcomes

## Installation

### Using Conda

```bash
conda install bioconda::twinspector
```

### From Source

```bash
git clone https://github.com/clementlab/TwInsPEctor.git
cd TwInsPEctor
pip install .
```

## Overview of command line interface

### Required arguments

| Argument | Description |
|--------|-------------|
| `-r1`, `--fastq_r1` | path to fastq r1 file. |
| `-w`, `--wt_seq` | full wildtype reference amplicon sequence including spacers. |
| `-t`, `--twin_seq` | full twin-pe reference amplicon sequence with 5′ and 3′ ends identical to the wildtype amplicon. |
| `-g`, `--peg_spacers` | comma-separated pegRNA spacer sequences: `<spacerA>,<spacerB>` |
---

### Optional Arguments

| Argument | Description | Default |
|--------|-------------|--------|
| `-r2`, `--fastq_r2` | path to fastq r2 file for paired-end data. | None |
| `-o`, `--output_root` | root output directory for TwInsPEctor results. If not provided, a folder is created in the working directory based on input fastq names. | auto |
| `-ne`, `--num_changes_to_check` | minimum number of programmed bases that must be edited for classification. | `2` |
| `-rcm`, `--recoding_mode` | enable recoding mode when edits consist only of base substitutions. | off |
| `-dmas`, `--default_min_aln_score` | minimum homology score for CRISPResso2 to align read to compound reference. | `50` |
| `-pfr`, `--plot_full_reads` | display full read sequences in allele tables. | off |
| `-ied`, `--ignore_extraspacer_deletions` | ignore deletions outside the edit window (beyond spacers). | off |
| `-nat`, `--no_allele_tables` | skip generation of allele tables to reduce runtime. | off |
| `-mfa`, `--min_frequency_alleles` | minimum percent read frequency required to report an allele. | `0.0` |
| `-mnr`, `--max_n_rows` | maximum number of allele rows displayed in tables. | `50` |
| `-nrr`, `--no_rerun` | do not rerun CRISPResso2 if the same parameters were already completed. | off |
| `-kco`, `--keep_crispresso_outputs` | preserve CRISPResso2 output folders after analysis. | off |

### Usage

```
TwInsPEctor -r1 <FASTQ_R1> [-r2 <FASTQ_R2>] -w <WT_SEQUENCE> -t <TWINPE_SEQUENCE> -g <PEG_SPACER_A>,<PEG_SPACER_B> [options]
```

After installation, use the CLI for help:

```bash
TwInsPEctor --help
```

Or run the main module directly:

```bash
python -m TwInsPEctor
```

## Requirements

- Python >=3.8
- [CRISPResso2](https://github.com/pinellolab/CRISPResso2) (installed via Bioconda)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Authors

- Nate Masson
- Kendell Clement

## Dependency Notice

This software requires [CRISPResso2](https://github.com/pinellolab/CRISPResso2) to be installed separately.

CRISPResso2 is distributed under its own license terms, which may
restrict commercial use.

Users are responsible for ensuring compliance with the CRISPResso2
license when using this software.

This project does not redistribute CRISPResso2 and does not grant
any rights to it.
