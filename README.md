# TwInsPEctor

TwInsPEctor is a tool for analyzing twin prime editing outcomes from next-generation sequencing reads, including data processing and visualization.

It utilizes CRISPResso2 for alignment of reads to wild type, edited, and wildtype-edited composite references.

Reads are classified into one of ten allele outcome categories:
- Perfect TPE - exact programmed edit.
- Dual Flap - meets both Flap A and Flap B requirements.
- Flap A - at least N contiguous bases from the start of only the pegRNA-a templated sequence.
- Flap B - at least N contiguous bases from the start of only the pegRNA-b templated sequence.
- Flap A Hybrid - meets Flap A requirements but retains wild type sequence from the opposing 5' flap.
- Flap B Hybrid - meets Flap B requirements but retains wild type sequence from the opposing 5' flap.
- Imperfect TPE - contains a partially edited sequence meeting none of the flap category requirements.
- Null - lacks any detectable wild type or edited sequence in between the two pegRNA nick sites.
- Imperfect WT - partial wildtype sequence with none of the edited sequence.
- WT - unedited wild type sequence.

## Features

- Command-line interface for streamlined analysis
- Data processing utilities for twin prime editing experiments
- Visualization tools for editing outcomes

## Installation

*Typical installation time on a normal desktop computer is less than 5 minutes.*

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

### Required Arguments

| Argument | Description |
| :--- | :--- |
| `-r1`, `--fastq_r1` | Path to FASTQ R1 file. |
| `-w`, `--wt_seq` | Full wild-type reference amplicon sequence including spacers. |
| `-t`, `--tpe_seq` | Full Twin prime edited reference amplicon sequence with 5' & 3' ends identical to wildtype reference amplicon. |
| `-g`, `--peg_spacers` | Comma-separated pegRNA spacer sequences: `<spacer A>,<spacer B>`. Should include bases immediately adjacent to but not including the PAM sequence (usually 20nt 5' of NGG). |

### Recommended Arguments

| Argument | Description |
| :--- | :--- |
| `-rt`, `--rt_templates` | Comma-separated pegRNA reverse transcriptase templates: `<RT template A>,<RT template B>`. Informs flap analysis and plotting. |

### Optional Arguments

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-r2`, `--fastq_r2` | Path to FASTQ R2 file for paired-end data. | `None` |
| `-o`, `--output_root` | Root output folder for CRISPResso2 and TwInsPEctor results. If not provided, a folder will be created in the current working directory based on the input FASTQ file names. | `auto` |
| `-rcm`, `--recoding_mode` | Run in recoding mode if the wild-type and twin prime edited sequences are the same length and should be evaluated as having only base substitutions. | `off` |
| `-ne`, `--min_num_base_edits` | Minimum number of base changes required for a read to be considered edited. | `3` (replacement)<br>`2` (recoding) |
| `-dmas`, `--default_min_aln_score` | Default minimum homology score for a read to align to the compound reference amplicon. | `30` |
| `-pfr`, `--plot_full_reads` | Display full read sequences in allele tables. | `off` |
| `-ncda`, `--no_collapse_displayed_alleles` | Do not combine alleles that become identical in the displayed allele table window. | `off` |
| `-naa`, `--no_alignment_adjustments` | Do not visually adjust homologies in the allele tables. This does not affect allele classification. | `off` |
| `-ied`, `--ignore_extraspacer_deletions` | Classification ignores deletions occurring beyond the spacers (outside edit window). | `off` |
| `-nf`, `--no_figures` | Skip all figures if only text outputs are desired. | `off` |
| `-nsf`, `--no_summary_figures` | Skip summary barplots if they are not desired. | `off` |
| `-nbf`, `--no_per_base_figures` | Skip per-base barplots if they are not desired. | `off` |
| `-nmf`, `--no_mutation_figures` | Skip mutation barplots if they are not desired. | `off` |
| `-pet`, `--plot_extended_tables` | Generates separate allele tables for each category. Controlled by `--max_n_rows` and `--min_frequency_alleles`. | `off` |
| `-pdf`, `--save_pdf` | Only save PDF versions of all plots. | `off` |
| `-npng`, `--no_save_png` | Do not save PNG versions of all plots. | `off` |
| `-mfa`, `--min_frequency_alleles` | Minimum percent read frequency required to report an allele in the allele tables. | `0.1` |
| `-mnr`, `--max_n_rows` | Maximum number of allele rows to display in the allele tables by category. | `25` |
| `-mna`, `--max_n_alleles_to_write` | Maximum number of alleles per category to write to the f7 text file. | `50` |
| `-nrr`, `--no_rerun` | Don't rerun CRISPResso2 if a run using the same parameters has already been finished. | `off` |
| `-kco`, `--keep_crispresso_outputs` | Don't delete CRISPResso2 output folders after analysis. | `off` |
| `--crispresso_args` | Additional arguments to pass to CRISPResso2 (wrapped in quotes); do not include `--n_processes`. | `""` |
| `-coa`, `--cleavage_offset_a` | Cleavage offset for pegRNA spacer A. | `-3` |
| `-cob`, `--cleavage_offset_b` | Cleavage offset for pegRNA spacer B. | `-3` |
| `-p`, `--n_processes`, `--n_threads` | Total process budget. CRISPResso2 divides it across reference runs; allele-table plotting uses up to `N` processes. | `1` |
| `-v`, `--verbose` | Print verbose CRISPResso2 output. | `off` |

### Usage

```
TwInsPEctor -r1 <FASTQ_R1> [-r2 <FASTQ_R2>] -w <WT_SEQUENCE> -t <TWINPE_SEQUENCE> -g <PEG_SPACER_A>,<PEG_SPACER_B> -rt RT_TEMPLATE_A>,<RT_TEMPLATE_B> [options]
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

### System Requirements
- **Operating Systems**: Linux, macOS
- **Tested on**: Rocky Linux 8.10
- **Hardware Requirements**: This software can run on a standard desktop computer and does not require any non-standard hardware.

## Demo

For instructions on running TwInsPEctor on the provided demo data, including expected outputs, please see [demo/demo.md](demo/demo.md).

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
