# TwInsPEctor Demo

For the demo, use this toy fastq file [here](https://github.com/clementlab/TwInsPEctor/blob/main/demo/demo.fastq.gz). You can download it with:

```
wget https://github.com/clementlab/TwInsPEctor/raw/refs/heads/main/demo/demo.fastq.gz
```

The wildtype reference amplicon is:

```
CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGCCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC
```

The twin-pe reference amplicon is:

```
CGCCGGGAACTGCCGCTGGCCCCCCACCGCCGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCCCTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC
```

The two guide spacer sequences are:

```
GCTGGCCCCCCACCGCCCCA,ACAGCGACTCCTGGAAGTGG 
```

The two RT templates are (optional):
```
TTGTACCGTACACCACTGAGACCGCGGTGGTTGACCAGACAAACC,GTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCC
```

Run the bioconda installation of TwInsPEctor using:

```
TwInsPEctor -r1 demo.fastq.gz -w CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGCCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC -t CGCCGGGAACTGCCGCTGGCCCCCCACCGCCGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCCCTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC -g GCTGGCCCCCCACCGCCCCA,ACAGCGACTCCTGGAAGTGG -rt TTGTACCGTACACCACTGAGACCGCGGTGGTTGACCAGACAAACC,GTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCC
```

Or the script using:

```
python TwInsPEctor.py -r1 demo.fastq.gz -w CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGCCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC -t CGCCGGGAACTGCCGCTGGCCCCCCACCGCCGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCCCTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTC -g GCTGGCCCCCCACCGCCCCA,ACAGCGACTCCTGGAAGTGG -rt TTGTACCGTACACCACTGAGACCGCGGTGGTTGACCAGACAAACC,GTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCCC
```

- modify path to TwInsPEctor.py and demo.fastq.gz as needed
- to specify output folder, include: --output_root path/to/folder/
- for shorter run times, include one of these: ``` --n_processes 4 ``` ``` --n_processes 8 ``` ``` --n_processes 12 ``` ``` --n_processes 16 ```

The following outputs are produced:
- TwInsPEctor/
  - Figures a1 to a9 - Summary bar plots detailing read inputs, stacked/categorized editing outcomes, and the most frequent alleles aligned to composite A/B, twinPE, and WT references.
    - For the demo, each bar in a1 should show 100 total reads, and each category in the categorized outcomes plots should show 10 reads.
  - b.base_plots/ (Figures b1 to b4) - Per-base bar plots analyzing contiguous editing, including 3' flap integration, 3' flap completion, and 5' flap removal.
  - c.mutation_plots/ (Figures c1 to c14) - Plots detailing insertion, deletion, and substitution positions and combinations, broken down by specific outcome categories (e.g., Perfect TPE, Dual Flap, Flap A, WT).
  - d.text_files/ (Files d1 to d10) - Raw text outputs containing base counts, mutation subtype counts, top alleles, reference sequences, and execution commands.
  - e.allele_tables/ (Figures e1 to e10) - (optional, use --plot_extended_tables) Allele tables for each outcome category, generating separate alignment plots against Composite A, Composite B, TwinPE, and WT references.
  - TwInsPEctor_report.html - The comprehensive interactive HTML report bundling all of the above plots and data.
