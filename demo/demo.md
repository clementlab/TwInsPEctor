# TwInsPEctor demo run

TwInsPEctor is a tool for analyzing twin prime editing outcomes from next-generation sequencing reads.

Each read is categorized into one of eight allele types with detailed visualizations:
- Perfect TPE:    complete programmed edit without indels.
- TPE Indel:      complete programmed edit with indels.
- Left Flap:      at least N consecutive programmed bases starting from the left but not from the right.
- Right Flap:     at least N consecutive programmed bases starting from the right but not from the left.
- Imperfect TPE:  incomplete programmed edit (neither or both flaps).
- Imperfect WT:   incomplete wildtype sequence and none of the programmed edit.
- WT Indel:       complete wildtype sequence with indels and none of the programmed edit.
- WT:             complete wildtype sequence without indels and none of the programmed edit.

# Processing overview
The tool requires four inputs:
- the input fastq file(s)
- the wildtype reference sequence
- the twin-pe edited sequence
- the two guide spacer sequences

# Demo overview
For this run, we'll be using the input fastq file [here](https://github.com/clementlab/TwInsPEctor/blob/main/demo/demo.fastq.gz). You can download it like this:
```
wget https://github.com/clementlab/TwInsPEctor/raw/refs/heads/main/demo/demo.fastq.gz
```
The wildtype reference sequence is:
```
CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGGAGTCCATTAGCAGAAGTGGCCCTTGGCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC
```
The twin-pe edited sequence is:
```
CGCCGGGAACTGCCGCTGGCCCCCCACCGCCAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC
```
The two guide spacer sequences are these:
```
GCTGGCCCCCCACCGCCCCA
```
```
GCGACTCCTGGAAGTGGCCA
```
We can run TwInsPEctor with the following command for bioconda installations:
```
TwInsPEctor -r1 demo.fastq.gz -w CGCCGGGAACTGCCGCTGGCCCCCCACCGCCCCAAGGATCTCCCGGTCCCCGCCCGGCGTGCTGACGTCACGGCGCTGCCCCAGGGTGTGCTGGGCAGGTCGCGGGGAGCGCTGGGAAATGGAGTCCATTAGCAGAAGTGGCCCTTGGCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC -t CGCCGGGAACTGCCGCTGGCCCCCCACCGCCAGGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACCTCCACTTCCAGGAGTCGCTGTGCCCCGATGCACACTGGGAAGTCCGCAGCTCCGAGGCGCCCAGTGGAAATCGCCAGATGAGGGCCTCCTC -g GCTGGCCCCCCACCGCCCCA,GCGACTCCTGGAAGTGGCCA
```
This produces the following outputs:
- TwInsPEctor_report.html
- TwInsPEctor_outputs
  - Figures a1 to a10 - barplots of the editing outcomes
  - Figures b1.a/b1.b to b8.a/b8.b - detailed allele tables
  - Figures c1 to c7 - raw outputs
